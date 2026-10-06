#!/usr/bin/env python3
"""
05g_compare_runs_orthology.py

Compares SATURN macrogene runs side by side and benchmarks each against an
external reference that does not depend on any embedding model: eggNOG
Metazoa-level orthologous groups (OGs).

Default runs (any missing run is skipped with a warning):
  canon5sp_ESM1b     20260617 canonical 5-species run (Spur genes dropped)
  ESM1b_4sp          20261002 macrogene-only, ESM-1b
  ESMC600M_raw_4sp   20261002 macrogene-only, ESMC-600M as converted
  ESMC600M_resc_4sp  20261005 macrogene-only, ESMC-600M length-normalized (05e)

All comparisons use the genes common to every run (HV genes; Spur excluded).
Hard assignment = argmax of each gene's macrogene weights (as in 04p).

Sections
  R1 run descriptives  : macrogene sizes, singletons, species tiers, assignment
                         confidence, genes per tier per species
  R2 run agreement     : ARI / AMI between every pair of runs (pooled + per species)
  R3 ortholog benchmark: for each species pair, cross-species gene pairs that
                         share a Metazoa OG ("OG pairs") vs. cross-species gene
                         pairs placed in the same macrogene ("co-assigned"):
                           recall    = OG pairs co-assigned / OG pairs
                           precision = co-assigned pairs sharing an OG /
                                       co-assigned pairs where both genes have an OG
                           null      = same, after shuffling macrogene labels
                                       within each species (sizes preserved)
                           OG recovery = multi-species OGs with >=1 cross-species
                                       pair co-assigned / multi-species OGs
Outputs: --out_dir (CSV tables, PNG figures, summary.txt)
"""
import argparse
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("/scratch/dark_genes/SATURN_Mnemi")
RUNS = BASE / "04_saturn_runs"
DEFAULT_RUNS = ",".join([
    "canon5sp_ESM1b=" + str(RUNS / "20260617_Mlei_Cgig_Spur_Crob_Drer_K3000"),
    "ESM1b_4sp=" + str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESM1b_macrogene"),
    "ESMC600M_raw_4sp=" + str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESMC600M_macrogene"),
    "ESMC600M_resc_4sp=" + str(RUNS / "20261005_Mlei_Cgig_Crob_Drer_K3000_ESMC600M_rescaled_macrogene"),
])

ap = argparse.ArgumentParser()
ap.add_argument("--runs", default=DEFAULT_RUNS, help="comma-separated label=run_dir")
ap.add_argument("--eggnog", default=str(BASE / "01_proteomes/eggnog_merged/ALL4_gene_protein_metazoaOG_slim.csv"))
ap.add_argument("--species", default="Mlei,Cgig,Crob,Drer")
ap.add_argument("--n_perm", type=int, default=20)
ap.add_argument("--out_dir", default=str(RUNS / "20261005_run_comparison_orthology"))
args = ap.parse_args()

SPECIES = args.species.split(",")
NS = len(SPECIES)
PAIRS = [(SPECIES[i], SPECIES[j]) for i in range(NS) for j in range(i + 1, NS)]
OUT = Path(args.out_dir)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
RCOL_LIST = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # categorical slots 1-4, fixed order
INK, MUTED = "#222222", "#8a8a85"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SUMMARY = []


def log(msg=""):
    print(msg, flush=True)
    SUMMARY.append(str(msg))


def save(df, name, index=False):
    df.to_csv(OUT / name, index=index)
    print("  wrote " + name, flush=True)


def to_np(v):
    if hasattr(v, "detach"):
        v = v.detach().cpu().numpy()
    return np.asarray(v, dtype=np.float32).ravel()


def load_run(run_dir):
    hits = sorted((Path(run_dir) / "saturn_results").glob("*genes_to_macrogenes*.pkl"))
    if len(hits) != 1:
        return None, "expected 1 genes_to_macrogenes pkl, found " + str(len(hits))
    raw = pickle.load(open(hits[0], "rb"))
    keys = [k for k in raw if k.split("_", 1)[0] in SPECIES]
    W = np.vstack([to_np(raw[k]) for k in keys])
    top2 = np.partition(W, -2, axis=1)[:, -2:]
    df = pd.DataFrame({"key": keys, "species": [k.split("_", 1)[0] for k in keys],
                       "gene": [k.split("_", 1)[1] for k in keys], "mg": W.argmax(1),
                       "w_max": top2[:, 1], "w_margin": top2[:, 1] - top2[:, 0]})
    return (df, W.shape[1], str(hits[0])), None


t0 = time.time()
log("=" * 70)
log("SATURN run comparison + eggNOG ortholog benchmark")
log("=" * 70)
RUNSD = {}
for item in args.runs.split(","):
    lab, d = item.split("=", 1)
    res, err = load_run(d)
    if res is None:
        log("SKIP " + lab + ": " + err + " (" + d + ")")
        continue
    RUNSD[lab] = res
    log(lab + ": " + res[2])
LABS = list(RUNSD)
# colour follows the run, not its position: known labels keep fixed slots
KNOWN = ["canon5sp_ESM1b", "ESM1b_4sp", "ESMC600M_raw_4sp", "ESMC600M_resc_4sp"]
RCOL, extra = {}, iter(c for c in RCOL_LIST + ["#e87ba4", "#008300"])
for l in KNOWN:
    RCOL[l] = RCOL_LIST[KNOWN.index(l)]
for l in LABS:
    if l not in RCOL:
        RCOL[l] = next(c for c in ["#e87ba4", "#008300", "#4a3aa7", "#e34948"] if c not in RCOL.values())
if len(LABS) < 2:
    raise SystemExit("Need at least two runs")

common = sorted(set.intersection(*[set(RUNSD[l][0]["key"]) for l in LABS]))
cg = pd.DataFrame({"key": common, "species": [k.split("_", 1)[0] for k in common],
                   "gene": [k.split("_", 1)[1] for k in common]})
N = len(cg)
spv = cg["species"].values
LAB = {l: RUNSD[l][0].set_index("key").loc[common, "mg"].values for l in LABS}
log("Genes common to all runs: " + str(N) + "  " + str(cg["species"].value_counts().reindex(SPECIES).to_dict()))

# ════════════════════════════════════════════════════════════
# R1. Descriptives
# ════════════════════════════════════════════════════════════
log("\n── R1. Run descriptives (each run's own gene universe, Spur excluded) ──")
rows, tier_rows, gt_rows = [], [], []
SIZES = {}
for l in LABS:
    df, K, _ = RUNSD[l]
    cnt = pd.crosstab(df["mg"], df["species"]).reindex(index=range(K), columns=SPECIES, fill_value=0)
    size = cnt.sum(axis=1).values
    used = size > 0
    nsp = (cnt >= 1).sum(axis=1).values
    nsp3 = (cnt >= 3).sum(axis=1).values
    SIZES[l] = size[used]
    rows.append({"run": l, "n_genes": len(df), "n_macrogenes": K, "n_nonempty": int(used.sum()),
                 "size_median": float(np.median(size[used])), "size_mean": float(size[used].mean()),
                 "size_p90": float(np.percentile(size[used], 90)), "size_max": int(size.max()),
                 "n_singleton": int((size == 1).sum()), "n_size_le2": int(((size >= 1) & (size <= 2)).sum()),
                 "w_max_median": float(df["w_max"].median()), "w_margin_median": float(df["w_margin"].median())})
    for t in range(1, NS + 1):
        tier_rows.append({"run": l, "n_species": t, "n_mg_ge1gene": int(((nsp == t) & used).sum()),
                          "n_mg_ge3genes": int(((nsp3 == t) & used).sum())})
    tg = nsp[df["mg"].values]
    for s in SPECIES:
        m = df["species"].values == s
        gt_rows.append({"run": l, "species": s, **{"genes_in_" + str(t) + "sp_mg": int((tg[m] == t).sum()) for t in range(1, NS + 1)},
                        "frac_in_single_species_mg": float((tg[m] == 1).mean())})
r1 = pd.DataFrame(rows); save(r1, "R1_run_descriptives.csv"); log(r1.round(3).to_string(index=False))
tiers = pd.DataFrame(tier_rows); save(tiers, "R1_tier_distribution.csv")
log("\nMacrogenes by number of species (>=1 gene):")
log(tiers.pivot(index="n_species", columns="run", values="n_mg_ge1gene")[LABS].to_string())
gt = pd.DataFrame(gt_rows); save(gt, "R1_genes_per_tier_by_species.csv")
log("\nFraction of each species' genes in single-species macrogenes:")
log(gt.pivot(index="species", columns="run", values="frac_in_single_species_mg").loc[SPECIES, LABS].round(3).to_string())

# ════════════════════════════════════════════════════════════
# R2. Agreement between runs
# ════════════════════════════════════════════════════════════
log("\n── R2. Pairwise agreement between runs (common genes) ──")
rows = []
for scope in ["ALL"] + SPECIES:
    msk = np.ones(N, bool) if scope == "ALL" else spv == scope
    for i, a in enumerate(LABS):
        for b in LABS[i + 1:]:
            rows.append({"scope": scope, "run_a": a, "run_b": b,
                         "ARI": adjusted_rand_score(LAB[a][msk], LAB[b][msk]),
                         "AMI": adjusted_mutual_info_score(LAB[a][msk], LAB[b][msk])})
r2 = pd.DataFrame(rows); save(r2, "R2_pairwise_agreement.csv")
log(r2[r2["scope"] == "ALL"].round(3).to_string(index=False))

# ════════════════════════════════════════════════════════════
# R3. Ortholog benchmark
# ════════════════════════════════════════════════════════════
log("\n── R3. eggNOG Metazoa-OG benchmark ──")
egg = pd.read_csv(args.eggnog).dropna(subset=["metazoa_OG"])
egg = egg[egg["species"].isin(SPECIES)].drop_duplicates(["species", "gene_id", "metazoa_OG"])
idx = pd.Series(np.arange(N), index=pd.MultiIndex.from_arrays([cg["species"], cg["gene"]]))
mi = pd.MultiIndex.from_arrays([egg["species"], egg["gene_id"]])
egg = egg[mi.isin(idx.index)].copy()
egg["i"] = idx.loc[pd.MultiIndex.from_arrays([egg["species"], egg["gene_id"]])].values
annotated = np.zeros(N, bool); annotated[egg["i"].unique()] = True
log("Common genes with a Metazoa OG: " + str(int(annotated.sum())) + " / " + str(N) + "  " +
    str({s: round(float(annotated[spv == s].mean()), 3) for s in SPECIES}))


def pair_codes(left, right, on):
    m = left.merge(right, on=on)
    c = m["i_x"].values.astype(np.int64) * N + m["i_y"].values.astype(np.int64)
    return np.unique(c), m


# true OG pairs per species pair, and OG membership of each pair
TRUE, OGPAIRS = {}, {}
for s, t in PAIRS:
    a = egg[egg["species"] == s][["i", "metazoa_OG"]]
    b = egg[egg["species"] == t][["i", "metazoa_OG"]]
    codes, m = pair_codes(a, b, "metazoa_OG")
    TRUE[(s, t)] = codes
    OGPAIRS[(s, t)] = m.assign(code=m["i_x"].astype(np.int64) * N + m["i_y"].astype(np.int64))[["metazoa_OG", "code"]]
    log("  " + s + "-" + t + ": " + str(codes.size) + " OG pairs across " +
        str(m["metazoa_OG"].nunique()) + " shared OGs")


def bench(labels, s, t):
    ann = annotated
    ia = np.where((spv == s))[0]; ib = np.where((spv == t))[0]
    a = pd.DataFrame({"i": ia, "mg": labels[ia]}); b = pd.DataFrame({"i": ib, "mg": labels[ib]})
    pred, _ = pair_codes(a, b, "mg")
    true = TRUE[(s, t)]
    hit = np.intersect1d(pred, true, assume_unique=True)
    pi, pj = pred // N, pred % N
    pred_ann = pred[ann[pi] & ann[pj]]
    og = OGPAIRS[(s, t)]
    og_hit = og[og["code"].isin(hit)]["metazoa_OG"].nunique()
    return {"n_og_pairs": true.size, "n_coassigned_pairs": pred.size,
            "n_coassigned_both_annotated": pred_ann.size, "n_hits": hit.size,
            "recall": hit.size / max(true.size, 1),
            "precision": hit.size / max(pred_ann.size, 1),
            "og_recovery": og_hit / max(og["metazoa_OG"].nunique(), 1)}


rng = np.random.default_rng(42)
rows = []
for l in LABS:
    lab = LAB[l]
    for s, t in PAIRS:
        obs = bench(lab, s, t)
        nulls = []
        for _ in range(args.n_perm):
            perm = lab.copy()
            for sp in (s, t):
                m = spv == sp
                perm[m] = rng.permutation(perm[m])
            nulls.append(bench(perm, s, t))
        nd = pd.DataFrame(nulls)
        rows.append({"run": l, "species_1": s, "species_2": t, **obs,
                     "null_recall": nd["recall"].mean(), "null_precision": nd["precision"].mean(),
                     "null_og_recovery": nd["og_recovery"].mean(),
                     "recall_fold_over_null": obs["recall"] / max(nd["recall"].mean(), 1e-12),
                     "precision_fold_over_null": obs["precision"] / max(nd["precision"].mean(), 1e-12)})
    print("  benchmarked " + l, flush=True)
r3 = pd.DataFrame(rows)
r3["F1"] = 2 * r3["precision"] * r3["recall"] / (r3["precision"] + r3["recall"]).replace(0, np.nan)
save(r3, "R3_ortholog_benchmark_by_species_pair.csv")
show = ["run", "species_1", "species_2", "n_og_pairs", "n_coassigned_pairs", "recall", "null_recall",
        "precision", "null_precision", "og_recovery", "F1"]
log(r3[show].round(4).to_string(index=False))

# pooled over species pairs, and Mlei-involving pairs only
summ = []
for l in LABS:
    for scope, msk in [("all_pairs", np.ones(len(r3), bool)),
                       ("Mlei_pairs", (r3["species_1"] == "Mlei") | (r3["species_2"] == "Mlei")),
                       ("non_Mlei_pairs", (r3["species_1"] != "Mlei") & (r3["species_2"] != "Mlei"))]:
        d = r3[(r3["run"] == l) & msk]
        if d.empty:
            continue
        rec = d["n_hits"].sum() / d["n_og_pairs"].sum()
        prec = d["n_hits"].sum() / max(d["n_coassigned_both_annotated"].sum(), 1)
        summ.append({"run": l, "scope": scope, "recall_pooled": rec, "precision_pooled": prec,
                     "F1_pooled": 2 * rec * prec / (rec + prec) if rec + prec else np.nan,
                     "null_recall_mean": d["null_recall"].mean(), "null_precision_mean": d["null_precision"].mean(),
                     "og_recovery_mean": d["og_recovery"].mean(), "n_coassigned_pairs": int(d["n_coassigned_pairs"].sum())})
summ = pd.DataFrame(summ)
save(summ, "R3_ortholog_benchmark_summary.csv")
log("\nPooled ortholog benchmark:")
log(summ.round(4).to_string(index=False))
log("Recall rewards putting OG partners together; precision rewards not mixing unrelated genes.")
log("A run that co-assigns many more cross-species pairs gains recall and loses precision - read both.")

# ════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.titleweight": "bold", "axes.titlesize": 11})
FIG = OUT / "figures"
nl = len(LABS); w = 0.8 / nl

# FigR1: tier distribution
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
x = np.arange(1, NS + 1)
for ax, col, ttl in [(axes[0], "n_mg_ge1gene", "Species per macrogene (>=1 gene)"),
                     (axes[1], "n_mg_ge3genes", "Species per macrogene (>=3 genes)")]:
    for j, l in enumerate(LABS):
        d = tiers[tiers["run"] == l].set_index("n_species").loc[list(x), col]
        ax.bar(x - 0.4 + w / 2 + j * w, d, w * 0.95, color=RCOL[l], label=l)
    ax.set_xticks(x, [str(t) + "/" + str(NS) for t in x]); ax.set_title(ttl); ax.set_ylabel("Macrogenes")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "FigR1_tier_distribution.png", dpi=200); plt.close(fig)

# FigR2: macrogene size distribution (ECDF)
fig, ax = plt.subplots(figsize=(5.5, 3.8))
for l in LABS:
    v = np.sort(SIZES[l])
    ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", lw=2, color=RCOL[l], label=l)
ax.set_xscale("log"); ax.set_xlabel("Genes per macrogene (log)"); ax.set_ylabel("Cumulative fraction of macrogenes")
ax.set_title("Macrogene size distribution"); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "FigR2_macrogene_size_ecdf.png", dpi=200); plt.close(fig)

# FigR3: ARI heatmap (pooled)
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list("seq_blue", SEQ)
M = np.eye(nl)
for _, r in r2[r2["scope"] == "ALL"].iterrows():
    i, j = LABS.index(r["run_a"]), LABS.index(r["run_b"])
    M[i, j] = M[j, i] = r["ARI"]
fig, ax = plt.subplots(figsize=(4.8, 4.2))
im = ax.imshow(M, cmap=cmap, vmin=0, vmax=1)
for i in range(nl):
    for j in range(nl):
        ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=9,
                color="white" if M[i, j] > 0.55 else INK)
ax.set_xticks(range(nl), LABS, rotation=35, ha="right"); ax.set_yticks(range(nl), LABS)
ax.set_title("Macrogene agreement between runs (ARI)")
for sp in ax.spines.values():
    sp.set_visible(False)
fig.colorbar(im, ax=ax, fraction=0.046)
fig.tight_layout(); fig.savefig(FIG / "FigR3_run_agreement_ARI.png", dpi=200); plt.close(fig)

# FigR4: ortholog recall / precision / OG recovery per species pair
pl = [s + " - " + t for s, t in PAIRS]
fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharex=True)
x = np.arange(len(PAIRS))
for ax, met, ttl in [(axes[0], "recall", "Recall: OG pairs placed together"),
                     (axes[1], "precision", "Precision: co-assigned pairs that share an OG"),
                     (axes[2], "og_recovery", "OGs with >=1 cross-species pair together")]:
    for j, l in enumerate(LABS):
        d = r3[r3["run"] == l].set_index(["species_1", "species_2"]).loc[PAIRS]
        ax.bar(x - 0.4 + w / 2 + j * w, d[met], w * 0.95, color=RCOL[l], label=l)
    nul = r3[r3["run"] == LABS[0]].set_index(["species_1", "species_2"]).loc[PAIRS, "null_" + met]
    ax.scatter(x, nul, marker="_", s=500, color=INK, lw=2, label="shuffled null", zorder=3)
    ax.set_xticks(x, pl, rotation=30, ha="right"); ax.set_title(ttl, fontsize=10)
h, lb = axes[0].get_legend_handles_labels()
fig.legend(h, lb, loc="lower center", ncol=len(lb), frameon=False, fontsize=9)
fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(FIG / "FigR4_ortholog_benchmark.png", dpi=200); plt.close(fig)

# FigR5: precision vs recall, pooled, one point per run (all / Mlei pairs)
fig, ax = plt.subplots(figsize=(5.2, 4))
for l in LABS:
    for scope, mk in [("all_pairs", "o"), ("Mlei_pairs", "s")]:
        d = summ[(summ["run"] == l) & (summ["scope"] == scope)]
        if not d.empty:
            ax.scatter(d["recall_pooled"], d["precision_pooled"], s=70, marker=mk, color=RCOL[l],
                       edgecolor="white", lw=1.5, zorder=3)
handles = [plt.Line2D([], [], marker="o", ls="", color=RCOL[l], label=l) for l in LABS] + \
          [plt.Line2D([], [], marker="o", ls="", color=MUTED, label="all species pairs"),
           plt.Line2D([], [], marker="s", ls="", color=MUTED, label="Mlei pairs only")]
ax.legend(handles=handles, frameon=False, fontsize=7)
ax.set_xlabel("Recall (pooled)"); ax.set_ylabel("Precision (pooled)")
ax.set_title("Ortholog benchmark, pooled")
fig.tight_layout(); fig.savefig(FIG / "FigR5_precision_recall.png", dpi=200); plt.close(fig)

log("\nDone in " + str(round(time.time() - t0)) + " s. Outputs: " + str(OUT))
(OUT / "summary.txt").write_text("\n".join(SUMMARY) + "\n")
