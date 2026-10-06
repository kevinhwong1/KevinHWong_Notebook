#!/usr/bin/env python3
"""
05h_norm_vs_macrogene.py

Does ESMC embedding vector length explain how SATURN assigned genes to macrogenes?

x-axis everywhere: each gene's raw ESMC-600M vector length (L2 norm of the
FANTASIA embedding as converted by 05, before the 05e rescaling).

Runs compared (any missing run is skipped):
  ESM1b_4sp          control: grouping built without ESMC at all
  ESMC600M_raw_4sp   grouping built from the raw (length-varying) ESMC vectors
  ESMC600M_resc_4sp  control: same ESMC directions, all lengths equal (05e)
If length drove the raw-ESMC grouping, effects appear for the raw run only.

Panels / sections (common genes across runs; Spur excluded)
  H1  length decile vs macrogene outcome: share of genes in one-gene
      macrogenes, share in single-species macrogenes, median macrogene size.
      Deciles are computed WITHIN each species, so Drer's systematically
      shorter vectors do not masquerade as a length effect.
  H2  share of variance in log10(length) explained by macrogene membership
      (eta^2 = between-macrogene SS / total SS), vs a null where macrogene
      labels are shuffled across genes (sizes kept). Pooled and per species.
  H3  spread of member lengths within each macrogene (SD of log10 length,
      macrogenes with >= 3 genes), vs shuffled null.
  H4  gene stability vs length: share of a gene's macrogene partners that are
      the same in the ESM-1b run and each ESMC run, by within-species decile.
Outputs: --out_dir (CSV tables, PNG figures, summary.txt)
"""
import argparse
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("/scratch/dark_genes/SATURN_Mnemi")
RUNS = BASE / "04_saturn_runs"
DEFAULT_RUNS = ",".join([
    "ESM1b_4sp=" + str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESM1b_macrogene"),
    "ESMC600M_raw_4sp=" + str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESMC600M_macrogene"),
    "ESMC600M_resc_4sp=" + str(RUNS / "20261005_Mlei_Cgig_Crob_Drer_K3000_ESMC600M_rescaled_macrogene"),
])

ap = argparse.ArgumentParser()
ap.add_argument("--runs", default=DEFAULT_RUNS, help="comma-separated label=run_dir; first run = reference for H4")
ap.add_argument("--esmc_raw_dir", default=str(BASE / "03_embeddings_ESMC600M"))
ap.add_argument("--species", default="Mlei,Cgig,Crob,Drer")
ap.add_argument("--n_perm", type=int, default=100)
ap.add_argument("--out_dir", default=str(RUNS / "20261006_norm_vs_macrogene"))
args = ap.parse_args()

SPECIES = args.species.split(",")
FOLDER = {"Mlei": "Mnemi"}
OUT = Path(args.out_dir)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
RCOL = {"canon5sp_ESM1b": "#2a78d6", "ESM1b_4sp": "#eb6834",
        "ESMC600M_raw_4sp": "#1baf7a", "ESMC600M_resc_4sp": "#eda100"}
INK, MUTED, NULLC = "#222222", "#8a8a85", "#b4b2a9"
SUMMARY = []


def log(msg=""):
    print(msg, flush=True)
    SUMMARY.append(str(msg))


def save(df, name, index=False):
    df.to_csv(OUT / name, index=index)
    print("  wrote " + name, flush=True)


def vec(v):
    if hasattr(v, "detach"):
        v = v.detach().cpu().numpy()
    return np.asarray(v, dtype=np.float32).ravel()


def load_pt(p):
    import torch
    return torch.load(p, map_location="cpu")


def load_labels(run_dir):
    hits = sorted((Path(run_dir) / "saturn_results").glob("*genes_to_macrogenes*.pkl"))
    if len(hits) != 1:
        return None
    raw = pickle.load(open(hits[0], "rb"))
    keys = [k for k in raw if k.split("_", 1)[0] in SPECIES]
    return pd.Series([int(np.argmax(vec(raw[k]))) for k in keys], index=keys)


t0 = time.time()
log("=" * 70)
log("ESMC vector length vs macrogene assignment")
log("=" * 70)
LAB = {}
for item in args.runs.split(","):
    lab, d = item.split("=", 1)
    s = load_labels(d)
    if s is None:
        log("SKIP " + lab + " (no unique genes_to_macrogenes pkl in " + d + ")")
        continue
    LAB[lab] = s
    log(lab + ": " + d)
RUNL = list(LAB)
REF = RUNL[0]
for l in RUNL:
    RCOL.setdefault(l, MUTED)

common = sorted(set.intersection(*[set(s.index) for s in LAB.values()]))
g = pd.DataFrame({"key": common, "species": [k.split("_", 1)[0] for k in common],
                  "gene": [k.split("_", 1)[1] for k in common]})

# ── raw ESMC vector length per gene (exact, then case-insensitive key match) ──
norm = np.full(len(g), np.nan)
for s in SPECIES:
    f = FOLDER.get(s, s)
    d = load_pt(Path(args.esmc_raw_dir) / f / (f + "_gene_embeddings.pt"))
    lower = {}
    for k in d:
        lower.setdefault(k.lower(), []).append(k)
    for i in np.where(g["species"].values == s)[0]:
        gn = g.at[i, "gene"]
        k = gn if gn in d else (lower[gn.lower()][0] if len(lower.get(gn.lower(), [])) == 1 else None)
        if k is not None:
            norm[i] = np.linalg.norm(vec(d[k]))
    del d
ok = ~np.isnan(norm)
if (~ok).sum():
    log("Genes without an ESMC vector (dropped): " + str(int((~ok).sum())))
g = g[ok].reset_index(drop=True)
norm = norm[ok]
g["esmc_norm"] = norm
g["log10_norm"] = np.log10(norm)
# within-species deciles (1 = shortest)
g["decile"] = g.groupby("species")["esmc_norm"].transform(lambda v: pd.qcut(v.rank(method="first"), 10, labels=False) + 1)
spv = g["species"].values
N = len(g)
L = {l: LAB[l].loc[g["key"]].values for l in RUNL}
log("Genes analysed: " + str(N) + "  " + str(g["species"].value_counts().reindex(SPECIES).to_dict()))
log("ESMC length by species (median): " + str(g.groupby("species")["esmc_norm"].median().reindex(SPECIES).round(0).to_dict()))

# per-gene macrogene properties in each run (on the common gene set)
for l in RUNL:
    lab = L[l]
    size = pd.Series(lab).map(pd.Series(lab).value_counts()).values
    nsp = pd.DataFrame({"mg": lab, "sp": spv}).groupby("mg")["sp"].nunique()
    g["mg_" + l] = lab
    g["mg_size_" + l] = size
    g["single_gene_" + l] = size == 1
    g["single_species_" + l] = pd.Series(lab).map(nsp).values == 1

# ════════════════════════════════════════════════════════════
# H1. Length decile vs macrogene outcome
# ════════════════════════════════════════════════════════════
log("\n── H1. Within-species ESMC length decile vs macrogene outcome ──")
rows = []
for l in RUNL:
    for dec, d in g.groupby("decile"):
        rows.append({"run": l, "decile": int(dec), "n_genes": len(d),
                     "frac_single_gene_mg": d["single_gene_" + l].mean(),
                     "frac_single_species_mg": d["single_species_" + l].mean(),
                     "median_mg_size": d["mg_size_" + l].median()})
h1 = pd.DataFrame(rows)
save(h1, "H1_decile_vs_macrogene_outcome.csv")
log(h1.pivot(index="decile", columns="run", values="frac_single_gene_mg")[RUNL].round(3).to_string())
rows = []
for l in RUNL:
    for (s, dec), d in g.groupby(["species", "decile"]):
        rows.append({"run": l, "species": s, "decile": int(dec),
                     "frac_single_gene_mg": d["single_gene_" + l].mean(),
                     "frac_single_species_mg": d["single_species_" + l].mean(),
                     "median_mg_size": d["mg_size_" + l].median()})
save(pd.DataFrame(rows), "H1_decile_vs_macrogene_outcome_by_species.csv")
# simple trend summary: Spearman of decile vs outcome per run
tr = []
for l in RUNL:
    for col in ["single_gene_", "single_species_", "mg_size_"]:
        tr.append({"run": l, "outcome": col.strip("_"),
                   "spearman_rho_vs_decile": spearmanr(g["decile"], g[col + l].astype(float)).correlation})
tr = pd.DataFrame(tr)
save(tr, "H1_trend_spearman.csv")
log("\nSpearman rho, within-species length decile vs outcome (per gene):")
log(tr.pivot(index="outcome", columns="run", values="spearman_rho_vs_decile")[RUNL].round(3).to_string())

# ════════════════════════════════════════════════════════════
# H2. Variance in log length explained by macrogene membership
# ════════════════════════════════════════════════════════════
log("\n── H2. eta^2: share of variance in log10(ESMC length) explained by macrogene ──")


def eta2(y, lab):
    tot = ((y - y.mean()) ** 2).sum()
    df = pd.DataFrame({"y": y, "m": lab})
    grp = df.groupby("m")["y"].agg(["mean", "size"])
    btw = (grp["size"] * (grp["mean"] - y.mean()) ** 2).sum()
    return btw / tot


rng = np.random.default_rng(42)
rows = []
for scope in ["ALL"] + SPECIES:
    m = np.ones(N, bool) if scope == "ALL" else spv == scope
    y = g["log10_norm"].values[m]
    for l in RUNL:
        lab = L[l][m]
        obs = eta2(y, lab)
        null = np.array([eta2(y, rng.permutation(lab)) for _ in range(args.n_perm)])
        rows.append({"scope": scope, "run": l, "eta2": obs, "null_mean": null.mean(),
                     "null_q975": np.quantile(null, 0.975), "eta2_minus_null": obs - null.mean(),
                     "n_genes": int(m.sum()), "n_macrogenes": len(np.unique(lab))})
h2 = pd.DataFrame(rows)
save(h2, "H2_eta2_length_by_macrogene.csv")
log(h2.round(3).to_string(index=False))
log("Many small macrogenes inflate eta^2 even by chance -> compare with null_mean, not with 0.")

# ════════════════════════════════════════════════════════════
# H3. Within-macrogene spread of length
# ════════════════════════════════════════════════════════════
log("\n── H3. SD of log10(ESMC length) within macrogenes (>= 3 genes) ──")


def within_sd(y, lab, min_n=3):
    df = pd.DataFrame({"y": y, "m": lab})
    s = df.groupby("m")["y"].agg(["std", "size"])
    return s.loc[s["size"] >= min_n, "std"].values


y_all = g["log10_norm"].values
SDS, NULLSD, rows = {}, {}, []
for l in RUNL:
    SDS[l] = within_sd(y_all, L[l])
    NULLSD[l] = within_sd(y_all, rng.permutation(L[l]))
    rows.append({"run": l, "n_mg_ge3": len(SDS[l]), "median_within_sd": np.median(SDS[l]),
                 "null_median_within_sd": np.median(NULLSD[l]),
                 "ratio_to_null": np.median(SDS[l]) / np.median(NULLSD[l])})
h3 = pd.DataFrame(rows)
save(h3, "H3_within_macrogene_length_spread.csv")
log(h3.round(4).to_string(index=False))
log("ratio_to_null well below 1 => macrogene members have unusually similar lengths")

# ════════════════════════════════════════════════════════════
# H4. Gene stability vs length
# ════════════════════════════════════════════════════════════
log("\n── H4. Co-member stability vs " + REF + ", by within-species ESMC length decile ──")


def comember_jaccard(la, lb):
    Ka, Kb = la.max() + 1, lb.max() + 1
    C = sparse.coo_matrix((np.ones(len(la)), (la, lb)), shape=(Ka, Kb)).tocsr()
    c = np.asarray(C[la, lb]).ravel()
    sa = np.bincount(la, minlength=Ka)[la]
    sb = np.bincount(lb, minlength=Kb)[lb]
    union = sa + sb - c - 1
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(union > 0, (c - 1) / union, np.nan)


rows = []
CMP = [l for l in RUNL if l != REF]
for l in CMP:
    j = comember_jaccard(L[REF], L[l])
    g["jaccard_" + REF + "_vs_" + l] = j
    for dec, idx in g.groupby("decile").groups.items():
        v = j[np.asarray(list(idx))]
        v = v[~np.isnan(v)]
        rows.append({"comparison": REF + " vs " + l, "run": l, "decile": int(dec), "n": len(v),
                     "mean_jaccard": v.mean(), "se": v.std(ddof=1) / np.sqrt(len(v)),
                     "frac_no_shared_partners": (v == 0).mean()})
h4 = pd.DataFrame(rows)
save(h4, "H4_stability_vs_length_decile.csv")
if len(h4):
    log(h4.pivot(index="decile", columns="run", values="mean_jaccard").round(3).to_string())
save(g, "H_per_gene.csv")

# ════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.titleweight": "bold", "axes.titlesize": 10})
FIG = OUT / "figures"
dec = np.arange(1, 11)

# FigH1
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
for ax, col, ttl, yl in [(axes[0], "frac_single_gene_mg", "Genes alone in a macrogene", "Fraction of genes"),
                         (axes[1], "frac_single_species_mg", "Genes in single-species macrogenes", "Fraction of genes"),
                         (axes[2], "median_mg_size", "Size of the gene's macrogene", "Median genes per macrogene")]:
    for l in RUNL:
        d = h1[h1["run"] == l].set_index("decile").loc[dec]
        ax.plot(dec, d[col], "-o", lw=2, ms=5, color=RCOL[l], label=l)
    ax.set_xticks(dec); ax.set_xlabel("ESMC vector length decile (within species; 1 = shortest)")
    ax.set_title(ttl); ax.set_ylabel(yl)
h, lb = axes[0].get_legend_handles_labels()
fig.legend(h, lb, loc="lower center", ncol=len(lb), frameon=False)
fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(FIG / "FigH1_length_decile_vs_macrogene.png", dpi=200); plt.close(fig)

# FigH2
scopes = ["ALL"] + SPECIES
fig, ax = plt.subplots(figsize=(7.5, 3.8))
x = np.arange(len(scopes)); w = 0.8 / len(RUNL)
for j, l in enumerate(RUNL):
    d = h2[h2["run"] == l].set_index("scope").loc[scopes]
    xs = x - 0.4 + w / 2 + j * w
    ax.bar(xs, d["eta2"], w * 0.95, color=RCOL[l], label=l)
    ax.scatter(xs, d["null_mean"], marker="_", s=150, color=INK, lw=2, zorder=3,
               label="shuffled null" if j == 0 else None)
ax.set_xticks(x, ["All species"] + SPECIES); ax.set_ylim(0, 1)
ax.set_ylabel("Variance in log10 length\nexplained by macrogene (η²)")
ax.set_title("How much of ESMC vector length the macrogenes explain")
ax.legend(frameon=False, fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig(FIG / "FigH2_eta2_length_by_macrogene.png", dpi=200); plt.close(fig)

# FigH3
fig, ax = plt.subplots(figsize=(5.8, 3.8))
for l in RUNL:
    v = np.sort(SDS[l]); ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", lw=2, color=RCOL[l], label=l)
v = np.sort(np.concatenate([NULLSD[l] for l in RUNL]))
ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", lw=2, ls="--", color=NULLC, label="shuffled null")
ax.set_xlabel("SD of log10 ESMC length within a macrogene"); ax.set_ylabel("Cumulative fraction of macrogenes")
ax.set_title("Spread of vector length inside macrogenes (≥3 genes)")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "FigH3_within_macrogene_length_spread.png", dpi=200); plt.close(fig)

# FigH4
if len(h4):
    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    for l in CMP:
        d = h4[h4["run"] == l].set_index("decile").loc[dec]
        ax.errorbar(dec, d["mean_jaccard"], yerr=1.96 * d["se"], fmt="-o", lw=2, ms=5, capsize=0,
                    color=RCOL[l], label=REF + " vs " + l)
    ax.set_xticks(dec); ax.set_xlabel("ESMC vector length decile (within species; 1 = shortest)")
    ax.set_ylabel("Mean share of macrogene partners kept")
    ax.set_title("Gene stability vs vector length")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "FigH4_stability_vs_length.png", dpi=200); plt.close(fig)

log("\nDone in " + str(round(time.time() - t0)) + " s. Outputs: " + str(OUT))
(OUT / "summary.txt").write_text("\n".join(SUMMARY) + "\n")
