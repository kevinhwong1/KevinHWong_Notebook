#!/usr/bin/env python3
"""
05p_score_kmeans_sweep.py  (runs locally on the Mac; no Pegasus needed)

Scores the 05o KMeans sweep (ESM-1b vs ESMC-600M rescaled; K = 500 ... 10,000) against eggNOG, with the same
metrics as 05g / 05n, so the macrogene-only SATURN runs at K = 3,000 can be placed on the same curves.

Download first (from the Mac project root):
  rsync -avz kxw755@pegasus2.ccs.miami.edu:/scratch/dark_genes/SATURN_Mnemi/04_saturn_runs/20261009_kmeans_sweep/ \
        output/SATURN/20261002_macrogene_comparison_ESM1b_vs_ESMC600M/20261009_kmeans_sweep/

Inputs
  <sweep>/labels/<model>_K<K>_seed42.csv.gz (+ .json timing), <sweep>/validate.csv (optional)
  05h per-gene table (SATURN K3000 final hard labels: ESM1b_4sp, ESMC600M_resc_4sp)
  output/Eggnog/ALL4_gene_protein_metazoaOG_slim.csv, output/Eggnog/<sp>_eggnog_annotations.csv (Pfam)
Universe: genes in the sweep AND the 05h table. Primary OG per gene (as 05n; pooled recall agrees with 05g
to within 0.002). Null = labels shuffled within species.

Metrics per labelling (one row per model x K, plus the two SATURN K3000 runs):
  recall / precision / F1     cross-species ortholog pairs (same Metazoa OG) placed together (05g definitions)
  recall_1to1                 pooled recall for OGs with exactly one gene per species present
  recall_small_og             mean per-OG recall for OGs with 2-8 genes (robust to the few giant families)
  pct_mg_same_OG / _Pfam      SATURN Fig. 2d analogue (05n N4): % of macrogenes containing both species of a
                              pair that hold >= 1 cross-species pair with the same OG / a shared Pfam domain
                              (mean over the 6 species pairs; also Mlei pairs only)
  pct_genes_single_species    % genes in macrogenes with only one species (all genes; Drer separately)
  match / OGs merged / OG subdivided / cross-cutting   05n N2 gene categories (% genes in OGs with >= 2 genes)
Outputs (OUT = <sweep>/scores): P1_scores.csv, P0_kmeans_meta.csv, figures/FigP1..P2, summary.txt
"""
import argparse
import json
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--n_perm", type=int, default=3)
ap.add_argument("--reference", choices=["eggnog", "of_og", "of_hog"], default="eggnog",
                help="orthology reference: eggNOG Metazoa OGs (default) or OrthoFinder orthogroups / N0 HOGs (06c tables "
                     "copied to output/OrthoFinder/)")
ap.add_argument("--eggnog_genes_only", action="store_true",
                help="with an OrthoFinder reference: score only genes that also have an eggNOG Metazoa OG")
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--sweep_dir", default="", help="default: output/SATURN/<comparison>/20261009_kmeans_sweep")
args = ap.parse_args()

ROOT = Path(__file__).resolve().parents[2] / "output"
CMP = ROOT / "SATURN" / "20261002_macrogene_comparison_ESM1b_vs_ESMC600M"
SWEEP = Path(args.sweep_dir) if args.sweep_dir else CMP / "20261009_kmeans_sweep"
HPG = CMP / "20261006_norm_vs_macrogene" / "H_per_gene.csv"
EGGDIR = ROOT / "Eggnog"
OUT = SWEEP / ("scores" + ("" if args.reference == "eggnog" else "_" + args.reference) + ("_eggnoggenes" if args.eggnog_genes_only and args.reference != "eggnog" else ""))
(OUT / "figures").mkdir(parents=True, exist_ok=True)

SPECIES = ["Mlei", "Cgig", "Crob", "Drer"]
PAIRS = list(combinations(SPECIES, 2))
MODELS = {"ESM1b": ("ESM-1b", "#eb6834", "ESM1b_4sp"),
          "ESMC600M_rescaled": ("ESMC-600M rescaled", "#eda100", "ESMC600M_resc_4sp")}
CATS = ["match", "OGs merged", "OG subdivided", "cross-cutting"]
INK, MUTED = "#222222", "#8a8a85"
LOG = []
rng = np.random.default_rng(args.seed)


def log(m=""):
    print(m); LOG.append(str(m))



def load_reference(root, reference, eggnog_genes_only):
    """Orthology reference as (species, gene_id, protein_id, metazoa_OG) rows; 'metazoa_OG' holds the group id.
    eggnog = eggNOG Metazoa OGs (genes without one are unannotated); of_og / of_hog = OrthoFinder orthogroups /
    N0 HOGs from 06c (every gene has a group; genes with no ortholog are singletons).
    eggnog_genes_only: keep only genes that have an eggNOG Metazoa OG (same gene set for every reference)."""
    eg = pd.read_csv(root / "Eggnog" / "ALL4_gene_protein_metazoaOG_slim.csv").dropna(subset=["metazoa_OG"])
    if reference == "eggnog":
        return eg
    f = {"of_og": "ALL4_gene_orthofinder_OG.csv", "of_hog": "ALL4_gene_orthofinder_N0HOG.csv"}[reference]
    r = pd.read_csv(root / "OrthoFinder" / f).rename(columns={"group": "metazoa_OG"})
    if eggnog_genes_only:
        r = r.merge(eg[["species", "gene_id"]].drop_duplicates(), on=["species", "gene_id"])
    return r[["species", "gene_id", "protein_id", "metazoa_OG"]]


# ── labels ──
files = sorted((SWEEP / "labels").glob("*_K*_seed*.csv.gz"))
if not files:
    raise SystemExit("No sweep labels in " + str(SWEEP / "labels") + " (rsync them from Pegasus first)")
LAB, META = {}, []
for f in files:
    stem = f.name.replace(".csv.gz", "")
    model, rest = stem.rsplit("_K", 1)
    k = int(rest.split("_seed")[0])
    LAB[(model, k)] = pd.read_csv(f).set_index("key")["label"]
    js = f.with_name(stem + ".json")
    if js.exists():
        META.append(json.loads(js.read_text()))
h = pd.read_csv(HPG, usecols=["key", "species", "gene", "mg_ESM1b_4sp", "mg_ESMC600M_resc_4sp"]).set_index("key")
keys = sorted(set(h.index).intersection(*[set(v.index) for v in LAB.values()]))
g = h.loc[keys].reset_index()
N = len(g)
sp = g["species"].values
log("=" * 70)
log("05p  KMeans K sweep scored against eggNOG")
log("=" * 70)
log(f"Labelings: {len(LAB)}  " + str(sorted(LAB)))
log(f"Genes scored: {N}  " + str(g['species'].value_counts().reindex(SPECIES).to_dict()))
if META:
    M0 = pd.DataFrame(META).sort_values(["model", "k"])
    M0.to_csv(OUT / "P0_kmeans_meta.csv", index=False)
    log("KMeans runs:\n" + M0.to_string(index=False))
vf = SWEEP / "validate.csv"
if vf.exists():
    log("\nValidation (is KMeans a fair stand-in for SATURN's final macrogenes?):\n" + pd.read_csv(vf).round(3).to_string(index=False))

# ── eggNOG: primary OG per gene; Pfam tokens ──
egg = load_reference(ROOT, args.reference, args.eggnog_genes_only)
log("Orthology reference: " + args.reference + (" (eggNOG-annotated genes only)" if args.eggnog_genes_only else ""))
egg = egg[egg["species"].isin(SPECIES)]
cnt = egg.groupby(["species", "gene_id", "metazoa_OG"]).size().reset_index(name="n")
cnt = cnt.sort_values(["species", "gene_id", "n", "metazoa_OG"], ascending=[True, True, False, True])
prim = cnt.drop_duplicates(["species", "gene_id"]).rename(columns={"gene_id": "gene", "metazoa_OG": "og"})
g = g.merge(prim[["species", "gene", "og"]], on=["species", "gene"], how="left")
ann = g["og"].notna().values
ogc = np.full(N, -1)
ogc[ann] = pd.factorize(g.loc[ann, "og"])[0]
log(f"Genes with a Metazoa OG: {int(ann.sum())}")
ogm = egg.drop_duplicates(["species", "gene_id", "metazoa_OG"]).rename(columns={"gene_id": "gene", "metazoa_OG": "tok"})
pf = []
for s in SPECIES:
    e = pd.read_csv(EGGDIR / f"{s}_eggnog_annotations.csv", usecols=["gene_id", "PFAMs"])
    e = e.assign(tok=e["PFAMs"].astype(str).str.split(",")).explode("tok")
    e = e[~e["tok"].isin(["-", "nan"]) & e["tok"].notna()]
    pf.append(e[["gene_id", "tok"]].drop_duplicates().assign(species=s).rename(columns={"gene_id": "gene"}))
TOK = {"OG": g[["species", "gene"]].reset_index().merge(ogm[["species", "gene", "tok"]], on=["species", "gene"])[["index", "species", "tok"]],
       "Pfam": g[["species", "gene"]].reset_index().merge(pd.concat(pf), on=["species", "gene"])[["index", "species", "tok"]]}

# OG reference quantities (annotated genes)
A_og, A_sp = ogc[ann], sp[ann]
og_sp = pd.crosstab(A_og, A_sp).reindex(columns=SPECIES, fill_value=0).values.astype(float)


def xpairs(c):
    t = c.sum(axis=1)
    return (t ** 2 - (c ** 2).sum(axis=1)) / 2


og_total = xpairs(og_sp)
og_size = og_sp.sum(axis=1)
has = og_total > 0
one2one = has & (og_sp <= 1).all(axis=1)
small = has & (og_size <= 8)


def shuffled(lab):
    out = lab.copy()
    for s in SPECIES:
        m = sp == s
        out[m] = rng.permutation(out[m])
    return out


def score(lab, full=True):
    lab = pd.factorize(lab)[0]
    r = {}
    la = lab[ann]
    cell = pd.crosstab([A_og, la], A_sp).reindex(columns=SPECIES, fill_value=0)
    hits = pd.Series(xpairs(cell.values.astype(float)), index=cell.index.get_level_values(0)).groupby(level=0).sum()
    hits = hits.reindex(range(len(og_total)), fill_value=0).values
    mg_ann = pd.crosstab(la, A_sp).reindex(columns=SPECIES, fill_value=0).values.astype(float)
    r["recall"] = hits[has].sum() / og_total[has].sum()
    r["precision"] = hits.sum() / xpairs(mg_ann).sum()
    r["F1"] = 2 * r["recall"] * r["precision"] / (r["recall"] + r["precision"])
    r["recall_1to1"] = hits[one2one].sum() / og_total[one2one].sum()
    r["recall_small_og"] = float(np.mean(hits[small] / og_total[small]))
    # Fig. 2d analogue
    pres = pd.crosstab(lab, sp).reindex(columns=SPECIES, fill_value=0) > 0
    for tk, d in TOK.items():
        d = d.assign(mg=lab[d["index"].values])[["species", "mg", "tok"]].drop_duplicates()
        vals, mvals = [], []
        for s, t in PAIRS:
            both = pres.index[pres[s] & pres[t]]
            hit = d[d.species == s].merge(d[d.species == t], on=["mg", "tok"])["mg"].unique()
            v = 100 * np.isin(both, hit).mean() if len(both) else np.nan
            vals.append(v)
            if "Mlei" in (s, t):
                mvals.append(v)
        r["pct_mg_same_" + tk] = np.nanmean(vals)
        r["pct_mg_same_" + tk + "_Mlei_pairs"] = np.nanmean(mvals)
    if not full:
        return r
    cnt_all = pd.crosstab(lab, sp).reindex(columns=SPECIES, fill_value=0)
    single = (cnt_all > 0).sum(axis=1) == 1
    r["n_mg"] = len(cnt_all)
    r["median_mg_size"] = float(cnt_all.sum(axis=1).median())
    r["pct_mg_single_species"] = 100 * single.mean()
    r["pct_genes_single_species"] = 100 * cnt_all[single].values.sum() / N
    r["pct_Drer_genes_single_species"] = 100 * cnt_all.loc[single, "Drer"].sum() / cnt_all["Drer"].sum()
    r["pct_mg_all_4_species"] = 100 * ((cnt_all > 0).sum(axis=1) == 4).mean()
    # N2 categories
    keep = og_size[A_og] >= 2
    df = pd.DataFrame({"og": A_og, "mg": la})
    split = df["og"].map(df.groupby("og")["mg"].nunique()).values > 1
    mixed = df["mg"].map(df.groupby("mg")["og"].nunique()).values > 1
    cat = np.where(~split & ~mixed, CATS[0], np.where(~split & mixed, CATS[1], np.where(split & ~mixed, CATS[2], CATS[3])))
    r.update((pd.Series(cat[keep]).value_counts(normalize=True).reindex(CATS, fill_value=0) * 100).to_dict())
    return r


rows = []
for (model, k), s in sorted(LAB.items(), key=lambda x: (x[0][0], x[0][1])):
    lab = s.loc[g["key"]].values
    r = score(lab)
    nl = pd.DataFrame([score(shuffled(lab), full=False) for _ in range(args.n_perm)]).mean()
    rows.append({"model": model, "K": k, "source": "KMeans (05o)", **r, **{"null_" + c: v for c, v in nl.items()}})
    log(f"scored {model} K={k}")
for model, (_, _, col) in MODELS.items():
    lab = g["mg_" + col].values
    r = score(lab)
    rows.append({"model": model, "K": 3000, "source": "SATURN K3000 run", **r})
P = pd.DataFrame(rows)
P.to_csv(OUT / "P1_scores.csv", index=False)
show = ["model", "K", "source", "n_mg", "median_mg_size", "recall", "precision", "F1", "recall_1to1", "recall_small_og",
        "pct_mg_same_OG", "pct_mg_same_Pfam", "pct_mg_same_Pfam_Mlei_pairs", "pct_genes_single_species",
        "pct_Drer_genes_single_species"] + CATS
log("\n" + P[show].round(3).to_string(index=False))
log("\nShuffled-label null (recall, precision, % MG same OG / Pfam):\n" +
    P[P.source != "SATURN K3000 run"][["model", "K", "null_recall", "null_precision", "null_pct_mg_same_OG", "null_pct_mg_same_Pfam"]].round(4).to_string(index=False))

# ── Figures ──
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.titleweight": "bold", "axes.titlesize": 10})


def save(fig, name):
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / "figures" / f"{name}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


km = P[P.source == "KMeans (05o)"]
sat = P[P.source == "SATURN K3000 run"]
# FigP1: recall-precision trade-off across K
fig, ax = plt.subplots(figsize=(5.6, 4.4))
for m, (nm, col, _) in MODELS.items():
    d = km[km.model == m].sort_values("K")
    ax.plot(d["recall"], d["precision"], "-o", color=col, label=nm + " (KMeans)")
    for _, r in d.iterrows():
        ax.annotate(f"{r.K:,}", (r.recall, r.precision), textcoords="offset points", xytext=(5, 3), fontsize=7, color=col)
    s = sat[sat.model == m]
    ax.scatter(s["recall"], s["precision"], marker="*", s=160, color=col, edgecolor=INK, linewidth=0.6, zorder=5,
               label=nm + " SATURN K=3,000")
ax.set_xlabel("Recall (ortholog pairs placed together)"); ax.set_ylabel("Precision (co-grouped pairs that are orthologs)")
ax.set_title("Ortholog recall vs precision across K")
ax.legend(frameon=False, fontsize=8)
save(fig, "FigP1_recall_precision_by_K")

# FigP2: metrics vs K
panels = [("recall", "Recall"), ("precision", "Precision"), ("recall_1to1", "Recall, 1:1 OGs"),
          ("pct_mg_same_Pfam", "% macrogenes with a\nshared-Pfam cross-species pair"),
          ("pct_genes_single_species", "% genes in single-species\nmacrogenes (dashed: zebrafish)"),
          ("OGs merged", "% genes: OGs merged (solid)\nvs OG subdivided (dashed)")]
fig, axes = plt.subplots(2, 3, figsize=(11, 6.2))
for ax, (col, ttl) in zip(axes.ravel(), panels):
    for m, (nm, c, _) in MODELS.items():
        d = km[km.model == m].sort_values("K")
        ax.plot(d["K"], d[col], "-o", color=c, label=nm, ms=4)
        if "null_" + col in d and d["null_" + col].notna().any():
            ax.plot(d["K"], d["null_" + col], ":", color=c, lw=1)
        if col == "pct_genes_single_species":
            ax.plot(d["K"], d["pct_Drer_genes_single_species"], "--", color=c, lw=1.2)
        if col == "OGs merged":
            ax.plot(d["K"], d["OG subdivided"], "--", color=c, lw=1.2)
        s = sat[sat.model == m]
        ax.scatter(s["K"], s[col], marker="*", s=110, color=c, edgecolor=INK, linewidth=0.5, zorder=5)
    ax.set_xscale("log"); ax.set_xticks([500, 1000, 2000, 3000, 5000, 10000])
    ax.set_xticklabels(["500", "1k", "2k", "3k", "5k", "10k"]); ax.set_xlabel("K (number of macrogenes)")
    ax.set_title(ttl, fontsize=9)
axes[0, 0].legend(frameon=False, fontsize=8)
fig.suptitle("K sweep: dotted = shuffled labels; star = SATURN K=3,000 macrogene-only run", fontsize=9, color=MUTED)
save(fig, "FigP2_metrics_vs_K")

(OUT / "summary.txt").write_text("\n".join(LOG) + "\n")
log("\nWrote " + str(OUT))
