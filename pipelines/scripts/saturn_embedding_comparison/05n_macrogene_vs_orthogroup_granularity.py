#!/usr/bin/env python3
"""
05n_macrogene_vs_orthogroup_granularity.py  (runs locally on the Mac; no Pegasus needed)

How do SATURN macrogenes relate to eggNOG Metazoa orthologous groups (OGs)?
Are macrogenes coarser than OGs (OGs merged), finer (OGs subdivided), or do they cut across them?

Inputs (already downloaded):
  - 05h per-gene table: hard macrogene labels for ESM1b_4sp, ESMC600M_raw_4sp, ESMC600M_resc_4sp
  - canonical 5sp run: {SP}_gene_to_macrogene.csv (Spur dropped)
  - eggNOG slim table: output/Eggnog/ALL4_gene_protein_metazoaOG_slim.csv (species, gene_id, protein_id, metazoa_OG)
Universe: genes present in all four runs (same 31,702 genes as 05g/05j/05k).

One OG per gene: a gene's proteins can map to different OGs (isoforms). For the partition analyses (N1, N2)
each gene gets its "primary" OG = the OG carried by most of its proteins (ties: alphabetical). N3 also reports
the pooled recall with all OG memberships, which reproduces the 05g benchmark.

Steps
  N1  counts and sizes: number of OGs vs macrogenes, group-size distributions, effective number of groups
  N2  nesting: homogeneity / completeness / V-measure / ARI vs OGs (annotated genes only), per-OG split,
      per-macrogene purity, and a per-gene category:
        match          OG kept whole, macrogene contains only that OG
        OGs merged     OG kept whole, macrogene also contains other OGs   (macrogenes coarser)
        OG subdivided  OG split over macrogenes, macrogene contains only that OG  (macrogenes finer)
        cross-cutting  OG split AND macrogene mixes OGs
      (genes whose OG has >= 2 genes in the universe; null = macrogene labels shuffled within species)
  N3  cross-species recall by OG size: is the low pooled recall driven by large OGs being split?
      Also "1:1 OGs" (exactly one gene per species present) as the cleanest ortholog sets.
  N4  SATURN paper Fig. 2d analogue: % of macrogenes (containing both species of a pair) that hold >= 1
      cross-species homolog pair, homology = same Metazoa OG (strict) or shared Pfam domain (permissive).
      Paper: frog-zebrafish, 2,000 macrogenes, BLASTP homologs, 56% (top-1 gene/species) / 91.2% (top-10).

Outputs (OUT folder):
  N1_counts.csv, N1_size_distribution.csv
  N2_nesting_scores.csv, N2_gene_categories.csv, N2_og_split_by_size.csv, N2_mg_purity.csv
  N3_recall_by_og_size.csv, N3_recall_check.csv, N4_paper_fig2d_analogue.csv
  figures/FigN1..N3 (png + pdf), summary.txt
Usage:  python3 scripts/SATURN/05n_macrogene_vs_orthogroup_granularity.py [--n_perm 20] [--seed 42]
"""
import argparse
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--n_perm", type=int, default=20, help="shuffled-label null replicates")
ap.add_argument("--reference", choices=["eggnog", "of_og", "of_hog"], default="eggnog",
                help="orthology reference: eggNOG Metazoa OGs (default) or OrthoFinder orthogroups / N0 HOGs (06c tables "
                     "copied to output/OrthoFinder/)")
ap.add_argument("--eggnog_genes_only", action="store_true",
                help="with an OrthoFinder reference: score only genes that also have an eggNOG Metazoa OG")
ap.add_argument("--seed", type=int, default=42)
args = ap.parse_args()
REFNAME = {"eggnog": "eggNOG Metazoa OGs", "of_og": "OrthoFinder orthogroups", "of_hog": "OrthoFinder N0 HOGs"}[args.reference]

ROOT = Path(__file__).resolve().parents[2] / "output"
CMP = ROOT / "SATURN" / "20261002_macrogene_comparison_ESM1b_vs_ESMC600M"
HPG = CMP / "20261006_norm_vs_macrogene" / "H_per_gene.csv"
CANON = ROOT / "SATURN" / "20260711_SATURN_FINAL" / "macrogene"
EGG = ROOT / "Eggnog" / "ALL4_gene_protein_metazoaOG_slim.csv"
OUT = CMP / ("20261009_macrogene_vs_orthogroup" + ("" if args.reference == "eggnog" else "_" + args.reference) + ("_eggnoggenes" if args.eggnog_genes_only and args.reference != "eggnog" else ""))
(OUT / "figures").mkdir(parents=True, exist_ok=True)

SPECIES = ["Mlei", "Cgig", "Crob", "Drer"]
RUNS = ["canon5sp_ESM1b", "ESM1b_4sp", "ESMC600M_raw_4sp", "ESMC600M_resc_4sp"]
RLAB = {"canon5sp_ESM1b": "Canonical 5sp ESM-1b", "ESM1b_4sp": "ESM-1b 4sp",
        "ESMC600M_raw_4sp": "ESMC raw 4sp", "ESMC600M_resc_4sp": "ESMC rescaled 4sp"}
RCOL = {"canon5sp_ESM1b": "#2a78d6", "ESM1b_4sp": "#eb6834", "ESMC600M_raw_4sp": "#1baf7a", "ESMC600M_resc_4sp": "#eda100"}
INK, MUTED = "#222222", "#8a8a85"
CATS = ["match", "OGs merged", "OG subdivided", "cross-cutting"]
CCOL = {"match": "#256abf", "OGs merged": "#9ec5f4", "OG subdivided": "#eda100", "cross-cutting": "#c9c9c4"}
BINS = [(2, 2), (3, 4), (5, 8), (9, 16), (17, 32), (33, 10**9)]
BLAB = ["2", "3-4", "5-8", "9-16", "17-32", "33+"]
LOG = []


def log(m=""):
    print(m); LOG.append(str(m))


def binlab(n):
    for (lo, hi), l in zip(BINS, BLAB):
        if lo <= n <= hi:
            return l
    return "1"



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


# ── Load macrogene labels (same universe as 05j/05k) ──
h = pd.read_csv(HPG, usecols=["key", "species", "gene", "mg_ESM1b_4sp", "mg_ESMC600M_raw_4sp", "mg_ESMC600M_resc_4sp"])
can = pd.concat([pd.read_csv(CANON / f"{s}_gene_to_macrogene.csv")[["gene", "macrogene"]].assign(species=s) for s in SPECIES])
can["key"] = can["species"] + "_" + can["gene"]
g = h.merge(can[["key", "macrogene"]].rename(columns={"macrogene": "mg_canon5sp_ESM1b"}), on="key", how="inner")
g = g.reset_index(drop=True)
N = len(g)
log("=" * 70)
log("05n  Macrogenes vs eggNOG Metazoa orthologous groups")
log("=" * 70)
log("Genes in all four runs: " + str(N) + "  " + str(g["species"].value_counts().reindex(SPECIES).to_dict()))

# ── eggNOG: primary OG per gene + all OG memberships ──
egg = load_reference(ROOT, args.reference, args.eggnog_genes_only)
log("Orthology reference: " + args.reference + (" (eggNOG-annotated genes only)" if args.eggnog_genes_only else ""))
egg = egg[egg["species"].isin(SPECIES)]
cnt = egg.groupby(["species", "gene_id", "metazoa_OG"]).size().reset_index(name="n_prot")
cnt = cnt.sort_values(["species", "gene_id", "n_prot", "metazoa_OG"], ascending=[True, True, False, True])
n_og = cnt.groupby(["species", "gene_id"]).size().rename("n_og_memberships").reset_index()
prim = cnt.drop_duplicates(["species", "gene_id"])[["species", "gene_id", "metazoa_OG"]]
prim = prim.merge(n_og, on=["species", "gene_id"]).rename(columns={"gene_id": "gene", "metazoa_OG": "og"})
g = g.merge(prim, on=["species", "gene"], how="left")
ann = g["og"].notna().values
log("Genes with a Metazoa OG: " + str(int(ann.sum())) + " / " + str(N) + "  " +
    str({s: round(float(ann[g["species"].values == s].mean()), 3) for s in SPECIES}))
multi = int((g["n_og_memberships"] > 1).sum())
log("Genes whose proteins map to >1 OG (primary OG used in N1/N2): " + str(multi) +
    " (" + str(round(100 * multi / ann.sum(), 1)) + "% of annotated)")

A = g[ann].copy().reset_index(drop=True)          # annotated genes
spA = A["species"].values
og_codes, og_names = pd.factorize(A["og"])
og_size = np.bincount(og_codes)
og_nsp = A.groupby("og")["species"].nunique().reindex(og_names).values
A["og_size"] = og_size[og_codes]
A["og_nsp"] = og_nsp[og_codes]
rng = np.random.default_rng(args.seed)


def shuffled(labels, species):
    out = labels.copy()
    for s in SPECIES:
        m = species == s
        out[m] = rng.permutation(out[m])
    return out


# ════════════════════════════════════════════════════════════
# N1. Counts and sizes
# ════════════════════════════════════════════════════════════
log("\n── N1. How many groups, and how big? ──")


def ent(codes):
    p = np.bincount(codes) / len(codes); p = p[p > 0]
    return float(-(p * np.log(p)).sum())


rows, size_rows = [], []
og_tab = pd.Series(og_size)
rows.append({"grouping": REFNAME, "universe": "annotated genes", "n_genes": len(A),
             "n_groups": len(og_names), "n_groups_ge2_genes": int((og_tab >= 2).sum()),
             "n_singletons": int((og_tab == 1).sum()),
             "n_groups_multi_species": int((og_nsp >= 2).sum()),
             "n_groups_all_4_species": int((og_nsp == 4).sum()),
             "mean_size": og_tab.mean(), "median_size": og_tab.median(), "max_size": int(og_tab.max()),
             "effective_n_groups": np.exp(ent(og_codes))})
size_rows += [{"grouping": REFNAME, "size": int(k), "n_groups": int(v)} for k, v in og_tab.value_counts().items()]
for r in RUNS:
    for uni, df in [("all genes", g), ("annotated genes", A)]:
        lab = df["mg_" + r].values
        sz = pd.Series(lab).value_counts()
        nsp = df.groupby("mg_" + r)["species"].nunique()
        rows.append({"grouping": RLAB[r], "universe": uni, "n_genes": len(df), "n_groups": len(sz),
                     "n_groups_ge2_genes": int((sz >= 2).sum()), "n_singletons": int((sz == 1).sum()),
                     "n_groups_multi_species": int((nsp >= 2).sum()), "n_groups_all_4_species": int((nsp == 4).sum()),
                     "mean_size": sz.mean(), "median_size": sz.median(), "max_size": int(sz.max()),
                     "effective_n_groups": np.exp(ent(pd.factorize(lab)[0]))})
        if uni == "annotated genes":
            size_rows += [{"grouping": RLAB[r], "size": int(k), "n_groups": int(v)} for k, v in sz.value_counts().items()]
N1 = pd.DataFrame(rows)
N1.to_csv(OUT / "N1_counts.csv", index=False)
pd.DataFrame(size_rows).to_csv(OUT / "N1_size_distribution.csv", index=False)
log(N1.round(1).to_string(index=False))
og_by_nsp = pd.Series(og_nsp).value_counts().sort_index()
log("eggNOG OGs by number of species present: " + str(og_by_nsp.to_dict()))
log("Share of annotated genes in OGs of size >= 9: " +
    str(round(100 * (A["og_size"] >= 9).mean(), 1)) + "%;  in OGs of size 1: " + str(round(100 * (A["og_size"] == 1).mean(), 1)) + "%")

# ════════════════════════════════════════════════════════════
# N2. Nesting: do macrogenes merge OGs, split them, or cut across?
# ════════════════════════════════════════════════════════════
log("\n── N2. Nesting of OGs within macrogenes (annotated genes) ──")


def nesting_scores(c, k):
    """homogeneity/completeness/V-measure (Rosenberg & Hirschberg 2007) and ARI. c = reference (OG), k = macrogene."""
    k = pd.factorize(k)[0]
    n = len(c)
    joint = pd.factorize(c.astype(np.int64) * (k.max() + 1) + k)[0]
    Hc, Hk, Hj = ent(c), ent(k), ent(joint)
    hom = 1 - (Hj - Hk) / Hc if Hc > 0 else 1.0   # each macrogene holds one OG
    comp = 1 - (Hj - Hc) / Hk if Hk > 0 else 1.0  # each OG sits in one macrogene
    v = 2 * hom * comp / (hom + comp)
    c2 = lambda x: (x * (x - 1) / 2).sum()
    nij = np.bincount(joint).astype(float)
    sij, sa, sb = c2(nij), c2(np.bincount(c).astype(float)), c2(np.bincount(k).astype(float))
    exp = sa * sb / (n * (n - 1) / 2)
    ari = (sij - exp) / (0.5 * (sa + sb) - exp)
    return {"homogeneity": hom, "completeness": comp, "v_measure": v, "ARI": ari}


def gene_categories(c, k, keep):
    """Per-gene category for genes whose OG has >= 2 genes (keep mask)."""
    k = pd.factorize(k)[0]
    df = pd.DataFrame({"og": c, "mg": k})
    og_nmg = df.groupby("og")["mg"].nunique()
    mg_nog = df.groupby("mg")["og"].nunique()
    split = df["og"].map(og_nmg).values > 1
    mixed = df["mg"].map(mg_nog).values > 1
    cat = np.where(~split & ~mixed, "match", np.where(~split & mixed, "OGs merged",
                   np.where(split & ~mixed, "OG subdivided", "cross-cutting")))
    vc = pd.Series(cat[keep]).value_counts(normalize=True).reindex(CATS, fill_value=0) * 100
    return vc, split, mixed


keep = A["og_size"].values >= 2
log("Genes with OG size >= 2 (used for categories): " + str(int(keep.sum())) + " / " + str(len(A)))
sc_rows, cat_rows, split_rows, pur_rows = [], [], [], []
for r in RUNS:
    k = A["mg_" + r].values
    s = nesting_scores(og_codes, k)
    vc, split, mixed = gene_categories(og_codes, k, keep)
    nulls_s, nulls_c = [], []
    for _ in range(args.n_perm):
        kp = shuffled(k, spA)
        nulls_s.append(nesting_scores(og_codes, kp))
        nulls_c.append(gene_categories(og_codes, kp, keep)[0])
    ns, nc = pd.DataFrame(nulls_s).mean(), pd.concat(nulls_c, axis=1).mean(axis=1)
    sc_rows.append({"run": r, **s, **{"null_" + key: v for key, v in ns.items()}})
    cat_rows.append({"run": r, "kind": "observed", **vc.to_dict()})
    cat_rows.append({"run": r, "kind": "shuffled null", **nc.to_dict()})
    # per-OG: how many macrogenes does it span, share in its main macrogene; per OG size bin and species span
    df = pd.DataFrame({"og": og_codes, "mg": pd.factorize(k)[0], "species": spA})
    per_og = df.groupby("og").agg(size=("mg", "size"), n_mg=("mg", "nunique"),
                                  modal=("mg", lambda v: v.value_counts().iloc[0]), nsp=("species", "nunique"))
    per_og = per_og[per_og["size"] >= 2]
    per_og["frac_modal"] = per_og["modal"] / per_og["size"]
    per_og["bin"] = per_og["size"].map(binlab)
    for (b, multi_sp), d in per_og.groupby(["bin", per_og["nsp"] >= 2]):
        split_rows.append({"run": r, "og_size_bin": b, "og_scope": "multi-species" if multi_sp else "single-species",
                           "n_og": len(d), "pct_og_whole": 100 * (d["n_mg"] == 1).mean(),
                           "median_n_mg_spanned": d["n_mg"].median(), "mean_frac_in_main_mg": d["frac_modal"].mean()})
    # per-macrogene purity (macrogenes with >= 2 annotated genes)
    per_mg = df.groupby("mg").agg(n_ann=("og", "size"), n_og=("og", "nunique"),
                                  modal=("og", lambda v: v.value_counts().iloc[0]))
    per_mg = per_mg[per_mg["n_ann"] >= 2]
    pur_rows.append({"run": r, "n_mg_ge2_annotated": len(per_mg), "pct_mg_one_og": 100 * (per_mg["n_og"] == 1).mean(),
                     "median_n_og_per_mg": per_mg["n_og"].median(),
                     "mean_frac_main_og": (per_mg["modal"] / per_mg["n_ann"]).mean(),
                     "median_mg_size_annotated": per_mg["n_ann"].median()})
N2s, N2c = pd.DataFrame(sc_rows), pd.DataFrame(cat_rows)
N2o, N2p = pd.DataFrame(split_rows), pd.DataFrame(pur_rows)
N2s.to_csv(OUT / "N2_nesting_scores.csv", index=False)
N2c.to_csv(OUT / "N2_gene_categories.csv", index=False)
N2o.to_csv(OUT / "N2_og_split_by_size.csv", index=False)
N2p.to_csv(OUT / "N2_mg_purity.csv", index=False)
log("Nesting scores (homogeneity: macrogenes hold one OG; completeness: OGs stay in one macrogene):")
log(N2s.round(3).to_string(index=False))
log("\nGene categories (% of genes in OGs with >= 2 genes):")
log(N2c.round(1).to_string(index=False))
log("\nMacrogene purity (macrogenes with >= 2 annotated genes):")
log(N2p.round(2).to_string(index=False))
log("\n% of multi-species OGs kept whole, by OG size:")
d = N2o[N2o.og_scope == "multi-species"]
log(d.pivot(index="og_size_bin", columns="run", values="pct_og_whole").reindex(BLAB)[RUNS].round(1).to_string())
log("Number of multi-species OGs per size bin: " + str(d[d.run == RUNS[0]].set_index("og_size_bin")["n_og"].reindex(BLAB).to_dict()))

# ════════════════════════════════════════════════════════════
# N3. Cross-species recall by OG size
# ════════════════════════════════════════════════════════════
log("\n── N3. Cross-species ortholog-pair recall by OG size (primary OG) ──")


def xpairs(counts):
    """cross-species pairs from a (rows x species) count matrix: ((sum)^2 - sum of squares) / 2"""
    t = counts.sum(axis=1)
    return (t ** 2 - (counts ** 2).sum(axis=1)) / 2


og_sp = pd.crosstab(og_codes, spA).reindex(columns=SPECIES, fill_value=0)
og_total = pd.Series(xpairs(og_sp.values.astype(float)), index=og_sp.index)
one_to_one = (og_sp.values <= 1).all(axis=1) & (og_sp.values.sum(axis=1) >= 2)


def og_hits(k):
    df = pd.DataFrame({"og": og_codes, "mg": k, "species": spA})
    cell = df.groupby(["og", "mg", "species"]).size().unstack(fill_value=0).reindex(columns=SPECIES, fill_value=0)
    hp = pd.Series(xpairs(cell.values.astype(float)), index=cell.index.get_level_values(0))
    return hp.groupby(level=0).sum().reindex(og_total.index, fill_value=0)


has = og_total > 0
og_bin = pd.Series(og_size, index=og_sp.index).map(binlab)
groups = [(b, has & (og_bin == b)) for b in BLAB] + [("1:1 OGs", has & one_to_one), ("all", has)]
rec_rows = []
for r in RUNS:
    k = A["mg_" + r].values
    hits = og_hits(k)
    null_hits = [og_hits(shuffled(k, spA)) for _ in range(args.n_perm)]
    for b, m in groups:
        tot = og_total[m].sum()
        per = hits[m] / og_total[m]
        nh = np.mean([nh_[m].sum() for nh_ in null_hits]) / tot
        rec_rows.append({"run": r, "og_group": b, "n_og": int(m.sum()), "n_pairs": int(tot),
                         "pct_of_all_pairs": 100 * tot / og_total[has].sum(),
                         "recall_pooled": hits[m].sum() / tot, "recall_mean_per_og": per.mean(),
                         "pct_og_any_pair_together": 100 * (hits[m] > 0).mean(),
                         "pct_og_all_pairs_together": 100 * (per == 1).mean(), "null_recall_pooled": nh})
N3 = pd.DataFrame(rec_rows)
N3.to_csv(OUT / "N3_recall_by_og_size.csv", index=False)
for col, ttl in [("recall_pooled", "Pooled pair recall"), ("recall_mean_per_og", "Mean recall per OG (each OG weighted equally)"),
                 ("pct_og_any_pair_together", "% OGs with >= 1 cross-species pair together")]:
    log(ttl + ":")
    log(N3.pivot(index="og_group", columns="run", values=col).reindex(BLAB + ["1:1 OGs", "all"])[RUNS].round(3).to_string())
log("Share of all cross-species OG pairs contributed by each size bin (%):")
log(N3[N3.run == RUNS[0]].set_index("og_group")[["n_og", "n_pairs", "pct_of_all_pairs"]].reindex(BLAB + ["1:1 OGs", "all"]).round(1).to_string())

# check: pooled recall with ALL OG memberships (as in 05g) on the same genes
log("\nCheck vs 05g (all OG memberships per gene, unique cross-species gene pairs):")
idx = pd.Series(np.arange(N), index=pd.MultiIndex.from_arrays([g["species"], g["gene"]]))
e2 = egg.drop_duplicates(["species", "gene_id", "metazoa_OG"])
e2 = e2[pd.MultiIndex.from_arrays([e2["species"], e2["gene_id"]]).isin(idx.index)].copy()
e2["i"] = idx.loc[pd.MultiIndex.from_arrays([e2["species"], e2["gene_id"]])].values
spv = g["species"].values
true = []
for s, t in combinations(SPECIES, 2):
    m = e2[e2.species == s][["i", "metazoa_OG"]].merge(e2[e2.species == t][["i", "metazoa_OG"]], on="metazoa_OG")
    true.append(np.unique(m["i_x"].values.astype(np.int64) * N + m["i_y"].values))
true = np.unique(np.concatenate(true))
chk = []
for r in RUNS:
    lab = g["mg_" + r].values
    pred = []
    for s, t in combinations(SPECIES, 2):
        a = pd.DataFrame({"i": np.where(spv == s)[0]}); a["mg"] = lab[a["i"]]
        b = pd.DataFrame({"i": np.where(spv == t)[0]}); b["mg"] = lab[b["i"]]
        m = a.merge(b, on="mg")
        pred.append(m["i_x"].values.astype(np.int64) * N + m["i_y"].values)
    pred = np.unique(np.concatenate(pred))
    rec_all = np.intersect1d(pred, true, assume_unique=True).size / true.size
    rec_prim = N3[(N3.run == r) & (N3.og_group == "all")]["recall_pooled"].iloc[0]
    chk.append({"run": r, "recall_all_OG_memberships": rec_all, "recall_primary_OG": rec_prim})
CHK = pd.DataFrame(chk)
CHK.to_csv(OUT / "N3_recall_check.csv", index=False)
log(CHK.round(3).to_string(index=False))

# ════════════════════════════════════════════════════════════
# N4. SATURN paper Fig. 2d analogue: share of macrogenes containing a cross-species homolog pair
# ════════════════════════════════════════════════════════════
# Paper (Rosen et al. 2024, frog-zebrafish, 2,000 macrogenes, ESM2-15B): 56% of macrogenes have a BLASTP homolog
# pair among the top-1 gene per species, 91.2% among the top-10 (random: 0.25% / 18.8%).
# Here: hard assignment, ALL member genes (macrogenes hold ~10 genes, so close to the paper's top-10 version);
# homology = same eggNOG Metazoa OG (strict) or >= 1 shared Pfam domain (permissive, closer to a BLASTP hit).
# Denominator: macrogenes containing >= 1 gene from each species of the pair. Null: labels shuffled within species.
log("\n── N4. SATURN Fig. 2d analogue: % of macrogenes with a cross-species homolog pair ──")
EGGDIR = ROOT / "Eggnog"
pf = []
for s in SPECIES:
    e = pd.read_csv(EGGDIR / f"{s}_eggnog_annotations.csv", usecols=["gene_id", "PFAMs"])
    e = e.assign(pfam=e["PFAMs"].astype(str).str.split(",")).explode("pfam")
    e = e[(e["pfam"] != "-") & (e["pfam"] != "nan") & e["pfam"].notna()]
    pf.append(e[["gene_id", "pfam"]].drop_duplicates().assign(species=s))
pf = pd.concat(pf).rename(columns={"gene_id": "gene"})
ogm = egg.drop_duplicates(["species", "gene_id", "metazoa_OG"]).rename(columns={"gene_id": "gene", "metazoa_OG": "tok"})[["species", "gene", "tok"]]
pfm = pf.rename(columns={"pfam": "tok"})
TOK = {"same_OG": g[["species", "gene"]].merge(ogm, on=["species", "gene"]),
       "shared_Pfam": g[["species", "gene"]].merge(pfm, on=["species", "gene"])}


def homolog_mg(lab):
    """per species pair: n macrogenes with both species, % with >=1 cross-species pair sharing OG / Pfam"""
    gl = g[["species", "gene"]].assign(mg=lab)
    out = {}
    pres = pd.crosstab(gl["mg"], gl["species"]).reindex(columns=SPECIES, fill_value=0) > 0
    tk = {k: v.merge(gl, on=["species", "gene"])[["species", "mg", "tok"]].drop_duplicates() for k, v in TOK.items()}
    for s, t in combinations(SPECIES, 2):
        both = pres.index[pres[s] & pres[t]]
        row = {"n_mg_both_species": len(both)}
        for k, d in tk.items():
            hit = d[d.species == s].merge(d[d.species == t], on=["mg", "tok"])["mg"].unique()
            row["pct_" + k] = 100 * np.isin(both, hit).mean() if len(both) else np.nan
        out[(s, t)] = row
    return out


n4 = []
for r in RUNS:
    lab = g["mg_" + r].values
    obs = homolog_mg(lab)
    nul = [homolog_mg(shuffled(lab, g["species"].values)) for _ in range(max(5, args.n_perm // 4))]
    for (s, t), row in obs.items():
        n4.append({"run": r, "pair": s + "-" + t, **row,
                   "null_pct_same_OG": np.mean([x[(s, t)]["pct_same_OG"] for x in nul]),
                   "null_pct_shared_Pfam": np.mean([x[(s, t)]["pct_shared_Pfam"] for x in nul])})
N4 = pd.DataFrame(n4)
N4.to_csv(OUT / "N4_paper_fig2d_analogue.csv", index=False)
for col in ["pct_same_OG", "null_pct_same_OG", "pct_shared_Pfam", "null_pct_shared_Pfam"]:
    log(col + " (% of macrogenes containing both species):")
    log(N4.pivot(index="pair", columns="run", values=col)[RUNS].round(1).to_string())

# ── Figures ──
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.titleweight": "bold", "axes.titlesize": 10})


def save(fig, name):
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / "figures" / f"{name}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


# FigN1: group-size distributions (share of annotated genes in groups of each size)
SD = pd.DataFrame(size_rows)
fig, ax = plt.subplots(figsize=(6.4, 3.8))
edges = [1, 2, 3, 5, 9, 17, 33, 65, 10**6]
elab = ["1", "2", "3-4", "5-8", "9-16", "17-32", "33-64", "65+"]
series = [(REFNAME, "#222222")] + [(RLAB[r], RCOL[r]) for r in RUNS]
w = 0.16
for j, (nm, col) in enumerate(series):
    d = SD[SD.grouping == nm]
    genes = d["size"] * d["n_groups"]
    share = [100 * genes[(d["size"] >= lo) & (d["size"] < hi)].sum() / genes.sum() for lo, hi in zip(edges[:-1], edges[1:])]
    ax.bar(np.arange(len(elab)) + (j - 2) * w, share, w, color=col, label=nm)
ax.set_xticks(range(len(elab))); ax.set_xticklabels(elab)
ax.set_xlabel("Group size (annotated genes per OG or macrogene)"); ax.set_ylabel("% of annotated genes")
ax.set_title("Where genes sit: OG sizes vs macrogene sizes")
ax.legend(frameon=False, fontsize=8)
save(fig, "FigN1_group_sizes")

# FigN2: gene categories, observed vs shuffled
fig, ax = plt.subplots(figsize=(6.8, 3.4))
rowsN2 = []
for r in RUNS:
    rowsN2.append((RLAB[r], N2c[(N2c.run == r) & (N2c.kind == "observed")].iloc[0]))
rowsN2.append(("Shuffled (ESM-1b 4sp)", N2c[(N2c.run == "ESM1b_4sp") & (N2c.kind == "shuffled null")].iloc[0]))
for i, (nm, row) in enumerate(rowsN2[::-1]):
    left = 0
    for c in CATS:
        ax.barh(i, row[c], left=left, color=CCOL[c], label=c if i == 0 else None, edgecolor="white", linewidth=0.5)
        if row[c] >= 6:
            ax.text(left + row[c] / 2, i, f"{row[c]:.0f}", ha="center", va="center", fontsize=8,
                    color="white" if c in ("match",) else INK)
        left += row[c]
ax.set_yticks(range(len(rowsN2))); ax.set_yticklabels([n for n, _ in rowsN2[::-1]])
ax.set_xlabel("% of genes (OGs with ≥2 genes)"); ax.set_xlim(0, 100)
ax.set_title("How each gene's OG relates to its macrogene")
ax.legend(frameon=False, fontsize=8, ncol=4, loc="upper center", bbox_to_anchor=(0.45, -0.22))
save(fig, "FigN2_gene_categories")

# FigN3: recall by OG size
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.6), gridspec_kw={"width_ratios": [2, 1.2]})
x = np.arange(len(BLAB))
for r in RUNS:
    d = N3[N3.run == r].set_index("og_group").reindex(BLAB)
    a1.plot(x, d["recall_mean_per_og"], marker="o", color=RCOL[r], label=RLAB[r])
a1.set_xticks(x); a1.set_xticklabels(BLAB); a1.set_ylim(0, 1)
a1.set_xlabel("OG size (annotated genes in the 4 species)"); a1.set_ylabel("Mean cross-species recall per OG")
a1.set_title("Small OGs vs large OGs"); a1.legend(frameon=False, fontsize=8)
d = N3[N3.run == RUNS[0]].set_index("og_group").reindex(BLAB)
a2.bar(x, d["pct_of_all_pairs"], color="#9ec5f4")
a2.set_xticks(x); a2.set_xticklabels(BLAB)
a2.set_xlabel("OG size"); a2.set_ylabel("% of all OG pairs"); a2.set_title("Who drives pooled recall")
save(fig, "FigN3_recall_by_og_size")

(OUT / "summary.txt").write_text("\n".join(LOG) + "\n")
log("\nWrote outputs to " + str(OUT))
