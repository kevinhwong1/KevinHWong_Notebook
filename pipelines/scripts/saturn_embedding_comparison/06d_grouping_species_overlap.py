#!/usr/bin/env python3
"""
06d_grouping_species_overlap.py  (runs locally on the Mac; no Pegasus needed)

Puts the two orthology references and the four macrogene sets on the same footing:
  eggNOG Metazoa OGs | OrthoFinder orthogroups | canonical 5sp ESM-1b | ESM-1b 4sp | ESMC raw 4sp | ESMC rescaled 4sp
and asks, for each grouping: which species end up together, and how do orthogroups split / merge with macrogenes?

Universe: the 31,702 SATURN genes present in all four macrogene runs (same as 05j / 05n / 05p).
Orthogroups were inferred from whole proteomes; only their labels on these genes are used here, so group sizes
and species counts are counted WITHIN the SATURN gene set. eggNOG: primary OG per gene (as 05n); genes without
an eggNOG OG are "unassigned". Two universes: all SATURN genes (eggNOG reported on its assigned genes) and
eggNOG-annotated genes only (all six groupings defined on the same 21,674 genes).

Outputs (OUT):
  D0_orthofinder_per_species.csv     OrthoFinder Statistics_PerSpecies (first block), if the stats folder was copied
  D0_orthofinder_orthologues.csv     ortholog counts per species pair and type (1:1, 1:many, many:1, many:many)
  D1_species_span_summary.csv        per grouping x universe: groups, sizes, genes / groups by no. of species spanned
  D2_species_combinations.csv        per grouping x universe: groups and genes per exact species combination
  D3_species_pair_cooccurrence.csv   per grouping x universe: % of species A's genes whose group contains species B
  D3b_per_species_multi.csv          per grouping x universe x species: % of genes in groups with >= 1 other species
  D4_pair_agreement.csv              per pair of groupings: cross-species gene pairs shared (containment, Jaccard)
  D5_nesting_vs_reference.csv        each grouping vs each orthology reference: homogeneity/completeness/ARI and
                                     % of genes match / OGs merged / OG subdivided / cross-cutting (as 05n)
  figures/  FigD1..FigD5 (.png for the notebook, .pdf)
  summary.txt

Optional input for D0 (from Pegasus):
  rsync -av kxw755@pegasus2.ccs.miami.edu:/scratch/dark_genes/SATURN_Mnemi/06_orthofinder/Results_20261015/Results_Oct15/Comparative_Genomics_Statistics \
      <repo>/output/OrthoFinder/
"""
import argparse
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2] / "output"
ap = argparse.ArgumentParser()
ap.add_argument("--out_dir", default="", help="default: output/SATURN/<comparison>/20261015_grouping_species_overlap")
args = ap.parse_args()

CMP = ROOT / "SATURN" / "20261002_macrogene_comparison_ESM1b_vs_ESMC600M"
HPG = CMP / "20261006_norm_vs_macrogene" / "H_per_gene.csv"
CANON = ROOT / "SATURN" / "20260711_SATURN_FINAL" / "macrogene"
EGG = ROOT / "Eggnog" / "ALL4_gene_protein_metazoaOG_slim.csv"
OFDIR = ROOT / "OrthoFinder"
OUT = Path(args.out_dir) if args.out_dir else CMP / "20261015_grouping_species_overlap"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SPECIES = ["Mlei", "Cgig", "Crob", "Drer"]
GROUPINGS = ["eggNOG", "OrthoFinder", "canon5sp_ESM1b", "ESM1b_4sp", "ESMC600M_raw_4sp", "ESMC600M_resc_4sp"]
LAB = {"eggNOG": "eggNOG Metazoa OGs", "OrthoFinder": "OrthoFinder orthogroups",
       "canon5sp_ESM1b": "Canonical 5sp ESM-1b", "ESM1b_4sp": "ESM-1b 4sp",
       "ESMC600M_raw_4sp": "ESMC raw 4sp", "ESMC600M_resc_4sp": "ESMC rescaled 4sp"}
COL = {"eggNOG": "#222222", "OrthoFinder": "#7b4fc9", "canon5sp_ESM1b": "#2a78d6", "ESM1b_4sp": "#eb6834",
       "ESMC600M_raw_4sp": "#1baf7a", "ESMC600M_resc_4sp": "#eda100"}
REFS = ["eggNOG", "OrthoFinder"]
CATS = ["match", "OGs merged", "OG subdivided", "cross-cutting"]
CATCOL = {"match": "#184f95", "OGs merged": "#6da7ec", "OG subdivided": "#eda100", "cross-cutting": "#d9d9d4"}
NSPCOL = {1: "#d9d9d4", 2: "#9ec5f4", 3: "#3987e5", 4: "#0d366b"}
INK, MUTED = "#222222", "#8a8a85"
COMBOS = [c for k in range(1, 5) for c in combinations(SPECIES, k)]
LOG = []


def log(m=""):
    print(m); LOG.append(str(m))


plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.titleweight": "bold", "axes.titlesize": 9})


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ═══════════════ load the six groupings on the SATURN gene universe ═══════════════
h = pd.read_csv(HPG, usecols=["key", "species", "gene", "mg_ESM1b_4sp", "mg_ESMC600M_raw_4sp", "mg_ESMC600M_resc_4sp"])
can = pd.concat([pd.read_csv(CANON / f"{s}_gene_to_macrogene.csv")[["gene", "macrogene"]].assign(species=s) for s in SPECIES])
can["key"] = can["species"] + "_" + can["gene"]
g = h.merge(can[["key", "macrogene"]].rename(columns={"macrogene": "mg_canon5sp_ESM1b"}), on="key", how="inner")

eg = pd.read_csv(EGG).dropna(subset=["metazoa_OG"])
eg = eg[eg["species"].isin(SPECIES)]
cnt = eg.groupby(["species", "gene_id", "metazoa_OG"]).size().reset_index(name="n")
cnt = cnt.sort_values(["species", "gene_id", "n", "metazoa_OG"], ascending=[True, True, False, True])
prim = cnt.drop_duplicates(["species", "gene_id"])[["species", "gene_id", "metazoa_OG"]]
g = g.merge(prim.rename(columns={"gene_id": "gene", "metazoa_OG": "mg_eggNOG"}), on=["species", "gene"], how="left")

of = pd.read_csv(OFDIR / "ALL4_gene_orthofinder_OG.csv", usecols=["species", "gene_id", "group"])
of = of.drop_duplicates(["species", "gene_id"]).rename(columns={"gene_id": "gene", "group": "mg_OrthoFinder"})
g = g.merge(of, on=["species", "gene"], how="left")
# OrthoFinder gene IDs can differ in case from SATURN names (06c: 99.9% exact, 100% case-insensitive)
miss = g["mg_OrthoFinder"].isna()
if miss.any():
    ofl = of.assign(gl=of["gene"].str.lower()).drop_duplicates(["species", "gl"])
    fill = g.loc[miss, ["species", "gene"]].assign(gl=lambda d: d["gene"].str.lower()).merge(
        ofl[["species", "gl", "mg_OrthoFinder"]], on=["species", "gl"], how="left")["mg_OrthoFinder"].values
    g.loc[miss, "mg_OrthoFinder"] = fill
g = g.reset_index(drop=True)
for G in GROUPINGS:
    g[G] = g["mg_" + G].astype("string")
log("=" * 78)
log("06d  Species overlap and orthogroup <-> macrogene nesting for six groupings")
log("=" * 78)
log(f"SATURN genes (all four runs): {len(g):,}  " + str(g['species'].value_counts().reindex(SPECIES).to_dict()))
log(f"with an eggNOG OG: {int(g['eggNOG'].notna().sum()):,};  with an OrthoFinder OG: {int(g['OrthoFinder'].notna().sum()):,}")
ANN = g["eggNOG"].notna().values
UNIVERSES = {"all SATURN genes": np.ones(len(g), bool), "eggNOG-annotated genes": ANN}

# ═══════════════ D0. OrthoFinder statistics ═══════════════
log("\n── D0. OrthoFinder statistics ──")
st = OFDIR / "Comparative_Genomics_Statistics"
if (st / "Statistics_PerSpecies.tsv").exists():
    lines = (st / "Statistics_PerSpecies.tsv").read_text().splitlines()
    block = []
    for ln in lines:
        if not ln.strip():
            if block: break
            continue
        block.append(ln.split("\t"))
    d0 = pd.DataFrame(block[1:], columns=["statistic"] + block[0][1:])
    d0.columns = ["statistic"] + [c.replace(".faa", "") for c in d0.columns[1:]]
    d0.to_csv(OUT / "D0_orthofinder_per_species.csv", index=False)
    log(d0.to_string(index=False))
    rows = []
    for typ in ["one-to-one", "one-to-many", "many-to-one", "many-to-many"]:
        f = st / f"OrthologuesStats_{typ}.tsv"
        if f.exists():
            m = pd.read_csv(f, sep="\t", index_col=0)
            m.index = [str(i).replace(".faa", "") for i in m.index]
            m.columns = [str(c).replace(".faa", "") for c in m.columns]
            for a in m.index:
                for b in m.columns:
                    if a != b:
                        rows.append({"type": typ, "species_A": a, "species_B": b, "n_genes_A": int(m.loc[a, b])})
    if rows:
        d0o = pd.DataFrame(rows)
        d0o.to_csv(OUT / "D0_orthofinder_orthologues.csv", index=False)
        log("Orthologue counts (genes in species A with orthologues of each type in species B):")
        log(d0o.pivot_table(index=["species_A", "species_B"], columns="type", values="n_genes_A").to_string())
else:
    log(f"(skipped: copy Comparative_Genomics_Statistics into {OFDIR} to include OrthoFinder's own statistics)")


# ═══════════════ helpers ═══════════════
def group_table(df, G):
    """per group: size, species set, n species (within the universe df)."""
    d = df[df[G].notna()]
    gt = d.groupby(G).agg(size=("species", "size"), sp=("species", lambda v: tuple(s for s in SPECIES if s in set(v))))
    gt["nsp"] = gt["sp"].map(len)
    return d, gt


def ent(codes):
    p = np.bincount(codes) / len(codes); p = p[p > 0]
    return float(-(p * np.log(p)).sum())


# ═══════════════ D1-D3. species overlap per grouping ═══════════════
log("\n── D1. Genes by number of species spanned by their group ──")
d1, d2, d3, d3b = [], [], [], []
for uni, mask in UNIVERSES.items():
    df = g[mask]
    for G in GROUPINGS:
        d, gt = group_table(df, G)
        gene_nsp = d[G].map(gt["nsp"]).values
        gene_sp = d[G].map(gt["sp"])
        row = {"grouping": LAB[G], "universe": uni, "genes_in_universe": len(df), "genes_assigned": len(d),
               "n_groups": len(gt), "n_single_gene_groups": int((gt["size"] == 1).sum()),
               "median_group_size": float(gt["size"].median()), "max_group_size": int(gt["size"].max()),
               "effective_n_groups": float(np.exp(ent(pd.factorize(d[G])[0])))}
        for k in range(1, 5):
            row[f"pct_genes_in_{k}sp_groups"] = 100 * float((gene_nsp == k).mean())
            row[f"n_groups_{k}sp"] = int((gt["nsp"] == k).sum())
        row["pct_genes_in_multisp_groups"] = 100 * float((gene_nsp >= 2).mean())
        d1.append(row)
        for c in COMBOS:
            sel = gt["sp"] == c
            d2.append({"grouping": LAB[G], "universe": uni, "combination": "+".join(c), "n_species": len(c),
                       "n_groups": int(sel.sum()), "n_genes": int(gt.loc[sel, "size"].sum())})
        for a in SPECIES:
            ma = (d["species"] == a).values
            if ma.sum() == 0:
                continue
            d3b.append({"grouping": LAB[G], "universe": uni, "species": a, "n_genes": int(ma.sum()),
                        "pct_in_multisp_group": 100 * float((gene_nsp[ma] >= 2).mean())})
            for b in SPECIES:
                if a != b:
                    d3.append({"grouping": LAB[G], "universe": uni, "species_A": a, "species_B": b,
                               "pct_A_genes_with_B_in_group": 100 * float(gene_sp[ma].map(lambda t: b in t).mean())})
D1, D2, D3, D3b = map(pd.DataFrame, (d1, d2, d3, d3b))
D1.to_csv(OUT / "D1_species_span_summary.csv", index=False)
D2.to_csv(OUT / "D2_species_combinations.csv", index=False)
D3.to_csv(OUT / "D3_species_pair_cooccurrence.csv", index=False)
D3b.to_csv(OUT / "D3b_per_species_multi.csv", index=False)
show = ["grouping", "universe", "genes_assigned", "n_groups", "median_group_size", "effective_n_groups",
        "pct_genes_in_1sp_groups", "pct_genes_in_2sp_groups", "pct_genes_in_3sp_groups", "pct_genes_in_4sp_groups",
        "n_groups_4sp"]
log(D1[show].round(1).to_string(index=False))
log("\n% of each species' genes in a group with >= 1 other species:")
log(D3b.pivot_table(index=["universe", "grouping"], columns="species", values="pct_in_multisp_group", sort=False)
    .reindex(columns=SPECIES).round(1).to_string())
log("\n% of species A's genes whose group contains species B (all SATURN genes; rows A, columns B):")
for G in GROUPINGS:
    t = D3[(D3.universe == "all SATURN genes") & (D3.grouping == LAB[G])].pivot(index="species_A", columns="species_B",
                                                                                values="pct_A_genes_with_B_in_group")
    log(f"  {LAB[G]}:\n" + t.reindex(index=SPECIES, columns=SPECIES).round(1).to_string())
log("\nGenes per species combination (eggNOG-annotated genes):")
log(D2[D2.universe == "eggNOG-annotated genes"].pivot_table(index="combination", columns="grouping", values="n_genes", sort=False)
    .reindex(columns=[LAB[G] for G in GROUPINGS]).to_string())


# ═══════════════ D4. pair agreement between groupings ═══════════════
log("\n── D4. Agreement on cross-species gene pairs ──")


def cross_pairs(df, G):
    """encoded (i<j) cross-species gene pairs placed in the same group; i, j = positional index in df."""
    d = df.reset_index(drop=True)
    d = d[d[G].notna()]
    sp = d["species"].values
    out = []
    for _, idx in d.groupby(G).indices.items():
        if len(idx) < 2:
            continue
        ii = d.index.values[idx]
        s = sp[idx]
        a, b = np.triu_indices(len(ii), 1)
        keep = s[a] != s[b]
        out.append(ii[a[keep]].astype(np.int64) * len(df) + ii[b[keep]])
    return np.unique(np.concatenate(out)) if out else np.array([], np.int64)


d4 = []
for uni, mask in UNIVERSES.items():
    df = g[mask].reset_index(drop=True)
    gs = [G for G in GROUPINGS if not (G == "eggNOG" and uni == "all SATURN genes")]
    P = {G: cross_pairs(df, G) for G in gs}
    for A in gs:
        for B in gs:
            inter = np.intersect1d(P[A], P[B], assume_unique=True).size
            d4.append({"universe": uni, "A": LAB[A], "B": LAB[B], "pairs_A": P[A].size, "pairs_B": P[B].size,
                       "shared": inter, "pct_A_pairs_in_B": 100 * inter / P[A].size if P[A].size else np.nan,
                       "jaccard": inter / (P[A].size + P[B].size - inter) if (P[A].size + P[B].size) else np.nan})
D4 = pd.DataFrame(d4)
D4.to_csv(OUT / "D4_pair_agreement.csv", index=False)
for uni in UNIVERSES:
    t = D4[D4.universe == uni]
    order = [LAB[G] for G in GROUPINGS if LAB[G] in set(t["A"])]
    log(f"\n{uni}: cross-species pairs per grouping: " + str(t.drop_duplicates("A").set_index("A")["pairs_A"].reindex(order).to_dict()))
    log("% of the row grouping's pairs also together in the column grouping:")
    log(t.pivot(index="A", columns="B", values="pct_A_pairs_in_B").reindex(index=order, columns=order).round(1).to_string())


# ═══════════════ D5. nesting vs each orthology reference ═══════════════
log("\n── D5. Orthogroups vs macrogenes (and vs each other): split / merge ──")


def nesting_scores(c, k):
    k = pd.factorize(k)[0]
    n = len(c)
    joint = pd.factorize(c.astype(np.int64) * (k.max() + 1) + k)[0]
    Hc, Hk, Hj = ent(c), ent(k), ent(joint)
    hom = 1 - (Hj - Hk) / Hc if Hc > 0 else 1.0
    comp = 1 - (Hj - Hc) / Hk if Hk > 0 else 1.0
    c2 = lambda x: (x * (x - 1) / 2).sum()
    nij = np.bincount(joint).astype(float)
    sij, sa, sb = c2(nij), c2(np.bincount(c).astype(float)), c2(np.bincount(k).astype(float))
    exp = sa * sb / (n * (n - 1) / 2)
    return {"homogeneity": hom, "completeness": comp, "v_measure": 2 * hom * comp / (hom + comp),
            "ARI": (sij - exp) / (0.5 * (sa + sb) - exp)}


def gene_categories(c, k, keep):
    k = pd.factorize(k)[0]
    df = pd.DataFrame({"og": c, "mg": k})
    split = df["og"].map(df.groupby("og")["mg"].nunique()).values > 1
    mixed = df["mg"].map(df.groupby("mg")["og"].nunique()).values > 1
    cat = np.where(~split & ~mixed, "match", np.where(~split & mixed, "OGs merged",
                   np.where(split & ~mixed, "OG subdivided", "cross-cutting")))
    return pd.Series(cat[keep]).value_counts(normalize=True).reindex(CATS, fill_value=0) * 100


d5 = []
for ref in REFS:
    for uni, mask in UNIVERSES.items():
        if ref == "eggNOG" and uni == "all SATURN genes":
            continue
        df = g[mask & g[ref].notna().values].reset_index(drop=True)
        c = pd.factorize(df[ref])[0]
        keep = np.bincount(c)[c] >= 2
        for G in GROUPINGS:
            if G == ref or (G == "eggNOG" and uni == "all SATURN genes"):
                continue
            sub = df[G].notna().values
            cc, kk, kp = pd.factorize(df.loc[sub, ref])[0], df.loc[sub, G].values, keep[sub]
            r = {"reference": LAB[ref], "universe": uni, "grouping": LAB[G], "n_genes": int(sub.sum()),
                 "n_genes_in_ref_groups_ge2": int(kp.sum())}
            r.update(nesting_scores(cc, kk))
            r.update(gene_categories(cc, kk, kp).to_dict())
            d5.append(r)
D5 = pd.DataFrame(d5)
D5.to_csv(OUT / "D5_nesting_vs_reference.csv", index=False)
log(D5[["reference", "universe", "grouping", "homogeneity", "completeness", "ARI"] + CATS].round(3).to_string(index=False))

# ═══════════════ figures ═══════════════
# FigD1: genes by no. of species spanned
fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), sharey=True)
for ax, uni in zip(axes, UNIVERSES):
    t = D1[D1.universe == uni].set_index("grouping").reindex([LAB[G] for G in GROUPINGS])
    y = np.arange(len(t))[::-1]
    left = np.zeros(len(t))
    for k in range(1, 5):
        v = t[f"pct_genes_in_{k}sp_groups"].values
        ax.barh(y, v, left=left, color=NSPCOL[k], label=f"{k} species", height=0.7)
        left += v
    ax.set_yticks(y); ax.set_yticklabels(t.index)
    ax.set_xlim(0, 100); ax.set_xlabel("% of genes (assigned)")
    ax.set_title(uni + (" (eggNOG: its 21,674 assigned genes)" if uni == "all SATURN genes" else ""))
axes[0].legend(frameon=False, ncol=4, loc="lower left", bbox_to_anchor=(0, 1.08), fontsize=8)
fig.suptitle("Genes by the number of species in their group", fontsize=9, color=MUTED, y=1.08)
save(fig, "FigD1_species_span")

# FigD2: genes per species combination (annotated genes; heatmap, % of genes)
t = D2[D2.universe == "eggNOG-annotated genes"].copy()
tot = t.groupby("grouping")["n_genes"].transform("sum")
t["pct"] = 100 * t["n_genes"] / tot
M = t.pivot_table(index="combination", columns="grouping", values="pct", sort=False).reindex(
    index=["+".join(c) for c in COMBOS], columns=[LAB[G] for G in GROUPINGS])
fig, ax = plt.subplots(figsize=(7.2, 5.6))
im = ax.imshow(M.values, cmap="Blues", aspect="auto", vmin=0)
ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(M.columns, rotation=30, ha="right")
ax.set_yticks(range(M.shape[0])); ax.set_yticklabels(M.index)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M.values[i, j]
        ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=7, color="white" if v > np.nanmax(M.values) * 0.6 else INK)
ax.set_title("% of genes in groups with each exact species combination (eggNOG-annotated genes)")
fig.colorbar(im, ax=ax, shrink=0.6, label="% of genes")
save(fig, "FigD2_species_combinations")

# FigD3: per-species % in multi-species groups
fig, axes = plt.subplots(1, 2, figsize=(10, 3.2), sharey=True)
for ax, uni in zip(axes, UNIVERSES):
    t = D3b[D3b.universe == uni]
    gs = [G for G in GROUPINGS if LAB[G] in set(t["grouping"])]
    w = 0.8 / len(gs)
    for i, G in enumerate(gs):
        v = t[t.grouping == LAB[G]].set_index("species").reindex(SPECIES)["pct_in_multisp_group"].values
        ax.bar(np.arange(4) + (i - (len(gs) - 1) / 2) * w, v, w, color=COL[G], label=LAB[G])
    ax.set_xticks(range(4)); ax.set_xticklabels(SPECIES); ax.set_ylim(0, 100)
    ax.set_title(uni); ax.set_ylabel("% of genes in a multi-species group")
axes[1].legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1, 1))
save(fig, "FigD3_per_species_multispecies")

# FigD4: pair agreement heatmap (annotated genes)
t = D4[D4.universe == "eggNOG-annotated genes"]
order = [LAB[G] for G in GROUPINGS]
M = t.pivot(index="A", columns="B", values="pct_A_pairs_in_B").reindex(index=order, columns=order)
fig, ax = plt.subplots(figsize=(6.4, 5.2))
im = ax.imshow(M.values, cmap="Purples", vmin=0, vmax=100)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=30, ha="right")
ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
for i in range(len(order)):
    for j in range(len(order)):
        v = M.values[i, j]
        ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=7.5, color="white" if v > 60 else INK)
ax.set_xlabel("... also together in this grouping"); ax.set_ylabel("% of this grouping's cross-species pairs ...")
ax.set_title("Agreement on cross-species gene pairs (eggNOG-annotated genes)")
fig.colorbar(im, ax=ax, shrink=0.7, label="%")
save(fig, "FigD4_pair_agreement")

# FigD5: split/merge categories vs each reference (annotated genes)
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.0))
fig.subplots_adjust(wspace=0.75, top=0.78)
for ax, ref in zip(axes, REFS):
    t = D5[(D5.reference == LAB[ref]) & (D5.universe == "eggNOG-annotated genes")].set_index("grouping")
    rows = [LAB[G] for G in GROUPINGS if LAB[G] in t.index]
    t = t.reindex(rows)
    y = np.arange(len(t))[::-1]
    left = np.zeros(len(t))
    for ccat in CATS:
        v = t[ccat].values
        ax.barh(y, v, left=left, color=CATCOL[ccat], label=ccat, height=0.7)
        left += v
    ax.set_yticks(y); ax.set_yticklabels(t.index); ax.set_xlim(0, 100)
    ax.set_xlabel("% of genes in reference groups of >= 2 genes")
    ax.set_title(f"Reference: {LAB[ref]}")
h_, l_ = axes[0].get_legend_handles_labels()
fig.legend(h_, l_, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.0), fontsize=8)
save(fig, "FigD5_nesting_categories")

(OUT / "summary.txt").write_text("\n".join(LOG) + "\n")
log(f"\nWrote tables and figures to {OUT}")
