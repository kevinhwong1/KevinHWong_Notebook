#!/usr/bin/env python3
"""
05j_shared_macrogenes_by_species.py  (runs locally on the Mac; no Pegasus needed)

How many macrogenes are shared between species, per model run?

Inputs (already downloaded):
  - 05h per-gene table: hard macrogene labels for ESM1b_4sp, ESMC600M_raw_4sp, ESMC600M_resc_4sp
  - canonical 5sp run: {SP}_gene_to_macrogene.csv (Spur dropped)
Universe: genes present in all four runs (so counts are comparable).
A species "is in" a macrogene if it contributes >= MIN_GENES genes (reported for 1 and 3).

Outputs (OUT folder):
  J1_species_pair_shared_mg.csv     per run x species pair: shared macrogenes, Jaccard of the
                                    two species' macrogene sets, % of each species' macrogenes shared
  J2_species_combination_counts.csv per run: macrogenes per exact species combination (UpSet-style)
  J5_single_species_macrogenes.csv per run x species: macrogenes containing ONLY that species,
                                    and how many genes sit in them
  J3_per_species_sharing.csv        per run x species: macrogenes containing the species, share of
                                    them that also contain >=1 other species, genes in shared macrogenes
  figures/FigJ1..J4
"""
import sys
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2] / "output" / "SATURN"
CMP = ROOT / "20261002_macrogene_comparison_ESM1b_vs_ESMC600M"
HPG = CMP / "20261006_norm_vs_macrogene" / "H_per_gene.csv"
CANON = ROOT / "20260711_SATURN_FINAL" / "macrogene"
OUT = CMP / "20261008_shared_macrogenes_by_species"
(OUT / "figures").mkdir(parents=True, exist_ok=True)

SPECIES = ["Mlei", "Cgig", "Crob", "Drer"]
RUNS = ["canon5sp_ESM1b", "ESM1b_4sp", "ESMC600M_raw_4sp", "ESMC600M_resc_4sp"]
RCOL = {"canon5sp_ESM1b": "#2a78d6", "ESM1b_4sp": "#eb6834", "ESMC600M_raw_4sp": "#1baf7a", "ESMC600M_resc_4sp": "#eda100"}
INK, MUTED = "#222222", "#8a8a85"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
PAIRS = list(combinations(SPECIES, 2))
LOG = []


def log(m=""):
    print(m); LOG.append(str(m))


h = pd.read_csv(HPG, usecols=["key", "species", "gene", "mg_ESM1b_4sp", "mg_ESMC600M_raw_4sp", "mg_ESMC600M_resc_4sp"])
can = pd.concat([pd.read_csv(CANON / f"{s}_gene_to_macrogene.csv")[["gene", "macrogene"]].assign(species=s) for s in SPECIES])
can["key"] = can["species"] + "_" + can["gene"]
g = h.merge(can[["key", "macrogene"]].rename(columns={"macrogene": "mg_canon5sp_ESM1b"}), on="key", how="inner")
log("Genes in all four runs: " + str(len(g)) + "  " + str(g["species"].value_counts().reindex(SPECIES).to_dict()))

pair_rows, combo_rows, sp_rows = [], [], []
for mn in (1, 3):
    for r in RUNS:
        cnt = pd.crosstab(g["mg_" + r], g["species"]).reindex(columns=SPECIES, fill_value=0)
        pres = cnt >= mn
        for a, b in PAIRS:
            both = int((pres[a] & pres[b]).sum()); na, nb = int(pres[a].sum()), int(pres[b].sum())
            pair_rows.append({"min_genes": mn, "run": r, "species_1": a, "species_2": b, "shared_mg": both,
                              "mg_with_sp1": na, "mg_with_sp2": nb, "jaccard": both / max(na + nb - both, 1),
                              "pct_of_sp1_mg_shared": 100 * both / max(na, 1), "pct_of_sp2_mg_shared": 100 * both / max(nb, 1)})
        combo = pres.apply(lambda row: "+".join(s for s in SPECIES if row[s]), axis=1)
        combo = combo[combo != ""]
        for c, n in combo.value_counts().items():
            combo_rows.append({"min_genes": mn, "run": r, "combination": c, "n_species": c.count("+") + 1, "n_macrogenes": int(n)})
        nsp = pres.sum(axis=1)
        for s in SPECIES:
            mine = pres[s]
            shared = mine & (nsp >= 2)
            genes_shared = int(cnt.loc[shared, s].sum()); genes_all = int(cnt[s].sum())
            sp_rows.append({"min_genes": mn, "run": r, "species": s, "mg_with_species": int(mine.sum()),
                            "mg_shared_with_other_species": int(shared.sum()),
                            "pct_mg_shared": 100 * shared.sum() / max(mine.sum(), 1),
                            "genes_in_shared_mg": genes_shared, "pct_genes_in_shared_mg": 100 * genes_shared / max(genes_all, 1),
                            "mg_shared_with_all_4": int((mine & (nsp == 4)).sum())})
J1, J2, J3 = pd.DataFrame(pair_rows), pd.DataFrame(combo_rows), pd.DataFrame(sp_rows)
# single-species macrogenes (species counted with >=1 gene, i.e. no other species present at all)
ss_rows = []
for r in RUNS:
    cnt = pd.crosstab(g["mg_" + r], g["species"]).reindex(columns=SPECIES, fill_value=0)
    only = (cnt > 0).sum(axis=1) == 1
    for s in SPECIES:
        m = only & (cnt[s] > 0)
        sz = cnt.loc[m, s]
        ss_rows.append({"run": r, "species": s, "n_single_species_mg": int(m.sum()),
                        "n_genes_in_them": int(sz.sum()), "pct_species_genes": 100 * sz.sum() / cnt[s].sum(),
                        "n_one_gene_mg": int((sz == 1).sum()), "median_size": float(sz.median()) if len(sz) else 0})
J5 = pd.DataFrame(ss_rows)
J5.to_csv(OUT / "J5_single_species_macrogenes.csv", index=False)
log("Single-species macrogenes (only one species present):")
log(J5.pivot(index="species", columns="run", values="n_single_species_mg")[RUNS].loc[SPECIES].to_string())
log("Genes in them (% of the species' genes):")
log(J5.pivot(index="species", columns="run", values="pct_species_genes")[RUNS].loc[SPECIES].round(1).to_string())
log("...of which one-gene macrogenes:")
log(J5.pivot(index="species", columns="run", values="n_one_gene_mg")[RUNS].loc[SPECIES].to_string())
J1.to_csv(OUT / "J1_species_pair_shared_mg.csv", index=False)
J2.to_csv(OUT / "J2_species_combination_counts.csv", index=False)
J3.to_csv(OUT / "J3_per_species_sharing.csv", index=False)

for mn in (1, 3):
    log(f"\n=== Species counted if >= {mn} gene(s) in the macrogene ===")
    d = J1[J1.min_genes == mn].assign(pair=lambda x: x.species_1 + "-" + x.species_2)
    log("Shared macrogenes per species pair:\n" + d.pivot(index="pair", columns="run", values="shared_mg")[RUNS].loc[[a + "-" + b for a, b in PAIRS]].to_string())
    log("Jaccard of the two species' macrogene sets:\n" + d.pivot(index="pair", columns="run", values="jaccard")[RUNS].loc[[a + "-" + b for a, b in PAIRS]].round(3).to_string())
    d3 = J3[J3.min_genes == mn]
    log("% of each species' genes in macrogenes shared with >=1 other species:\n" + d3.pivot(index="species", columns="run", values="pct_genes_in_shared_mg")[RUNS].loc[SPECIES].round(1).to_string())
    d2 = J2[J2.min_genes == mn]
    piv = d2.pivot_table(index="combination", columns="run", values="n_macrogenes", fill_value=0)[RUNS]
    piv = piv.loc[sorted(piv.index, key=lambda c: (-c.count("+"), c))]
    log("Macrogenes per exact species combination:\n" + piv.astype(int).to_string())

# ── Figures ──
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.titleweight": "bold", "axes.titlesize": 10})
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list("b", SEQ)
w = 0.8 / len(RUNS)
for mn in (1, 3):
    tag = "ge%dgene" % mn
    # J1: species x species heatmaps
    d = J1[J1.min_genes == mn]; d3 = J3[J3.min_genes == mn]
    vmax = d["shared_mg"].max()
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.9))
    for ax, r in zip(axes, RUNS):
        M = np.zeros((4, 4))
        for _, x in d[d.run == r].iterrows():
            i, j = SPECIES.index(x.species_1), SPECIES.index(x.species_2); M[i, j] = M[j, i] = x.shared_mg
        for i, s in enumerate(SPECIES):
            M[i, i] = np.nan
        im = ax.imshow(M, cmap=cmap, vmin=0, vmax=vmax)
        for i in range(4):
            for j in range(4):
                if i == j:
                    n = int(d3[(d3.run == r) & (d3.species == SPECIES[i])]["mg_with_species"].iloc[0])
                    ax.text(j, i, f"{n}\n(total)", ha="center", va="center", fontsize=7, color=MUTED)
                else:
                    ax.text(j, i, int(M[i, j]), ha="center", va="center", fontsize=9, color="white" if M[i, j] > 0.55 * vmax else INK)
        ax.set_xticks(range(4), SPECIES); ax.set_yticks(range(4), SPECIES); ax.set_title(r)
        for sp in ax.spines.values():
            sp.set_visible(False)
    fig.colorbar(im, ax=axes, fraction=0.015, label="Shared macrogenes")
    fig.suptitle(f"Macrogenes shared by each species pair (species counted if ≥{mn} gene{'s' if mn > 1 else ''})", fontweight="bold")
    fig.savefig(OUT / "figures" / f"FigJ1_pair_heatmaps_{tag}.png", dpi=200, bbox_inches="tight"); plt.close(fig)

    # J2: grouped bars per pair (% of the smaller species' macrogenes shared -> use Jaccard)
    fig, axes = plt.subplots(1, 2, figsize=(13, 3.9))
    x = np.arange(len(PAIRS)); lab = [a + " – " + b for a, b in PAIRS]
    for ax, col, ttl in [(axes[0], "shared_mg", "Shared macrogenes"), (axes[1], "jaccard", "Jaccard of the two species' macrogene sets")]:
        for k, r in enumerate(RUNS):
            v = d[d.run == r].set_index(["species_1", "species_2"]).loc[PAIRS, col]
            ax.bar(x - 0.4 + w / 2 + k * w, v, w * 0.95, color=RCOL[r], label=r)
        ax.set_xticks(x, lab, rotation=30, ha="right"); ax.set_title(ttl)
    h_, l_ = axes[0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="lower center", ncol=4, frameon=False)
    fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig(OUT / "figures" / f"FigJ2_pair_bars_{tag}.png", dpi=200); plt.close(fig)

    # J3: per species, % of genes in shared macrogenes
    fig, ax = plt.subplots(figsize=(7, 3.8))
    x = np.arange(4)
    for k, r in enumerate(RUNS):
        v = d3[d3.run == r].set_index("species").loc[SPECIES, "pct_genes_in_shared_mg"]
        ax.bar(x - 0.4 + w / 2 + k * w, v, w * 0.95, color=RCOL[r], label=r)
    ax.set_xticks(x, SPECIES); ax.set_ylim(0, 100)
    ax.set_ylabel("% of the species' genes in\nmacrogenes shared with other species")
    ax.set_title(f"Genes in cross-species macrogenes (≥{mn} gene per species)")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout(); fig.savefig(OUT / "figures" / f"FigJ3_per_species_{tag}.png", dpi=200); plt.close(fig)

    # J4: UpSet-style combination counts
    d2 = J2[J2.min_genes == mn]
    piv = d2.pivot_table(index="combination", columns="run", values="n_macrogenes", fill_value=0)[RUNS]
    order = sorted(piv.index, key=lambda c: (-c.count("+"), -piv.loc[c].mean()))
    piv = piv.loc[order]
    fig, (ax, axm) = plt.subplots(2, 1, figsize=(13, 5.2), gridspec_kw={"height_ratios": [3, 1.2]}, sharex=True)
    x = np.arange(len(order))
    for k, r in enumerate(RUNS):
        ax.bar(x - 0.4 + w / 2 + k * w, piv[r], w * 0.95, color=RCOL[r], label=r)
    ax.set_ylabel("Macrogenes"); ax.set_title(f"Macrogenes per exact species combination (≥{mn} gene per species)")
    ax.legend(frameon=False, fontsize=8, ncol=4)
    for i, c in enumerate(order):
        mem = c.split("+")
        ys = [3 - SPECIES.index(s) for s in mem]
        axm.scatter([i] * 4, range(4), s=40, color="#e1e0d9", zorder=1)
        axm.scatter([i] * len(ys), ys, s=40, color=INK, zorder=2)
        if len(ys) > 1:
            axm.plot([i, i], [min(ys), max(ys)], color=INK, lw=2, zorder=2)
    axm.set_yticks(range(4), SPECIES[::-1]); axm.set_xticks([])
    for sp in ["left", "bottom"]:
        axm.spines[sp].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT / "figures" / f"FigJ4_species_combinations_{tag}.png", dpi=200); plt.close(fig)

# J5: single-species macrogenes per species
fig, axes = plt.subplots(1, 2, figsize=(12, 3.9))
x = np.arange(4)
for ax, col, ttl in [(axes[0], "n_single_species_mg", "Macrogenes containing only one species"),
                     (axes[1], "pct_species_genes", "% of the species' genes in those macrogenes")]:
    for k, r in enumerate(RUNS):
        v = J5[J5.run == r].set_index("species").loc[SPECIES, col]
        ax.bar(x - 0.4 + w / 2 + k * w, v, w * 0.95, color=RCOL[r], label=r)
    ax.set_xticks(x, SPECIES); ax.set_title(ttl)
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "figures" / "FigJ5_single_species_macrogenes.png", dpi=200); plt.close(fig)

(OUT / "summary.txt").write_text("\n".join(LOG) + "\n")
print("Outputs:", OUT)
