#!/usr/bin/env python3
"""
20261030_notebook_figs_macrogene_2x2.py
Figures for the 2026-10-30 notebook entry (macrogene composition, 2x2: ESM-1b / ESMC post-norm x averaged /
longest isoform). Reads the downloaded 05g / 08c outputs; writes PNGs to the notebook images folder.
Run from mnemi_phagocyte_scRNAseq/:
  python3 output/SATURN/20261030_macrogene_2x2_prelim/20261030_notebook_figs_macrogene_2x2.py \
      --data output/SATURN/20261030_macrogene_2x2_prelim --img ../KevinHWong_Notebook/images
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--data", required=True)
ap.add_argument("--img", required=True)
a = ap.parse_args()
D, IMG = Path(a.data), Path(a.img)
IMG.mkdir(parents=True, exist_ok=True)

INK, INK2, GRID, SURF = "#1f1f1e", "#5f5e58", "#e6e5e0", "#fcfcfb"
MUTED = "#b9b8b2"
COL = {"ESM1b": "#2a78d6", "ESMC": "#eb6834"}          # validated categorical pair (light surface)
LAB = {"ESM1b": "ESM-1b", "ESMC": "ESMC post-norm"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": GRID,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK,
                     "figure.facecolor": SURF, "axes.facecolor": SURF})
SP = ["Mlei", "Cgig", "Dmel", "Crob", "Drer", "Hsap"]
NM = {"Mlei": "Mnemiopsis", "Cgig": "oyster", "Crob": "Ciona", "Drer": "zebrafish", "Dmel": "fly", "Hsap": "human"}
CONDS = ["ESM1b_avg", "ESMC_avg", "ESM1b_longest", "ESMC_longest"]
CLAB = {"ESM1b_avg": "ESM-1b\naveraged", "ESMC_avg": "ESMC\naveraged",
        "ESM1b_longest": "ESM-1b\nlongest", "ESMC_longest": "ESMC\nlongest"}


def clean(ax, keep_left=False):
    for s in ["top", "right"] + ([] if keep_left else ["left"]):
        ax.spines[s].set_visible(False)


def pair_name(p):
    s1, s2 = p.split("-")
    return f"{NM[s1]} – {NM[s2]}"


def dumbbell(ax, rows, avg, lon, color, label_min=0.05, fmt="{:+.2f}", off=0.012):
    """rows: labels top->bottom; avg/lon: arrays. Open marker = averaged, filled = longest."""
    y = np.arange(len(rows))
    for i in range(len(rows)):
        ax.plot([avg[i], lon[i]], [y[i], y[i]], color=color, lw=2, alpha=0.55, solid_capstyle="round", zorder=1)
    ax.scatter(avg, y, s=56, facecolor=SURF, edgecolor=color, linewidth=2, zorder=3)
    ax.scatter(lon, y, s=56, facecolor=color, edgecolor=SURF, linewidth=2, zorder=4)
    for i in range(len(rows)):
        d = lon[i] - avg[i]
        if abs(d) >= label_min:
            ax.text(max(avg[i], lon[i]) + off, y[i], fmt.format(d), va="center", fontsize=8.5, color=INK2)
    ax.set_yticks(y, rows)
    ax.invert_yaxis()
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    clean(ax)
    ax.tick_params(axis="y", length=0)


legend_handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=SURF, markeredgecolor=INK2,
                         markeredgewidth=2, markersize=8, label="all isoforms averaged"),
                  Line2D([0], [0], marker="o", color="w", markerfacecolor=INK2, markersize=8,
                         label="longest isoform")]

# ── Fig 1: recall per species pair, avg -> longest, both models ──────────────
p = pd.read_csv(D / "08c" / "C1_effects_by_pair_OG.csv")
r = p[p["metric"] == "recall"].copy()
r = r.sort_values("ESM1b_longest", ascending=False)
rows = [pair_name(x) for x in r["pair"]]
fig, axs = plt.subplots(1, 2, figsize=(10.4, 6.0), sharey=True, gridspec_kw={"wspace": 0.06})
for ax, m in zip(axs, ["ESM1b", "ESMC"]):
    dumbbell(ax, rows, r[f"{m}_avg"].values, r[f"{m}_longest"].values, COL[m], off=0.015)
    ax.set_title(LAB[m], loc="left", fontsize=11, color=INK)
    ax.set_xlabel("macrogene ortholog recall (OrthoFinder orthogroups)")
    ax.set_xlim(0, 0.52)
axs[1].tick_params(axis="y", labelleft=False)
fig.legend(handles=legend_handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.58, 1.0))
fig.subplots_adjust(left=0.2, right=0.98, top=0.88, bottom=0.1, wspace=0.06)
fig.savefig(IMG / "20261030_macrogene_recall_by_pair_2x2.png", dpi=200)
plt.close(fig)

# ── Fig 2: run agreement (ARI) between the four runs ─────────────────────────
r2 = pd.read_csv(D / "05g_OG" / "R2_pairwise_agreement.csv")
r2 = r2[r2["scope"] == "ALL"]
M = np.eye(4)
for _, x in r2.iterrows():
    i, j = CONDS.index(x["run_a"]), CONDS.index(x["run_b"])
    M[i, j] = M[j, i] = x["ARI"]
RAMP = ["#e8f1fc", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
cmap = LinearSegmentedColormap.from_list("blue", RAMP)
fig, ax = plt.subplots(figsize=(5.6, 4.8))
Mm = np.ma.masked_where(np.eye(4, dtype=bool), M)
im = ax.imshow(Mm, cmap=cmap, vmin=0, vmax=0.5)
for i in range(4):
    for j in range(4):
        if i == j:
            ax.text(j, i, "–", ha="center", va="center", color=INK2)
        else:
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=10,
                    color="white" if M[i, j] > 0.3 else INK)
ax.set_xticks(range(4), [CLAB[c] for c in CONDS], fontsize=9)
ax.set_yticks(range(4), [CLAB[c].replace("\n", " ") for c in CONDS], fontsize=9)
for s in ax.spines.values():
    s.set_visible(False)
ax.set_xticks(np.arange(-.5, 4), minor=True)
ax.set_yticks(np.arange(-.5, 4), minor=True)
ax.grid(which="minor", color=SURF, lw=3)
ax.tick_params(which="both", length=0)
cb = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.03)
cb.set_label("adjusted Rand index (macrogene assignment)", color=INK2)
cb.outline.set_visible(False)
fig.tight_layout()
fig.savefig(IMG / "20261030_macrogene_ARI_2x2.png", dpi=200)
plt.close(fig)

# ── Fig 3: partner stability (avg vs longest) by isoform count ───────────────
s = pd.read_csv(D / "08c" / "C2_partner_stability_summary.csv")
s = s[s["comparison"] == "isoform"]
CLS = ["1", "2-3", "4-9", "10+"]
ref = pd.read_csv(D / "08c" / "C2_partner_stability_summary.csv")
ref = ref[(ref["comparison"] == "model") & (ref["isoform_class"] == "ALL") & (ref["species"] == "ALL")]
fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.4), sharey=True, gridspec_kw={"wspace": 0.08})
for ax, m in zip(axs, ["ESM1b", "ESMC"]):
    d = s[s["run_a"] == f"{m}_avg"]
    labs = []
    for sp in SP:
        v = d[d["species"] == sp].set_index("isoform_class").reindex(CLS)["J_all_mean"]
        hl = sp == "Hsap"
        ax.plot(range(4), v.values, color=COL[m] if hl else MUTED, lw=2 if hl else 1.2,
                marker="o", markersize=6 if hl else 4, markeredgecolor=SURF, zorder=3 if hl else 2)
        if sp == "Mlei":   # one class only (all single-protein genes): label to the left of the point
            ax.text(-0.1, v.values[0], NM[sp], fontsize=8.5, va="center", ha="right", color=INK2)
        else:
            labs.append([v.values[-1], NM[sp], INK if hl else INK2])
    v = d[d["species"] == "ALL"].set_index("isoform_class").reindex(CLS)["J_all_mean"]
    ax.plot(range(4), v.values, color=INK, lw=2, ls=(0, (4, 2)), marker="s", markersize=5, zorder=4)
    labs.append([v.values[-1], "all species", INK])
    labs.sort(key=lambda t: t[0])
    for i in range(1, len(labs)):                      # keep end labels >= 0.017 apart
        labs[i][0] = max(labs[i][0], labs[i - 1][0] + 0.017)
    for yv, txt, c in labs:
        ax.text(3.12, yv, txt, fontsize=8.5, va="center", color=c)
    mv = ref[ref["run_a"].str.endswith("avg")]["J_all_mean"].values[0]
    ax.axhline(mv, color=INK2, lw=1, ls=":")
    ax.text(-0.15, mv - 0.018, f"ESM-1b vs ESMC, same isoform rule ({mv:.2f})", fontsize=8, color=INK2)
    ax.set_xticks(range(4), ["1\n(same protein)", "2–3", "4–9", "10+"])
    ax.set_xlabel("isoforms per gene")
    ax.set_xlim(-0.75, 3.95)
    ax.set_title(LAB[m] + ": averaged vs longest", loc="left", fontsize=11, color=INK)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    clean(ax, keep_left=True)
axs[0].set_ylabel("mean Jaccard of macrogene partners")
axs[0].set_ylim(0, 0.45)
fig.tight_layout()
fig.savefig(IMG / "20261030_partner_stability_by_isoforms.png", dpi=200)
plt.close(fig)

# ── Fig 4: share of genes in single-species macrogenes ───────────────────────
c3 = pd.read_csv(D / "08c" / "C3_single_species_fraction.csv").set_index("species")
fig, axs = plt.subplots(1, 2, figsize=(9.6, 4.2), sharey=True)
rows = [NM[x] for x in SP]
for ax, m in zip(axs, ["ESM1b", "ESMC"]):
    dumbbell(ax, rows, 100 * c3.loc[SP, f"{m}_avg"].values, 100 * c3.loc[SP, f"{m}_longest"].values, COL[m],
             label_min=1.0, fmt="{:+.1f}", off=0.8)
    ax.set_title(LAB[m], loc="left", fontsize=11, color=INK)
    ax.set_xlabel("% of the species' genes in single-species macrogenes")
    ax.set_xlim(0, 27)
axs[1].tick_params(axis="y", labelleft=False)
fig.legend(handles=legend_handles, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.56, 1.0))
fig.subplots_adjust(left=0.12, right=0.98, top=0.82, bottom=0.16, wspace=0.06)
fig.savefig(IMG / "20261030_single_species_macrogenes_2x2.png", dpi=200)
plt.close(fig)
print("figures written to", IMG)
