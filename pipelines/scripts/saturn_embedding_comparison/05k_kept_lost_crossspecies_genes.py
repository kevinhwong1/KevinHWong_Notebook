#!/usr/bin/env python3
"""
05k_kept_lost_crossspecies_genes.py  (runs locally; no Pegasus needed)

Which genes keep or lose their cross-species macrogene partners when the model
changes from ESM-1b (4sp) to ESMC-600M (rescaled 4sp; raw also reported)?

For every gene (common genes of the 4sp runs, from 05h's H_per_gene.csv):
  xsp partners = genes from OTHER species in the same macrogene
  status (ESM-1b -> ESMC rescaled):
     kept   : has cross-species partners in both runs AND shares >=1 of its ESM-1b partners
     moved  : has cross-species partners in both runs but none of the same ones
     lost   : cross-species partners in ESM-1b, none in ESMC (gene now in a single-species macrogene)
     gained : none in ESM-1b, some in ESMC
     species-specific : no cross-species partners in either run (single-species macrogene in both)
Drer genes additionally get a vertebrate<->invertebrate version (partners from Mlei/Cgig/Crob only).

Annotation: eggNOG (output/Eggnog/{SP}_eggnog_annotations.csv): Preferred_name, Description,
COG category, Metazoa OG; immune flag = GO:0002376 (immune system process) or GO:0006952
(defense response) among the gene's propagated GO terms.

Outputs (OUT folder):
  K1_gene_status.csv              per gene: status, partners, annotation
  K2_drer_invert_pairs.csv        every Drer-invertebrate co-assigned pair in either run, with
                                  in_ESM1b / in_ESMC flags and whether the pair shares a Metazoa OG
  K3_status_summary.csv           counts per species x status (+ immune, + OG-supported)
  K4_COG_enrichment.csv           COG categories over-represented among lost vs kept (Fisher, BH)
  K5_drer_macrogene_fate.csv      each ESM-1b macrogene containing Drer + invertebrates: genes,
                                  names, and how many Drer genes keep invertebrate partners in ESMC
  figures/FigK1_status_by_species.png, FigK2_COG_lost_vs_kept_all_species.png
"""
from math import comb, log10
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2] / "output"
CMP = ROOT / "SATURN" / "20261002_macrogene_comparison_ESM1b_vs_ESMC600M"
HPG = CMP / "20261006_norm_vs_macrogene" / "H_per_gene.csv"
EGG = ROOT / "Eggnog"
OUT = CMP / "20261008_kept_lost_crossspecies_genes"
(OUT / "figures").mkdir(parents=True, exist_ok=True)
SPECIES = ["Mlei", "Cgig", "Crob", "Drer"]
INVERT = ["Mlei", "Cgig", "Crob"]
A, B, R = "ESM1b_4sp", "ESMC600M_resc_4sp", "ESMC600M_raw_4sp"
LOG = []
STATUSES = ["kept", "moved", "lost", "gained", "species-specific"]


def log(m=""):
    print(m); LOG.append(str(m))


g = pd.read_csv(HPG, usecols=["key", "species", "gene", "mg_" + A, "mg_" + B, "mg_" + R])
g = g.reset_index(drop=True)
N = len(g)
log("Genes: " + str(N))

# ── eggNOG annotation, one row per gene (first annotated protein) ──
ann = []
for s in SPECIES:
    e = pd.read_csv(EGG / f"{s}_eggnog_annotations.csv",
                    usecols=["gene_id", "Preferred_name", "Description", "COG_category", "eggNOG_OGs", "GOs"])
    e["metazoa_OG"] = e["eggNOG_OGs"].astype(str).str.extract(r"([^,]+)@33208\|Metazoa")[0]
    e["immune"] = e["GOs"].astype(str).str.contains("GO:0002376|GO:0006952")
    agg = e.groupby("gene_id").agg(Preferred_name=("Preferred_name", "first"), Description=("Description", "first"),
                                   COG_category=("COG_category", "first"),
                                   metazoa_OG=("metazoa_OG", lambda v: ";".join(sorted(set(v.dropna())))),
                                   immune=("immune", "any")).reset_index()
    ann.append(agg.assign(species=s))
ann = pd.concat(ann).rename(columns={"gene_id": "gene"})
g = g.merge(ann, on=["species", "gene"], how="left")
g["immune"] = g["immune"].fillna(False).astype(bool)
g["Preferred_name"] = g["Preferred_name"].replace("-", np.nan)
sp = g["species"].values


def partners(labels, allowed_species=None):
    """dict gene_index -> set of cross-species gene indices in the same macrogene"""
    df = pd.DataFrame({"i": np.arange(N), "mg": labels, "sp": sp})
    out = {i: set() for i in range(N)}
    for _, grp in df.groupby("mg"):
        if grp["sp"].nunique() < 2:
            continue
        idx, sps = grp["i"].values, grp["sp"].values
        for i, si in zip(idx, sps):
            m = sps != si
            if allowed_species is not None:
                m &= np.isin(sps, allowed_species)
            out[i] = set(idx[m])
    return out


PA, PB, PR = partners(g["mg_" + A].values), partners(g["mg_" + B].values), partners(g["mg_" + R].values)


def status(pa, pb):
    if pa and pb:
        return "kept" if pa & pb else "moved"
    if pa:
        return "lost"
    if pb:
        return "gained"
    return "species-specific"


g["status_" + B] = [status(PA[i], PB[i]) for i in range(N)]
g["status_" + R] = [status(PA[i], PR[i]) for i in range(N)]
g["n_xsp_partners_" + A] = [len(PA[i]) for i in range(N)]
g["n_xsp_partners_" + B] = [len(PB[i]) for i in range(N)]
g["n_partners_kept"] = [len(PA[i] & PB[i]) for i in range(N)]
# OG support: does the gene share a Metazoa OG with any of its ESM-1b partners?
ogs = g["metazoa_OG"].fillna("").str.split(";").apply(lambda v: set(x for x in v if x))
g["ESM1b_partner_shares_OG"] = [any(ogs[i] & ogs[j] for j in PA[i]) if ogs[i] else np.nan for i in range(N)]
g["ESMC_partner_shares_OG"] = [any(ogs[i] & ogs[j] for j in PB[i]) if ogs[i] else np.nan for i in range(N)]


def names(idx_set, k=8):
    v = [f"{sp[j]}:{g.at[j, 'Preferred_name'] if isinstance(g.at[j, 'Preferred_name'], str) else g.at[j, 'gene']}" for j in sorted(idx_set)]
    return "|".join(v[:k]) + ("|..." if len(v) > k else "")


g["ESM1b_partners_example"] = [names(PA[i]) for i in range(N)]
g["ESMC_resc_partners_example"] = [names(PB[i]) for i in range(N)]

# Drer: vertebrate <-> invertebrate links only
PAi, PBi = partners(g["mg_" + A].values, INVERT), partners(g["mg_" + B].values, INVERT)
d = sp == "Drer"
g["drer_invert_status"] = np.where(d, [status(PAi[i], PBi[i]) for i in range(N)], "")
cols = ["species", "gene", "Preferred_name", "Description", "COG_category", "metazoa_OG", "immune",
        "status_" + B, "status_" + R, "drer_invert_status", "n_xsp_partners_" + A, "n_xsp_partners_" + B,
        "n_partners_kept", "ESM1b_partner_shares_OG", "ESMC_partner_shares_OG",
        "ESM1b_partners_example", "ESMC_resc_partners_example", "mg_" + A, "mg_" + B, "mg_" + R]
g[cols].to_csv(OUT / "K1_gene_status.csv", index=False)

# ── K3 summary ──
rows = []
for s in SPECIES:
    for st in STATUSES:
        m = (sp == s) & (g["status_" + B] == st)
        rows.append({"species": s, "status": st, "n_genes": int(m.sum()), "pct": 100 * m.sum() / (sp == s).sum(),
                     "n_immune": int(g.loc[m, "immune"].sum()),
                     "pct_with_OG_partner_in_ESM1b": 100 * g.loc[m, "ESM1b_partner_shares_OG"].dropna().astype(bool).mean()
                     if g.loc[m, "ESM1b_partner_shares_OG"].notna().any() else np.nan})
k3 = pd.DataFrame(rows)
k3.to_csv(OUT / "K3_status_summary.csv", index=False)
log("\nStatus of each gene's cross-species partners, ESM-1b 4sp -> ESMC rescaled 4sp (% of species' genes):")
log(k3.pivot(index="species", columns="status", values="pct")[STATUSES].loc[SPECIES].round(1).to_string())
log("\n% of genes whose ESM-1b cross-species partners include an eggNOG Metazoa-OG partner, by status:")
log(k3.pivot(index="species", columns="status", values="pct_with_OG_partner_in_ESM1b")[["kept", "moved", "lost"]].loc[SPECIES].round(1).to_string())
dd = g[d]
log("\nDrer vertebrate<->invertebrate links (ESM-1b -> ESMC rescaled):")
log(dd["drer_invert_status"].value_counts().to_string())
lost_imm = dd[(dd["drer_invert_status"] == "lost") & dd["immune"]]
log(f"\nDrer immune-annotated genes (GO immune system process / defense response): {int(dd['immune'].sum())}; "
    f"of these, invertebrate links lost under ESMC: {len(lost_imm)}")
log("Examples (Drer gene: ESM-1b invertebrate partners):")
for _, x in lost_imm.sort_values("n_xsp_partners_" + A, ascending=False).head(25).iterrows():
    nm = x["Preferred_name"] if isinstance(x["Preferred_name"], str) else ""
    log(f"  {x['gene']} ({nm}; {str(x['Description'])[:60]}): {x['ESM1b_partners_example']}")

# ── K2 Drer-invertebrate pairs ──
def pairs(lab):
    df = pd.DataFrame({"i": np.arange(N), "mg": lab, "sp": sp})
    a = df[df.sp == "Drer"]; b = df[df.sp.isin(INVERT)]
    m = a.merge(b, on="mg")
    return set(zip(m.i_x, m.i_y))


pA, pB = pairs(g["mg_" + A].values), pairs(g["mg_" + B].values)
allp = sorted(pA | pB)
k2 = pd.DataFrame(allp, columns=["i", "j"])
k2["drer_gene"] = g.loc[k2.i, "gene"].values
k2["drer_name"] = g.loc[k2.i, "Preferred_name"].values
k2["invert_species"] = g.loc[k2.j, "species"].values
k2["invert_gene"] = g.loc[k2.j, "gene"].values
k2["invert_name"] = g.loc[k2.j, "Preferred_name"].values
k2["in_ESM1b"] = [p in pA for p in allp]
k2["in_ESMC_resc"] = [p in pB for p in allp]
k2["shares_metazoa_OG"] = [bool(ogs[i] & ogs[j]) if ogs[i] and ogs[j] else np.nan for i, j in allp]
k2["drer_immune"] = g.loc[k2.i, "immune"].values
k2.drop(columns=["i", "j"]).to_csv(OUT / "K2_drer_invert_pairs.csv", index=False)
log("\nDrer-invertebrate co-assigned pairs: ESM-1b " + str(len(pA)) + ", ESMC rescaled " + str(len(pB)) + ", both " + str(len(pA & pB)))
for lab, m in [("only ESM-1b", k2.in_ESM1b & ~k2.in_ESMC_resc), ("both", k2.in_ESM1b & k2.in_ESMC_resc),
               ("only ESMC", ~k2.in_ESM1b & k2.in_ESMC_resc)]:
    v = k2.loc[m, "shares_metazoa_OG"].dropna()
    log(f"  {lab:12s}: {int(m.sum()):6d} pairs; {100 * v.astype(bool).mean():.1f}% share a Metazoa OG (of {len(v)} with both annotated)")

# ── K4 COG enrichment: lost vs kept ──
def fisher_greater(a, b, c, d_):
    # one-sided P(X >= a) for 2x2 [[a,b],[c,d]] (hypergeometric)
    n1, n2, k = a + b, c + d_, a + c
    tot = comb(n1 + n2, k)
    return sum(comb(n1, x) * comb(n2, k - x) for x in range(a, min(n1, k) + 1)) / tot


rows = []
for s in SPECIES + ["ALL"]:
    sub = g if s == "ALL" else g[sp == s]
    sub = sub[sub["status_" + B].isin(["lost", "kept"]) & sub["COG_category"].notna() & (sub["COG_category"] != "-")]
    cats = sorted(set("".join(sub["COG_category"].astype(str))))
    L, Kp = sub[sub["status_" + B] == "lost"], sub[sub["status_" + B] == "kept"]
    for c in cats:
        a = int(L["COG_category"].str.contains(c, regex=False).sum()); b = len(L) - a
        cc = int(Kp["COG_category"].str.contains(c, regex=False).sum()); dd_ = len(Kp) - cc
        rows.append({"species": s, "COG": c, "lost_with": a, "lost_total": len(L), "kept_with": cc, "kept_total": len(Kp),
                     "pct_lost": 100 * a / max(len(L), 1), "pct_kept": 100 * cc / max(len(Kp), 1),
                     "p_enriched_in_lost": fisher_greater(a, b, cc, dd_), "p_enriched_in_kept": fisher_greater(cc, dd_, a, b)})
k4 = pd.DataFrame(rows)
for col in ["p_enriched_in_lost", "p_enriched_in_kept"]:
    k4[col.replace("p_", "padj_")] = np.nan
    for s, idx in k4.groupby("species").groups.items():
        p = k4.loc[idx, col].values; n = len(p); o = np.argsort(p)
        adj = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
        out = np.empty(n); out[o] = np.minimum(adj, 1)
        k4.loc[idx, col.replace("p_", "padj_")] = out
k4.to_csv(OUT / "K4_COG_enrichment.csv", index=False)
COG = {"A": "RNA processing", "B": "Chromatin", "C": "Energy production", "D": "Cell cycle", "E": "Amino acid metab.",
       "F": "Nucleotide metab.", "G": "Carbohydrate metab.", "H": "Coenzyme metab.", "I": "Lipid metab.",
       "J": "Translation", "K": "Transcription", "L": "Replication/repair", "M": "Cell wall/membrane",
       "N": "Motility", "O": "PTM/chaperones", "P": "Inorganic ion transport", "Q": "Secondary metab.",
       "S": "Unknown function", "T": "Signal transduction", "U": "Trafficking/secretion", "V": "Defense",
       "W": "Extracellular structures", "Y": "Nuclear structure", "Z": "Cytoskeleton"}
sig = k4[(k4.padj_enriched_in_lost < 0.05) | (k4.padj_enriched_in_kept < 0.05)].copy()
sig["name"] = sig["COG"].map(COG)
log("\nCOG categories differing between lost and kept genes (BH < 0.05):")
log(sig[["species", "COG", "name", "pct_lost", "pct_kept", "padj_enriched_in_lost", "padj_enriched_in_kept"]]
    .round(4).to_string(index=False) if len(sig) else "  none")

# ── K5 fate of ESM-1b macrogenes with Drer + invertebrates ──
rows = []
for mg, grp in g.groupby("mg_" + A):
    if not ((grp.species == "Drer").any() and grp.species.isin(INVERT).any()):
        continue
    dr = grp[grp.species == "Drer"].index
    kept = sum(1 for i in dr if PAi[i] & PBi[i])
    nm = grp["Preferred_name"].dropna()
    rows.append({"mg_ESM1b": mg, "n_genes": len(grp), **{"n_" + s: int((grp.species == s).sum()) for s in SPECIES},
                 "n_drer_keeping_invert_partners": kept, "frac_drer_kept": kept / len(dr),
                 "n_immune": int(grp["immune"].sum()),
                 "top_names": "|".join(nm.value_counts().index[:8]),
                 "drer_genes": "|".join(grp.loc[dr, "gene"].astype(str)[:10])})
k5 = pd.DataFrame(rows).sort_values(["frac_drer_kept", "n_genes"], ascending=[True, False])
k5.to_csv(OUT / "K5_drer_macrogene_fate.csv", index=False)
log(f"\nESM-1b macrogenes with Drer + invertebrates: {len(k5)}; Drer fully separated from its invertebrate partners in ESMC: "
    f"{int((k5.frac_drer_kept == 0).sum())}; fully kept: {int((k5.frac_drer_kept == 1).sum())}")
log("Largest fully-separated examples:")
log(k5[k5.frac_drer_kept == 0].sort_values("n_genes", ascending=False).head(12)[["mg_ESM1b", "n_genes", "n_Mlei", "n_Cgig", "n_Crob", "n_Drer", "n_immune", "top_names"]].to_string(index=False))

# ── Figures ──
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#8a8a85", "axes.titleweight": "bold", "axes.titlesize": 10})
SC = {"kept": "#2a78d6", "moved": "#9ec5f4", "lost": "#e34948", "gained": "#1baf7a", "species-specific": "#d3d1c7"}
fig, ax = plt.subplots(figsize=(7, 3.8))
left = np.zeros(4)
for st in STATUSES:
    v = k3[k3.status == st].set_index("species").loc[SPECIES, "pct"].values
    ax.barh(SPECIES[::-1], v[::-1], left=left[::-1], color=SC[st], label=st, edgecolor="white", linewidth=1.5)
    left += v
ax.set_xlabel("% of the species' genes"); ax.set_xlim(0, 100)
ax.set_title("Cross-species macrogene partners: ESM-1b → ESMC (rescaled)")
ax.legend(frameon=False, ncol=5, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18))
fig.tight_layout(); fig.savefig(OUT / "figures" / "FigK1_status_by_species.png", dpi=200); plt.close(fig)

fig, axes = plt.subplots(2, 2, figsize=(13, 11))
for ax, s_ in zip(axes.ravel(), SPECIES):
    dk = k4[k4.species == s_].copy(); dk["name"] = dk["COG"].map(COG).fillna(dk["COG"])
    dk = dk[(dk.lost_with + dk.kept_with) >= 20].sort_values("pct_lost")
    y = np.arange(len(dk))
    ax.barh(y + 0.2, dk["pct_kept"], 0.4, color="#2a78d6", label=f"kept (n={int(dk.kept_total.iloc[0]) if len(dk) else 0})")
    ax.barh(y - 0.2, dk["pct_lost"], 0.4, color="#e34948", label=f"lost (n={int(dk.lost_total.iloc[0]) if len(dk) else 0})")
    for yi, (_, x) in zip(y, dk.iterrows()):
        if x.padj_enriched_in_lost < 0.05 or x.padj_enriched_in_kept < 0.05:
            ax.text(max(x.pct_lost, x.pct_kept) + 0.3, yi, "*", va="center", fontsize=12)
    ax.set_yticks(y, dk["COG"] + "  " + dk["name"], fontsize=8); ax.set_xlabel("% of genes in the group")
    ax.set_title(s_); ax.legend(frameon=False, fontsize=8, loc="lower right")
fig.suptitle("Genes losing vs keeping cross-species partners (ESM-1b → ESMC rescaled), by COG category  (* BH < 0.05)", fontweight="bold")
fig.tight_layout(); fig.savefig(OUT / "figures" / "FigK2_COG_lost_vs_kept_all_species.png", dpi=200); plt.close(fig)

(OUT / "summary.txt").write_text("\n".join(LOG) + "\n")
print("Outputs:", OUT)

# ── K6: macrogene-level stability ESM-1b -> ESMC rescaled, with names (examples of what stays / changes) ──
la, lb = g["mg_" + A].values, g["mg_" + B].values
ct = pd.crosstab(la, lb)
size_a, size_b = ct.sum(axis=1), ct.sum(axis=0)
best = ct.idxmax(axis=1)
rows = []
for mg in ct.index:
    b_ = best[mg]; ov = ct.at[mg, b_]
    jac = ov / (size_a[mg] + size_b[b_] - ov)
    grp = g[la == mg]
    nsp_a = grp["species"].nunique()
    nsp_b = g.loc[lb == b_, "species"].nunique()
    kept_genes = grp[lb[la == mg] == b_]
    rows.append({"mg_ESM1b": mg, "n_genes": len(grp), "n_species_ESM1b": nsp_a,
                 **{"n_" + s: int((grp.species == s).sum()) for s in SPECIES},
                 "best_mg_ESMC_resc": b_, "n_genes_best": int(size_b[b_]), "n_species_best": nsp_b,
                 "overlap": int(ov), "jaccard": jac,
                 "species_kept_together": "+".join(s for s in SPECIES if (kept_genes.species == s).any()),
                 "n_immune": int(grp["immune"].sum()),
                 "names": "|".join(grp["Preferred_name"].dropna().value_counts().index[:10])})
k6 = pd.DataFrame(rows)
k6["status"] = np.where(k6.jaccard >= 0.5, "stable", np.where(k6.jaccard >= 0.2, "partial", "reorganized"))
k6.to_csv(OUT / "K6_macrogene_stability_ESM1b_to_ESMCresc.csv", index=False)
log("\nESM-1b macrogenes by fate in ESMC rescaled (best-match Jaccard: stable >= 0.5, partial 0.2-0.5, reorganized < 0.2):")
log(pd.crosstab(k6["n_species_ESM1b"], k6["status"]).to_string())
log("\nExamples: stable 4-species macrogenes (largest):")
log(k6[(k6.status == "stable") & (k6.n_species_ESM1b == 4) & (k6.species_kept_together.str.count(r"\+") == 3)]
    .sort_values("n_genes", ascending=False).head(12)[["mg_ESM1b", "n_genes", "n_Mlei", "n_Cgig", "n_Crob", "n_Drer", "jaccard", "names"]].round(2).to_string(index=False))
log("\nExamples: 4-species macrogenes where Drer split off but invertebrates stayed together:")
inv = k6[(k6.n_species_ESM1b == 4) & (k6.species_kept_together == "Mlei+Cgig+Crob") & (k6.jaccard >= 0.3)]
log(inv.sort_values("n_genes", ascending=False).head(12)[["mg_ESM1b", "n_genes", "n_Mlei", "n_Cgig", "n_Crob", "n_Drer", "jaccard", "names"]].round(2).to_string(index=False))
log("\nImmune-annotated genes by status (ESM-1b -> ESMC rescaled):")
log(pd.crosstab(g.loc[g.immune, "species"], g.loc[g.immune, "status_" + B]).reindex(SPECIES)[STATUSES].to_string())
(OUT / "summary.txt").write_text("\n".join(LOG) + "\n")
