#!/usr/bin/env python3
"""
05m_esmc_gained_links.py  (runs locally; no Pegasus needed)

Are the cross-species links that only ESMC makes real, judged by evidence that
does not depend on eggNOG orthologous groups?

Cross-species gene pairs = two genes from different species in the same macrogene.
Categories (ESM-1b 4sp vs ESMC rescaled 4sp; raw ESMC reported too):
  ESM1b_only, both, ESMC_only, random (cross-species pairs drawn at random = baseline)

Evidence per pair (from eggNOG-mapper annotations, gene level = union over proteins):
  share_OG        same Metazoa orthologous group (the benchmark used before)
  share_Pfam      at least one Pfam domain in common (domain-level homology; weaker than orthology,
                  but can catch distant relationships that OG assignment misses)
  same_COG        same COG functional category (broad function; excludes S = unknown)
  same_name_root  gene names share a family root (letters before the first digit, >= 3 letters: KLHL7 / KLHL10 -> KLHL; SLC5A1 / SLC16A2 -> SLC)
Each metric is computed only over pairs where both genes have that annotation.

Outputs (OUT folder): M1 evidence summary, M2 annotation coverage, M3 gained macrogenes with names,
M4 all ESMC-only pairs with annotations, figures/FigM1, FigM2.
"""
import re
from itertools import combinations
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
OUT = CMP / "20261008_esmc_gained_links"
(OUT / "figures").mkdir(parents=True, exist_ok=True)
SPECIES = ["Mlei", "Cgig", "Crob", "Drer"]
A, B, R = "ESM1b_4sp", "ESMC600M_resc_4sp", "ESMC600M_raw_4sp"
LOG = []


def log(m=""):
    print(m); LOG.append(str(m))


g = pd.read_csv(HPG, usecols=["species", "gene", "mg_" + A, "mg_" + B, "mg_" + R]).reset_index(drop=True)
N = len(g)
sp = g["species"].values

# ── gene-level annotation ──
ann = []
for s in SPECIES:
    e = pd.read_csv(EGG / f"{s}_eggnog_annotations.csv",
                    usecols=["gene_id", "Preferred_name", "Description", "COG_category", "eggNOG_OGs", "PFAMs", "GOs"])
    e["og"] = e["eggNOG_OGs"].astype(str).str.extract(r"([^,]+)@33208\|Metazoa")[0]
    e["immune"] = e["GOs"].astype(str).str.contains("GO:0002376|GO:0006952")
    agg = e.groupby("gene_id").agg(
        og=("og", lambda v: frozenset(v.dropna())),
        pfam=("PFAMs", lambda v: frozenset(x for s_ in v.dropna().astype(str) for x in s_.split(",") if x and x != "-")),
        cog=("COG_category", lambda v: frozenset(c for s_ in v.dropna().astype(str) for c in s_ if c not in "-S")),
        name=("Preferred_name", lambda v: next((x for x in v.dropna() if x != "-"), np.nan)),
        desc=("Description", "first"), immune=("immune", "any")).reset_index()
    ann.append(agg.assign(species=s))
ann = pd.concat(ann).rename(columns={"gene_id": "gene"})
g = g.merge(ann, on=["species", "gene"], how="left")
g["annotated"] = g["og"].notna()
for c in ["og", "pfam", "cog"]:
    g[c] = g[c].apply(lambda v: v if isinstance(v, frozenset) else frozenset())
g["immune"] = g["immune"].astype("boolean").fillna(False).astype(bool)


def root(n):
    if not isinstance(n, str):
        return None
    m2 = re.match(r"^([A-Za-z]{3,})", n.upper())   # letters-only root: KLHL7 -> KLHL, SLC5A1 -> SLC, TBX20 -> TBX
    return m2.group(1) if m2 else None


g["root"] = g["name"].apply(root)
OG, PF, CG, RT = g["og"].values, g["pfam"].values, g["cog"].values, g["root"].values


def pairs(lab):
    df = pd.DataFrame({"i": np.arange(N), "mg": lab, "sp": sp})
    out = []
    for a, b in combinations(SPECIES, 2):
        m = df[df.sp == a].merge(df[df.sp == b], on="mg")
        out.append(m["i_x"].values.astype(np.int64) * N + m["i_y"].values.astype(np.int64))
    return set(np.concatenate(out).tolist())


PA, PB, PR = pairs(g["mg_" + A].values), pairs(g["mg_" + B].values), pairs(g["mg_" + R].values)
rng = np.random.default_rng(42)
rand = set()
idx_sp = {s: np.where(sp == s)[0] for s in SPECIES}
for a, b in combinations(SPECIES, 2):
    ia = rng.choice(idx_sp[a], 40000); ib = rng.choice(idx_sp[b], 40000)
    rand |= set((ia.astype(np.int64) * N + ib).tolist())
CATS = {"ESM-1b only": PA - PB, "both": PA & PB, "ESMC only": PB - PA, "random pairs": rand}
CATS_RAW = {"ESMC raw only": PR - PA, "both (raw)": PA & PR}


def evidence(codes):
    c = np.fromiter(codes, dtype=np.int64)
    i, j = c // N, c % N
    res = {"n_pairs": len(c)}
    both_ann = g["annotated"].values[i] & g["annotated"].values[j]
    res["pct_both_annotated"] = 100 * both_ann.mean()
    for key, arr in [("OG", OG), ("Pfam", PF), ("COG", CG)]:
        has = np.array([bool(arr[x]) and bool(arr[y]) for x, y in zip(i, j)])
        sh = np.array([bool(arr[x] & arr[y]) for x, y in zip(i[has], j[has])])
        res[f"pct_share_{key}"] = 100 * sh.mean() if len(sh) else np.nan
        res[f"n_with_{key}"] = int(has.sum())
    hasr = np.array([RT[x] is not None and RT[y] is not None for x, y in zip(i, j)])
    res["pct_same_name_root"] = 100 * np.mean([RT[x] == RT[y] for x, y in zip(i[hasr], j[hasr])]) if hasr.any() else np.nan
    mlei = (sp[i] == "Mlei") | (sp[j] == "Mlei")
    res["pct_pairs_with_Mlei"] = 100 * mlei.mean()
    drer = (sp[i] == "Drer") | (sp[j] == "Drer")
    res["pct_pairs_with_Drer"] = 100 * drer.mean()
    return res


rows = []
for scope, sel in [("all pairs", None), ("Mlei pairs", "Mlei"), ("Drer pairs", "Drer"), ("invertebrate pairs", "inv")]:
    for k, v in {**CATS, **CATS_RAW}.items():
        c = np.fromiter(v, dtype=np.int64)
        i, j = c // N, c % N
        if sel == "Mlei":
            c = c[(sp[i] == "Mlei") | (sp[j] == "Mlei")]
        elif sel == "Drer":
            c = c[(sp[i] == "Drer") | (sp[j] == "Drer")]
        elif sel == "inv":
            c = c[(sp[i] != "Drer") & (sp[j] != "Drer")]
        rows.append({"scope": scope, "category": k, **evidence(set(c.tolist()))})
M1 = pd.DataFrame(rows)
M1.to_csv(OUT / "M1_pair_evidence.csv", index=False)
cols = ["category", "n_pairs", "pct_both_annotated", "pct_share_OG", "pct_share_Pfam", "pct_share_COG", "pct_same_name_root"]
for scope in M1.scope.unique():
    log(f"\n=== {scope}: cross-species pairs in the same macrogene ===")
    log(M1[M1.scope == scope][cols].round(1).to_string(index=False))

# ── M2: genes that gained partners — annotation coverage ──
def partners(lab):
    df = pd.DataFrame({"i": np.arange(N), "mg": lab, "sp": sp})
    out = np.zeros(N, bool)
    for _, grp in df.groupby("mg"):
        if grp.sp.nunique() > 1:
            out[grp.i.values] = True
    return out


xa, xb = partners(g["mg_" + A].values), partners(g["mg_" + B].values)
g["status"] = np.select([xa & xb, xa & ~xb, ~xa & xb], ["has partners in both", "lost", "gained"], "species-specific")
M2 = g.groupby(["species", "status"]).agg(n=("gene", "size"), pct_eggNOG_annotated=("annotated", "mean"),
                                         pct_named=("name", lambda v: v.notna().mean()),
                                         n_immune=("immune", "sum")).reset_index()
M2[["pct_eggNOG_annotated", "pct_named"]] *= 100
M2.to_csv(OUT / "M2_annotation_by_status.csv", index=False)
log("\n=== Annotation coverage by gene status (ESM-1b -> ESMC rescaled) ===")
log(M2.round(1).to_string(index=False))

# ── M3: gained macrogenes = ESMC macrogenes whose cross-species links are mostly new ──
lb = g["mg_" + B].values
rows = []
for mg, grp in g.groupby("mg_" + B):
    if grp.species.nunique() < 2:
        continue
    ii = grp.index.values
    cp = [(min(x, y), max(x, y)) for x, y in combinations(ii, 2) if sp[x] != sp[y]]
    # pair codes in PA are oriented by species order; test both orientations
    new = sum(1 for x, y in cp if (x * N + y) not in PA and (y * N + x) not in PA)
    sh_og = np.mean([bool(OG[x] & OG[y]) for x, y in cp if OG[x] and OG[y]]) if any(OG[x] and OG[y] for x, y in cp) else np.nan
    sh_pf = np.mean([bool(PF[x] & PF[y]) for x, y in cp if PF[x] and PF[y]]) if any(PF[x] and PF[y] for x, y in cp) else np.nan
    pf = pd.Series([d for x in ii for d in PF[x]]).value_counts()
    rows.append({"mg_ESMC_resc": mg, "n_genes": len(grp), **{"n_" + s: int((grp.species == s).sum()) for s in SPECIES},
                 "n_xsp_pairs": len(cp), "pct_pairs_new": 100 * new / len(cp),
                 "pct_pairs_share_OG": 100 * sh_og if sh_og == sh_og else np.nan,
                 "pct_pairs_share_Pfam": 100 * sh_pf if sh_pf == sh_pf else np.nan,
                 "n_unannotated": int((~grp.annotated).sum()), "n_immune": int(grp.immune.sum()),
                 "top_Pfam": "|".join(pf.index[:5]),
                 "names": "|".join(grp["name"].dropna().value_counts().index[:10]),
                 "genes": "|".join((grp.species + ":" + grp.gene.astype(str)).values[:15])})
M3 = pd.DataFrame(rows).sort_values(["pct_pairs_new", "n_genes"], ascending=[False, False])
M3.to_csv(OUT / "M3_esmc_macrogenes_by_novelty.csv", index=False)
gained = M3[M3.pct_pairs_new >= 80]
log(f"\n=== ESMC-rescaled cross-species macrogenes: {len(M3)}; with >=80% new cross-species pairs: {len(gained)} ===")
log("Among the mostly-new macrogenes: median % pairs sharing a Pfam domain = "
    f"{gained.pct_pairs_share_Pfam.median():.1f}; sharing an OG = {gained.pct_pairs_share_OG.median():.1f}")
log("Largest mostly-new macrogenes with high Pfam support (>= 50% of annotated pairs share a domain):")
log(gained[gained.pct_pairs_share_Pfam >= 50].sort_values("n_genes", ascending=False).head(15)[
    ["mg_ESMC_resc", "n_genes", "n_Mlei", "n_Cgig", "n_Crob", "n_Drer", "pct_pairs_share_Pfam", "pct_pairs_share_OG", "n_unannotated", "top_Pfam", "names"]].round(0).to_string(index=False))
log("Largest mostly-new macrogenes with low Pfam support (< 10%):")
log(gained[gained.pct_pairs_share_Pfam < 10].sort_values("n_genes", ascending=False).head(10)[
    ["mg_ESMC_resc", "n_genes", "n_Mlei", "n_Cgig", "n_Crob", "n_Drer", "pct_pairs_share_Pfam", "top_Pfam", "names"]].round(0).to_string(index=False))

# ── M4: ESMC-only pairs table ──
c = np.fromiter(CATS["ESMC only"], dtype=np.int64)
i, j = c // N, c % N
M4 = pd.DataFrame({"species_1": sp[i], "gene_1": g.gene.values[i], "name_1": g.name.values[i],
                   "species_2": sp[j], "gene_2": g.gene.values[j], "name_2": g.name.values[j],
                   "mg_ESMC_resc": lb[i],
                   "share_OG": [bool(OG[x] & OG[y]) if OG[x] and OG[y] else np.nan for x, y in zip(i, j)],
                   "share_Pfam": [bool(PF[x] & PF[y]) if PF[x] and PF[y] else np.nan for x, y in zip(i, j)],
                   "shared_Pfam": ["|".join(sorted(PF[x] & PF[y])) for x, y in zip(i, j)]})
M4.to_csv(OUT / "M4_esmc_only_pairs.csv", index=False)

# ── figures ──
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#8a8a85", "axes.titleweight": "bold", "axes.titlesize": 10})
CC = {"ESM-1b only": "#eb6834", "both": "#2a78d6", "ESMC only": "#eda100", "random pairs": "#b4b2a9"}
mets = [("pct_share_OG", "Same orthologous group"), ("pct_share_Pfam", "Share a Pfam domain"),
        ("pct_share_COG", "Same COG category"), ("pct_same_name_root", "Same gene-name family")]
scopes = ["all pairs", "Mlei pairs", "invertebrate pairs", "Drer pairs"]
fig, axes = plt.subplots(1, 4, figsize=(16, 3.8), sharey=True)
x = np.arange(len(mets)); w = 0.2
for ax, scope in zip(axes, scopes):
    d = M1[M1.scope == scope].set_index("category")
    for k, cat in enumerate(CC):
        ax.bar(x - 0.3 + k * w, d.loc[cat, [m for m, _ in mets]], w * 0.95, color=CC[cat], label=cat)
    ax.set_xticks(x, [t for _, t in mets], rotation=30, ha="right"); ax.set_title(scope)
axes[0].set_ylabel("% of annotated pairs")
h_, l_ = axes[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="lower center", ncol=4, frameon=False)
fig.suptitle("Evidence for cross-species gene pairs placed in the same macrogene", fontweight="bold")
fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(OUT / "figures" / "FigM1_pair_evidence.png", dpi=200); plt.close(fig)

(OUT / "summary.txt").write_text("\n".join(LOG) + "\n")
print("Outputs:", OUT)
