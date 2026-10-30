#!/usr/bin/env python3
"""
08c_macrogene_2x2.py   (2026-10-30)

Macrogene composition across the 2x2 design: protein embedding model (ESM-1b / ESMC-600M post-norm) x
protein per gene (all isoforms averaged / longest isoform). Complements 05g (run descriptives R1, run
agreement R2, ortholog benchmark R3), whose output folders it reads.

Macrogene membership is fixed at the end of SATURN pretraining (the genes_to_macrogenes pkl is written then),
so macrogene-only runs and full runs with the same inputs are comparable. Hard assignment = argmax of each
gene's macrogene weights (as 04p / 05g). Gene-level parts use the genes common to all four runs.

Run labels must be <MODEL>_<ISOFORM> with MODEL in {ESM1b, ESMC} and ISOFORM in {avg, longest}.

C1  effects        from 05g R3 (OG and HOG): per scope and per species pair, the four conditions plus
                   isoform effect within each model (longest - avg), model effect within each isoform rule
                   (ESMC - ESM1b) and the interaction ((ESMC_longest - ESMC_avg) - (ESM1b_longest - ESM1b_avg)).
C2  partner stability  per gene, Jaccard of its macrogene partners (other genes in its macrogene) between two
                   runs: all partners and cross-species partners only. Comparisons: isoform rule within each
                   model (main), plus model within each isoform rule (reference). Stratified by species and by
                   number of isoforms (ESM-1b protein set: 03_embeddings/{SP}/{SP}_gene_to_protein.filtered.json).
                   Single-isoform genes have the same protein under both rules, so their change between avg and
                   longest is the run-to-run baseline. Optional: cosine(avg, longest) of each gene's embedding.
C3  species structure  fraction of each species' genes in single-species macrogenes per condition (+ effects);
                   human genes that leave a human-only macrogene with longest isoforms (per model), with symbol
                   and isoform count.

Outputs (--out_dir): C1_effects_summary_{OG,HOG}.csv, C1_effects_by_pair_{OG,HOG}.csv,
  C2_partner_stability_by_gene.csv, C2_partner_stability_summary.csv,
  C3_single_species_fraction.csv, C3_Hsap_genes_leaving_human_only_mg_{ESM1b,ESMC}.csv, C_summary.txt

Usage:
  python 08c_macrogene_2x2.py --runs ESM1b_avg=DIR,ESMC_avg=DIR,ESM1b_longest=DIR,ESMC_longest=DIR \
      --og_dir <05g OG out> --hog_dir <05g HOG out> --out_dir DIR \
      [--emb ESM1b_avg=EMBDIR,ESMC_avg=EMBDIR,...]
"""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--runs", required=True)
ap.add_argument("--og_dir", required=True)
ap.add_argument("--hog_dir", default="")
ap.add_argument("--out_dir", required=True)
ap.add_argument("--emb", default="", help="optional label=embedding_dir for the embedding-change column")
ap.add_argument("--base", default="/scratch/dark_genes/SATURN_Mnemi")
ap.add_argument("--species", default="Mlei,Cgig,Crob,Drer,Dmel,Hsap")
args = ap.parse_args()

BASE = Path(args.base)
OUT = Path(args.out_dir)
OUT.mkdir(parents=True, exist_ok=True)
SPECIES = args.species.split(",")
CONDS = ["ESM1b_avg", "ESMC_avg", "ESM1b_longest", "ESMC_longest"]
LOG = []


def log(m=""):
    print(m, flush=True)
    LOG.append(str(m))


RUNS = dict(item.split("=", 1) for item in args.runs.split(","))
if sorted(RUNS) != sorted(CONDS):
    raise SystemExit(f"--runs labels must be exactly {CONDS}; got {sorted(RUNS)}")
log("=" * 70 + "\n08c macrogene composition, 2x2 (model x protein per gene)\n" + "=" * 70)
for c in CONDS:
    log(f"{c:14s} {RUNS[c]}")


# ── effects helper ───────────────────────────────────────────────────────────
def add_effects(w):
    w["iso_effect_ESM1b"] = w["ESM1b_longest"] - w["ESM1b_avg"]
    w["iso_effect_ESMC"] = w["ESMC_longest"] - w["ESMC_avg"]
    w["model_effect_avg"] = w["ESMC_avg"] - w["ESM1b_avg"]
    w["model_effect_longest"] = w["ESMC_longest"] - w["ESM1b_longest"]
    w["interaction"] = w["iso_effect_ESMC"] - w["iso_effect_ESM1b"]
    return w


# ════════════════════════════════════════════════════════════
# C1. effects on the 05g ortholog benchmark
# ════════════════════════════════════════════════════════════
log("\n── C1. Isoform and model effects on the ortholog benchmark (05g R3) ──")
for ref, d in [("OG", args.og_dir), ("HOG", args.hog_dir)]:
    if not d or not (Path(d) / "R3_ortholog_benchmark_summary.csv").exists():
        log(f"  {ref}: no 05g output at '{d}', skipped")
        continue
    s = pd.read_csv(Path(d) / "R3_ortholog_benchmark_summary.csv")
    rows = []
    for metric in ["recall_pooled", "precision_pooled", "F1_pooled", "og_recovery_mean"]:
        w = s.pivot(index="scope", columns="run", values=metric).reindex(columns=CONDS)
        w = add_effects(w)
        w.insert(0, "metric", metric)
        rows.append(w.reset_index())
    eff = pd.concat(rows, ignore_index=True)
    eff.to_csv(OUT / f"C1_effects_summary_{ref}.csv", index=False)
    log(f"\n  [{ref}] pooled over all species pairs:")
    log(eff[eff["scope"] == "all_pairs"].drop(columns="scope").round(3).to_string(index=False))
    p = pd.read_csv(Path(d) / "R3_ortholog_benchmark_by_species_pair.csv")
    p["pair"] = p["species_1"] + "-" + p["species_2"]
    rows = []
    for metric in ["recall", "precision", "F1", "og_recovery"]:
        w = p.pivot(index="pair", columns="run", values=metric).reindex(columns=CONDS)
        w = add_effects(w)
        w.insert(0, "metric", metric)
        rows.append(w.reset_index())
    pe = pd.concat(rows, ignore_index=True)
    pe.to_csv(OUT / f"C1_effects_by_pair_{ref}.csv", index=False)
    if ref == "OG":
        r = pe[pe["metric"] == "recall"].sort_values("iso_effect_ESM1b", ascending=False)
        log("\n  [OG] recall by species pair, sorted by the ESM-1b isoform effect:")
        log(r[["pair"] + CONDS + ["iso_effect_ESM1b", "iso_effect_ESMC"]].round(3).to_string(index=False))


# ════════════════════════════════════════════════════════════
# load macrogene assignments
# ════════════════════════════════════════════════════════════
def to_np(v):
    if hasattr(v, "detach"):
        v = v.detach().cpu().numpy()
    return np.asarray(v, dtype=np.float32).ravel()


def load_mg(run_dir):
    hits = sorted((Path(run_dir) / "saturn_results").glob("*genes_to_macrogenes*.pkl"))
    if len(hits) != 1:
        raise SystemExit(f"{run_dir}: expected 1 genes_to_macrogenes pkl, found {len(hits)}")
    raw = pickle.load(open(hits[0], "rb"))
    keys = [k for k in raw if k.split("_", 1)[0] in SPECIES]
    W = np.vstack([to_np(raw[k]) for k in keys])
    return pd.Series(W.argmax(1), index=keys), W.shape[1]


MG, K = {}, {}
for c in CONDS:
    MG[c], K[c] = load_mg(RUNS[c])
common = sorted(set.intersection(*[set(MG[c].index) for c in CONDS]))
sp_of = np.array([k.split("_", 1)[0] for k in common])
gene_of = [k.split("_", 1)[1] for k in common]
A = {c: MG[c].loc[common].values for c in CONDS}
N = len(common)
log(f"\nGenes common to all four runs: {N:,}  " +
    str(pd.Series(sp_of).value_counts().reindex(SPECIES).to_dict()))

# isoform counts (ESM-1b protein set; the matched ESMC averages the same set)
niso = {}
for sp in SPECIES:
    f = BASE / "03_embeddings" / sp / f"{sp}_gene_to_protein.filtered.json"
    if f.exists():
        for g, ps in json.load(open(f)).items():
            niso[f"{sp}_{g}"] = len(ps)
n_iso = np.array([niso.get(k, np.nan) for k in common], dtype=float)
log(f"isoform counts found for {np.isfinite(n_iso).sum():,} of {N:,} common genes")


def iso_class(n):
    if not np.isfinite(n):
        return "unknown"
    return "1" if n == 1 else "2-3" if n <= 3 else "4-9" if n <= 9 else "10+"


iso_cls = np.array([iso_class(n) for n in n_iso])

# optional: embedding change between avg and longest
emb_cos = {}
if args.emb:
    import torch
    EMB = dict(item.split("=", 1) for item in args.emb.split(","))

    def load_emb(d, sp):
        for f in [sp, "Mnemi"] if sp == "Mlei" else [sp]:
            p = Path(d) / f / f"{f}_gene_embeddings.pt"
            if p.exists():
                return torch.load(p, map_location="cpu")
        return {}
    for model in ["ESM1b", "ESMC"]:
        a, b = EMB.get(f"{model}_avg"), EMB.get(f"{model}_longest")
        if not (a and b):
            continue
        cos = np.full(N, np.nan)
        for sp in SPECIES:
            da, db = load_emb(a, sp), load_emb(b, sp)
            idx = np.where(sp_of == sp)[0]
            for i in idx:
                g = gene_of[i]
                if g in da and g in db:
                    va, vb = da[g].float().flatten(), db[g].float().flatten()
                    cos[i] = float(torch.nn.functional.cosine_similarity(va, vb, dim=0))
        emb_cos[model] = cos
        log(f"embedding cos(avg, longest) {model}: median {np.nanmedian(cos):.4f} "
            f"(multi-isoform genes {np.nanmedian(cos[n_iso > 1]):.4f})")


# ════════════════════════════════════════════════════════════
# C2. per-gene partner stability
# ════════════════════════════════════════════════════════════
log("\n── C2. Per-gene macrogene partner stability (Jaccard of macrogene partners) ──")
sp_codes = pd.Categorical(sp_of, categories=SPECIES).codes


def partner_jaccard(a, b):
    """all partners and cross-species partners; a, b = macrogene index per gene (same gene order)."""
    Ka, Kb = a.max() + 1, b.max() + 1
    joint = pd.Series(1, index=pd.MultiIndex.from_arrays([a, b])).groupby(level=[0, 1]).sum()
    inter = joint.reindex(pd.MultiIndex.from_arrays([a, b])).values - 1
    sa = np.bincount(a, minlength=Ka)[a] - 1
    sb = np.bincount(b, minlength=Kb)[b] - 1
    union = sa + sb - inter
    j_all = np.where(union > 0, inter / np.maximum(union, 1), 1.0)
    # same-species counts, to remove them
    js = pd.Series(1, index=pd.MultiIndex.from_arrays([a, b, sp_codes])).groupby(level=[0, 1, 2]).sum()
    inter_s = js.reindex(pd.MultiIndex.from_arrays([a, b, sp_codes])).values - 1
    sa_s = pd.Series(1, index=pd.MultiIndex.from_arrays([a, sp_codes])).groupby(level=[0, 1]).sum() \
        .reindex(pd.MultiIndex.from_arrays([a, sp_codes])).values - 1
    sb_s = pd.Series(1, index=pd.MultiIndex.from_arrays([b, sp_codes])).groupby(level=[0, 1]).sum() \
        .reindex(pd.MultiIndex.from_arrays([b, sp_codes])).values - 1
    ic, ac, bc = inter - inter_s, sa - sa_s, sb - sb_s
    uc = ac + bc - ic
    j_cross = np.where(uc > 0, ic / np.maximum(uc, 1), np.nan)   # NaN: no cross-species partner in either run
    return j_all, j_cross


COMPS = [("isoform", "ESM1b_avg", "ESM1b_longest"), ("isoform", "ESMC_avg", "ESMC_longest"),
         ("model", "ESM1b_avg", "ESMC_avg"), ("model", "ESM1b_longest", "ESMC_longest")]
gene_df = pd.DataFrame({"key": common, "species": sp_of, "gene": gene_of, "n_isoforms": n_iso,
                        "isoform_class": iso_cls})
for c in CONDS:
    gene_df[f"mg_{c}"] = A[c]
for m, cos in emb_cos.items():
    gene_df[f"emb_cos_avg_longest_{m}"] = cos
summ = []
for kind, ra, rb in COMPS:
    ja, jc = partner_jaccard(A[ra], A[rb])
    tag = f"{ra}_vs_{rb}"
    gene_df[f"J_all_{tag}"] = ja
    gene_df[f"J_cross_{tag}"] = jc
    for sp in SPECIES + ["ALL"]:
        msp = np.ones(N, bool) if sp == "ALL" else sp_of == sp
        for cls in ["ALL", "1", "2-3", "4-9", "10+"]:
            m = msp & ((iso_cls == cls) if cls != "ALL" else True)
            if m.sum() == 0:
                continue
            summ.append({"comparison": kind, "run_a": ra, "run_b": rb, "species": sp, "isoform_class": cls,
                         "n_genes": int(m.sum()),
                         "J_all_median": float(np.median(ja[m])), "J_all_mean": float(ja[m].mean()),
                         "frac_identical_partners": float((ja[m] == 1).mean()),
                         "J_cross_median": float(np.nanmedian(jc[m])) if np.isfinite(jc[m]).any() else np.nan,
                         "J_cross_mean": float(np.nanmean(jc[m])) if np.isfinite(jc[m]).any() else np.nan})
gene_df.to_csv(OUT / "C2_partner_stability_by_gene.csv", index=False)
summ = pd.DataFrame(summ)
summ.to_csv(OUT / "C2_partner_stability_summary.csv", index=False)
log("  isoform rule (avg vs longest) within each model; class '1' = single-isoform genes = baseline:")
log(summ[summ["comparison"] == "isoform"].pivot_table(
    index=["species", "isoform_class"], columns="run_a", values="J_all_mean").round(3).to_string())
log("\n  model (ESM-1b vs ESMC) within each isoform rule, all genes:")
log(summ[(summ["comparison"] == "model") & (summ["isoform_class"] == "ALL")]
    [["run_a", "run_b", "species", "J_all_mean", "J_cross_mean"]].round(3).to_string(index=False))


# ════════════════════════════════════════════════════════════
# C3. species structure: single-species macrogenes
# ════════════════════════════════════════════════════════════
log("\n── C3. Fraction of genes in single-species macrogenes (common genes) ──")
nsp_gene = {}
rows = []
for c in CONDS:
    ct = pd.crosstab(A[c], sp_of)
    nsp = (ct > 0).sum(axis=1)
    nsp_gene[c] = nsp.reindex(A[c]).values
    for sp in SPECIES + ["ALL"]:
        m = np.ones(N, bool) if sp == "ALL" else sp_of == sp
        rows.append({"run": c, "species": sp, "frac_single_species_mg": float((nsp_gene[c][m] == 1).mean()),
                     "frac_in_6sp_mg": float((nsp_gene[c][m] == len(SPECIES)).mean())})
c3 = pd.DataFrame(rows)
w = c3.pivot(index="species", columns="run", values="frac_single_species_mg").reindex(columns=CONDS)
w = add_effects(w).reindex(SPECIES + ["ALL"])
w.reset_index().to_csv(OUT / "C3_single_species_fraction.csv", index=False)
log(w.round(3).to_string())

sym = {}
try:
    import anndata as ad
    a = ad.read_h5ad(BASE / "02_processed_h5ad" / "Hsap_saturn.h5ad", backed="r")
    if "gene_symbol" in a.var.columns:
        sym = a.var["gene_symbol"].astype(str).to_dict()
    a.file.close()
except Exception as e:  # symbols are a convenience only
    log(f"  (no Hsap gene symbols: {e})")
for model in ["ESM1b", "ESMC"]:
    m = (sp_of == "Hsap") & (nsp_gene[f"{model}_avg"] == 1) & (nsp_gene[f"{model}_longest"] > 1)
    back = (sp_of == "Hsap") & (nsp_gene[f"{model}_avg"] > 1) & (nsp_gene[f"{model}_longest"] == 1)
    d = gene_df.loc[m, ["gene", "n_isoforms", f"mg_{model}_avg", f"mg_{model}_longest"]].copy()
    d.insert(1, "symbol", d["gene"].map(sym))
    d["n_species_in_longest_mg"] = nsp_gene[f"{model}_longest"][m]
    if f"emb_cos_avg_longest_{model}" in gene_df:
        d["emb_cos_avg_longest"] = gene_df.loc[m, f"emb_cos_avg_longest_{model}"].values
    d.sort_values("n_isoforms", ascending=False).to_csv(
        OUT / f"C3_Hsap_genes_leaving_human_only_mg_{model}.csv", index=False)
    hs = sp_of == "Hsap"
    log(f"  {model}: Hsap genes human-only (avg) -> shared (longest): {m.sum():,}; reverse: {back.sum():,}; "
        f"median isoforms of movers {np.nanmedian(n_iso[m]) if m.any() else float('nan'):.0f} "
        f"vs all Hsap {np.nanmedian(n_iso[hs]):.0f}")

(OUT / "C_summary.txt").write_text("\n".join(LOG) + "\n")
log(f"\nwritten -> {OUT}")
