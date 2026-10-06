#!/usr/bin/env python3
"""
05c_compare_macrogenes_ESM1b_vs_ESMC.py

Label-invariant comparison of SATURN macrogene assignments between two
macrogene-only runs (default: 4-species ESM-1b vs ESMC-600M, 2026-10-02),
plus an anchor comparison against the canonical 5-species ESM-1b run.

Macrogene IDs are arbitrary per run (KMeans init), so nothing here compares
macrogene numbers directly. Everything is based on which genes are grouped
together.

Sections
  0  Diagnostics      : gene universes, weight-matrix stats, macrogene sizes
  A  Partition        : ARI/AMI/NMI, macrogene best-match + Hungarian matching,
                        per-gene co-member stability (all + cross-species),
                        soft-weight kNN neighbourhood overlap
  B  Species structure: n_species tier distributions, tier transitions per gene,
                        species-pair shared macrogenes, cross-species
                        co-assigned gene pairs
  R  Reference anchor : where canonical 5sp (ESM-1b) macrogenes land in each
                        4sp run; focus macrogenes traced gene by gene

Hard assignment = argmax of each gene's macrogene weight vector (same
convention as 04p_macrogene_summary.sh).

Inputs  : {run}/saturn_results/*genes_to_macrogenes*.pkl
Outputs : --out_dir (CSV tables, PNG figures, summary.txt)
"""
import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import (adjusted_mutual_info_score, adjusted_rand_score,
                             normalized_mutual_info_score)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUNS = Path("/scratch/dark_genes/SATURN_Mnemi/04_saturn_runs")

ap = argparse.ArgumentParser()
ap.add_argument("--run_a", default=str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESM1b_macrogene"))
ap.add_argument("--run_b", default=str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESMC600M_macrogene"))
ap.add_argument("--label_a", default="ESM1b")
ap.add_argument("--label_b", default="ESMC600M")
ap.add_argument("--ref_run", default=str(RUNS / "20260617_Mlei_Cgig_Spur_Crob_Drer_K3000"),
                help="Canonical run used as an anchor; skipped if no pkl is found")
ap.add_argument("--ref_label", default="canon5sp_ESM1b")
ap.add_argument("--focus_ref_mgs", default="1089,2003",
                help="Canonical macrogene IDs to trace gene by gene")
ap.add_argument("--species", default="Mlei,Cgig,Crob,Drer")
ap.add_argument("--knn_k", type=int, default=10)
ap.add_argument("--skip_knn", action="store_true")
ap.add_argument("--out_dir", default=str(RUNS / "20261002_macrogene_comparison_ESM1b_vs_ESMC600M"))
args = ap.parse_args()

SPECIES = args.species.split(",")
LA, LB = args.label_a, args.label_b
OUT = Path(args.out_dir)
(OUT / "figures").mkdir(parents=True, exist_ok=True)

COL_A, COL_B = "#2a78d6", "#eb6834"   # categorical slots 1, 2
INK, MUTED = "#222222", "#8a8a85"

SUMMARY = []


def log(msg=""):
    print(msg, flush=True)
    SUMMARY.append(msg)


def save(df, name, index=False):
    df.to_csv(OUT / name, index=index)
    print("  wrote " + name, flush=True)


# ════════════════════════════════════════════════════════════
# Loading
# ════════════════════════════════════════════════════════════
def to_np(v):
    if hasattr(v, "detach"):
        v = v.detach().cpu().numpy()
    return np.asarray(v, dtype=np.float32).ravel()


def find_pkl(run_dir):
    return sorted((Path(run_dir) / "saturn_results").glob("*genes_to_macrogenes*.pkl"))


def load_g2m(run_dir, label, keep_species=None):
    hits = find_pkl(run_dir)
    if len(hits) == 0:
        raise FileNotFoundError("No genes_to_macrogenes pkl in " + str(Path(run_dir) / "saturn_results"))
    if len(hits) > 1:
        raise RuntimeError("Multiple genes_to_macrogenes pkls for " + label + ": " + str(hits))
    log(label + ": " + str(hits[0]))
    raw = pickle.load(open(hits[0], "rb"))
    keys = list(raw.keys())
    sp = [k.split("_", 1)[0] if "_" in k else "unknown" for k in keys]
    gene = [k.split("_", 1)[1] if "_" in k else k for k in keys]
    df = pd.DataFrame({"key": keys, "species": sp, "gene": gene})
    W = np.vstack([to_np(raw[k]) for k in keys])
    if keep_species is not None:
        m = df["species"].isin(keep_species).values
        df, W = df.loc[m].reset_index(drop=True), W[m]
    df["mg"] = W.argmax(1)
    top2 = np.partition(W, -2, axis=1)[:, -2:]
    df["w_max"] = top2[:, 1]
    df["w_margin"] = top2[:, 1] - top2[:, 0]
    return df, W


def tier_tables(df, K, species):
    """Per-macrogene gene counts by species and n_species tiers (own universe)."""
    cnt = pd.crosstab(df["mg"], df["species"]).reindex(index=range(K), columns=species, fill_value=0)
    nsp = (cnt >= 1).sum(axis=1).values
    nsp3 = (cnt >= 3).sum(axis=1).values
    return cnt, nsp, nsp3


# ════════════════════════════════════════════════════════════
# Partition helpers
# ════════════════════════════════════════════════════════════
def contingency(lx, ly, Kx, Ky):
    return sparse.coo_matrix((np.ones(len(lx), dtype=np.int32), (lx, ly)),
                             shape=(Kx, Ky)).toarray()


def agreement(lx, ly):
    return {"ARI": adjusted_rand_score(lx, ly),
            "AMI": adjusted_mutual_info_score(lx, ly),
            "NMI": normalized_mutual_info_score(lx, ly)}


def bin_status(j):
    return np.where(j >= 0.5, "stable", np.where(j >= 0.2, "partial", "reorganized"))


def match_table(C, row_label, col_label):
    """Best-match partner of every non-empty row macrogene."""
    size_r, size_c = C.sum(1), C.sum(0)
    rows = np.where(size_r > 0)[0]
    Cr = C[rows]
    best = Cr.argmax(1)
    ov = Cr[np.arange(len(rows)), best]
    p = Cr / size_r[rows, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        H = -np.nansum(np.where(p > 0, p * np.log(p), 0.0), axis=1)
    jac = ov / (size_r[rows] + size_c[best] - ov)
    return pd.DataFrame({
        "mg_" + row_label: rows,
        "n_genes": size_r[rows],
        "best_mg_" + col_label: best,
        "best_mg_n_genes": size_c[best],
        "overlap": ov,
        "frac_retained": ov / size_r[rows],
        "jaccard": jac,
        "n_partners": (Cr > 0).sum(1),
        "eff_n_partners": np.exp(H),
        "status": bin_status(jac),
    })


def hungarian(C):
    size_r, size_c = C.sum(1), C.sum(0)
    r_idx, c_idx = np.where(size_r > 0)[0], np.where(size_c > 0)[0]
    sub = C[np.ix_(r_idx, c_idx)]
    r, c = linear_sum_assignment(sub, maximize=True)
    ov = sub[r, c]
    keep = ov > 0
    r, c, ov = r[keep], c[keep], ov[keep]
    out = pd.DataFrame({"mg_a": r_idx[r], "mg_b": c_idx[c], "overlap": ov,
                        "n_a": size_r[r_idx[r]], "n_b": size_c[c_idx[c]]})
    out["jaccard"] = out["overlap"] / (out["n_a"] + out["n_b"] - out["overlap"])
    return out.sort_values("overlap", ascending=False), ov.sum() / C.sum()


# ════════════════════════════════════════════════════════════
# 0. Load + diagnostics
# ════════════════════════════════════════════════════════════
t0 = time.time()
log("=" * 70)
log("SATURN macrogene comparison: " + LA + " vs " + LB)
log("=" * 70)
A, WA = load_g2m(args.run_a, LA)
B, WB = load_g2m(args.run_b, LB)
if WA.shape[1] != WB.shape[1]:
    sys.exit("Different numbers of macrogenes: " + str(WA.shape[1]) + " vs " + str(WB.shape[1]))
K = WA.shape[1]
for lab, df in [(LA, A), (LB, B)]:
    extra = sorted(set(df["species"]) - set(SPECIES))
    if extra:
        log("WARNING: " + lab + " has species not in --species: " + str(extra))

log("\n── 0. Diagnostics ──")
diag = []
for lab, df, W in [(LA, A, WA), (LB, B, WB)]:
    sizes = df["mg"].value_counts().reindex(range(K), fill_value=0).values
    nz = sizes[sizes > 0]
    diag.append({
        "run": lab, "n_genes": len(df), "n_macrogenes": K,
        **{"n_" + s: int((df["species"] == s).sum()) for s in SPECIES},
        "W_min": float(W.min()), "W_max": float(W.max()), "W_mean": float(W.mean()),
        "W_frac_zero": float((W == 0).mean()),
        "W_rowsum_mean": float(W.sum(1).mean()), "W_rowsum_sd": float(W.sum(1).std()),
        "w_max_median": float(df["w_max"].median()),
        "w_margin_median": float(df["w_margin"].median()),
        "n_nonempty_mg": int((sizes > 0).sum()),
        "mg_size_median": float(np.median(nz)), "mg_size_p90": float(np.percentile(nz, 90)),
        "mg_size_max": int(nz.max()), "n_singleton_mg": int((nz == 1).sum()),
    })
diag = pd.DataFrame(diag)
save(diag, "00_diagnostics_runs.csv")
log(diag.T.to_string(header=False))

# Gene universes
ka, kb = set(A["key"]), set(B["key"])
shared = sorted(ka & kb)
uni = []
for s in SPECIES:
    a = {k for k in ka if k.startswith(s + "_")}
    b = {k for k in kb if k.startswith(s + "_")}
    uni.append({"species": s, "n_" + LA: len(a), "n_" + LB: len(b), "n_shared": len(a & b),
                "only_" + LA: len(a - b), "only_" + LB: len(b - a)})
uni = pd.DataFrame(uni)
save(uni, "00_gene_universe_by_species.csv")
log("\nGene universe (HV genes in pkl):")
log(uni.to_string(index=False))
only = pd.concat([A.loc[~A["key"].isin(kb), ["species", "gene", "mg"]].assign(only_in=LA),
                  B.loc[~B["key"].isin(ka), ["species", "gene", "mg"]].assign(only_in=LB)])
save(only, "00_genes_only_in_one_run.csv")

# Aligned shared-gene views
SA = A.set_index("key").loc[shared].reset_index()
SB = B.set_index("key").loc[shared].reset_index()
ia = A.reset_index().set_index("key").loc[shared, "index"].values
ib = B.reset_index().set_index("key").loc[shared, "index"].values
la, lb = SA["mg"].values, SB["mg"].values
spv = SA["species"].values
N = len(shared)
log("\nShared genes used for A/B comparisons: " + str(N))

# Own-universe tier tables
cntA, nspA, nsp3A = tier_tables(A, K, SPECIES)
cntB, nspB, nsp3B = tier_tables(B, K, SPECIES)

# ════════════════════════════════════════════════════════════
# A. Partition agreement
# ════════════════════════════════════════════════════════════
log("\n── A1. Global agreement (shared genes) ──")
rows = []
for scope in ["ALL"] + SPECIES:
    m = np.ones(N, bool) if scope == "ALL" else (spv == scope)
    r = agreement(la[m], lb[m])
    rows.append({"scope": scope, "n_genes": int(m.sum()),
                 "n_mg_" + LA: len(np.unique(la[m])), "n_mg_" + LB: len(np.unique(lb[m])), **r})
agree = pd.DataFrame(rows)
save(agree, "A1_agreement_ARI_AMI_NMI.csv")
log(agree.round(3).to_string(index=False))

log("\n── A2. Macrogene matching ──")
C = contingency(la, lb, K, K)
mAB = match_table(C, LA, LB)
mBA = match_table(C.T, LB, LA)
for df_, nsp_own, nsp_other, lab_own, lab_other in [(mAB, nspA, nspB, LA, LB), (mBA, nspB, nspA, LB, LA)]:
    own_col, best_col = "mg_" + lab_own, "best_mg_" + lab_other
    df_["n_species_" + lab_own] = nsp_own[df_[own_col].values]
    df_["n_species_best_" + lab_other] = nsp_other[df_[best_col].values]
# per-species composition of each A macrogene (shared universe)
compA = pd.crosstab(la, spv).reindex(columns=SPECIES, fill_value=0).add_prefix("n_")
mAB = mAB.merge(compA, left_on="mg_" + LA, right_index=True, how="left")
compB = pd.crosstab(lb, spv).reindex(columns=SPECIES, fill_value=0).add_prefix("n_")
mBA = mBA.merge(compB, left_on="mg_" + LB, right_index=True, how="left")
save(mAB, "A2_bestmatch_" + LA + "_to_" + LB + ".csv")
save(mBA, "A2_bestmatch_" + LB + "_to_" + LA + ".csv")

hung, frac_matched = hungarian(C)
hung = hung.rename(columns={"mg_a": "mg_" + LA, "mg_b": "mg_" + LB, "n_a": "n_" + LA, "n_b": "n_" + LB})
save(hung, "A2_hungarian_1to1_matching.csv")

for lab, df_ in [(LA + "->" + LB, mAB), (LB + "->" + LA, mBA)]:
    st = df_["status"].value_counts().reindex(["stable", "partial", "reorganized"], fill_value=0)
    log(lab + ": " + str(len(df_)) + " non-empty MGs | median Jaccard " +
        str(round(df_["jaccard"].median(), 3)) + " | median frac retained " +
        str(round(df_["frac_retained"].median(), 3)) + " | stable/partial/reorganized = " +
        "/".join(str(x) for x in st.values))
    big = df_[df_["n_genes"] >= 5]
    log("   (MGs with >=5 genes: n=" + str(len(big)) + ", median Jaccard " +
        str(round(big["jaccard"].median(), 3)) + ")")
log("Hungarian 1:1 matching places " + str(round(100 * frac_matched, 1)) +
    "% of shared genes on matched macrogene pairs")
log("Status bins are descriptive only: stable J>=0.5, partial 0.2-0.5, reorganized <0.2")

log("\n── A3. Per-gene co-member stability ──")
size_a, size_b = C.sum(1), C.sum(0)
c = C[la, lb]
sa, sb = size_a[la], size_b[lb]
inter, union = c - 1, sa + sb - c - 1
with np.errstate(divide="ignore", invalid="ignore"):
    jac = np.where(union > 0, inter / union, np.nan)

# cross-species co-members: genes of other species sharing the macrogene
xa = np.zeros(N); xb = np.zeros(N); xi = np.zeros(N)
for s in SPECIES:
    m = spv == s
    Cs = contingency(la[m], lb[m], K, K)
    sa_s, sb_s = Cs.sum(1), Cs.sum(0)
    xa[m] = sa[m] - sa_s[la[m]]
    xb[m] = sb[m] - sb_s[lb[m]]
    xi[m] = c[m] - Cs[la[m], lb[m]]
xu = xa + xb - xi
with np.errstate(divide="ignore", invalid="ignore"):
    xjac = np.where(xu > 0, xi / xu, np.nan)

gene_tab = pd.DataFrame({
    "species": spv, "gene": SA["gene"].values,
    "mg_" + LA: la, "mg_" + LB: lb,
    "mg_size_" + LA: sa, "mg_size_" + LB: sb, "n_shared_comembers": inter,
    "comember_jaccard": jac,
    "n_xsp_comembers_" + LA: xa.astype(int), "n_xsp_comembers_" + LB: xb.astype(int),
    "n_xsp_shared": xi.astype(int), "xsp_comember_jaccard": xjac,
    "tier_" + LA: nspA[la], "tier_" + LB: nspB[lb],
    "w_max_" + LA: SA["w_max"].values, "w_max_" + LB: SB["w_max"].values,
    "w_margin_" + LA: SA["w_margin"].values, "w_margin_" + LB: SB["w_margin"].values,
})

# ── A4. soft-weight kNN overlap ──
def knn(W, sp_codes, k, chunk=1024):
    Wn = W / np.clip(np.linalg.norm(W, axis=1, keepdims=True), 1e-12, None)
    n = Wn.shape[0]
    nn_all = np.empty((n, k), np.int32)
    nn_x = np.empty((n, k), np.int32)
    for st in range(0, n, chunk):
        en = min(st + chunk, n)
        S = Wn[st:en] @ Wn.T
        S[np.arange(en - st), np.arange(st, en)] = -np.inf
        nn_all[st:en] = np.argpartition(-S, k, axis=1)[:, :k]
        S[sp_codes[st:en, None] == sp_codes[None, :]] = -np.inf
        nn_x[st:en] = np.argpartition(-S, k, axis=1)[:, :k]
    return nn_all, nn_x


def overlap_frac(x, y):
    return (x[:, :, None] == y[:, None, :]).any(2).sum(1) / x.shape[1]


if not args.skip_knn:
    log("\n── A4. Soft-weight kNN overlap (cosine on macrogene weight vectors, k=" + str(args.knn_k) + ") ──")
    sp_codes = pd.Categorical(spv, categories=SPECIES).codes
    t1 = time.time()
    nnA, nxA = knn(WA[ia], sp_codes, args.knn_k)
    nnB, nxB = knn(WB[ib], sp_codes, args.knn_k)
    gene_tab["knn_overlap_all"] = overlap_frac(nnA, nnB)
    gene_tab["knn_overlap_xsp"] = overlap_frac(nxA, nxB)
    log("  kNN done in " + str(round(time.time() - t1)) + " s")
    # random expectation for reference
    log("  random-neighbour expectation ~ k/N = " + str(round(args.knn_k / N, 5)))

save(gene_tab, "A3_per_gene_stability.csv")

agg = {"n_genes": ("gene", "size"),
       "comember_J_median": ("comember_jaccard", "median"),
       "frac_comembers_identical": ("comember_jaccard", lambda v: float((v == 1).mean())),
       "frac_comembers_all_new": ("comember_jaccard", lambda v: float((v == 0).mean())),
       "frac_singleton_both": ("comember_jaccard", lambda v: float(v.isna().mean())),
       "xsp_J_median": ("xsp_comember_jaccard", "median"),
       "frac_xsp_partners_" + LA: ("n_xsp_comembers_" + LA, lambda v: float((v > 0).mean())),
       "frac_xsp_partners_" + LB: ("n_xsp_comembers_" + LB, lambda v: float((v > 0).mean()))}
if not args.skip_knn:
    agg["knn_all_median"] = ("knn_overlap_all", "median")
    agg["knn_xsp_median"] = ("knn_overlap_xsp", "median")
per_sp = gene_tab.groupby("species").agg(**agg).reindex(SPECIES)
allrow = gene_tab.assign(species="ALL").groupby("species").agg(**agg)
per_sp = pd.concat([per_sp, allrow])
save(per_sp, "A3_per_gene_stability_by_species.csv", index=True)
log(per_sp.round(3).to_string())

# Does stability depend on assignment confidence? (margin quartiles in A)
q = pd.qcut(gene_tab["w_margin_" + LA], 4, labels=["Q1_low", "Q2", "Q3", "Q4_high"], duplicates="drop")
conf = gene_tab.groupby(q, observed=True)["comember_jaccard"].median()
log("\nMedian co-member Jaccard by " + LA + " assignment margin quartile:")
log(conf.round(3).to_string())

# ════════════════════════════════════════════════════════════
# B. Species structure
# ════════════════════════════════════════════════════════════
log("\n── B1. Macrogene species-tier distribution (each run's own gene universe) ──")
nsp_max = len(SPECIES)
tiers = pd.DataFrame({"n_species": range(1, nsp_max + 1)})
for lab, nsp, nsp3, df in [(LA, nspA, nsp3A, A), (LB, nspB, nsp3B, B)]:
    used = np.bincount(df["mg"], minlength=K) > 0
    tiers["n_mg_" + lab] = [int(((nsp == t) & used).sum()) for t in tiers["n_species"]]
    tiers["n_mg_min3_" + lab] = [int(((nsp3 == t) & used).sum()) for t in tiers["n_species"]]
save(tiers, "B1_tier_distribution.csv")
log(tiers.to_string(index=False))
log("(min3 = species counted only if it contributes >=3 genes)")

gene_tiers = []
for lab, df, nsp in [(LA, A, nspA), (LB, B, nspB)]:
    t = pd.crosstab(df["species"], nsp[df["mg"].values]).reindex(SPECIES).fillna(0).astype(int)
    t.columns = ["genes_in_tier" + str(x) for x in t.columns]
    gene_tiers.append(t.assign(run=lab).reset_index())
gene_tiers = pd.concat(gene_tiers)
save(gene_tiers, "B1_genes_per_tier_by_species.csv")
log("\nGenes per tier, by species:")
log(gene_tiers.to_string(index=False))

log("\n── B2. Tier transitions per gene (shared genes; rows " + LA + ", cols " + LB + ") ──")
trans = []
for s in SPECIES + ["ALL"]:
    m = np.ones(N, bool) if s == "ALL" else (spv == s)
    ct = pd.crosstab(pd.Series(nspA[la[m]], name="tier_" + LA),
                     pd.Series(nspB[lb[m]], name="tier_" + LB))
    ct = ct.reindex(index=range(1, nsp_max + 1), columns=range(1, nsp_max + 1), fill_value=0)
    log(s + ":\n" + ct.to_string())
    long = ct.stack().rename("n_genes").reset_index()
    trans.append(long.assign(species=s))
save(pd.concat(trans), "B2_tier_transitions.csv")

log("\n── B3. Species-pair shared macrogenes (own universe) ──")
pairs = [(SPECIES[i], SPECIES[j]) for i in range(nsp_max) for j in range(i + 1, nsp_max)]
pr = []
for s, t in pairs:
    row = {"species_1": s, "species_2": t}
    for lab, cnt in [(LA, cntA), (LB, cntB)]:
        row["shared_mg_" + lab] = int(((cnt[s] > 0) & (cnt[t] > 0)).sum())
    pr.append(row)

log("\n── B4. Cross-species co-assigned gene pairs (shared genes) ──")


def pair_codes(labels, s, t):
    i_s, i_t = np.where(spv == s)[0], np.where(spv == t)[0]
    a = pd.DataFrame({"i": i_s, "mg": labels[i_s]})
    b = pd.DataFrame({"j": i_t, "mg": labels[i_t]})
    m = a.merge(b, on="mg")
    return m["i"].values.astype(np.int64) * N + m["j"].values.astype(np.int64)


for row in pr:
    s, t = row["species_1"], row["species_2"]
    pa, pb = pair_codes(la, s, t), pair_codes(lb, s, t)
    both = np.intersect1d(pa, pb, assume_unique=True).size
    row.update({"gene_pairs_" + LA: pa.size, "gene_pairs_" + LB: pb.size, "gene_pairs_both": both,
                "pair_jaccard": both / max(pa.size + pb.size - both, 1),
                "frac_" + LA + "_pairs_kept": both / max(pa.size, 1)})
pr = pd.DataFrame(pr)
save(pr, "B3_B4_species_pairs.csv")
log(pr.round(3).to_string(index=False))

# ════════════════════════════════════════════════════════════
# R. Reference anchor: canonical 5-species ESM-1b run
# ════════════════════════════════════════════════════════════
ref_ok = len(find_pkl(args.ref_run)) == 1
if not ref_ok:
    log("\n── R. Reference anchor skipped (no unique pkl in " + args.ref_run + ") ──")
else:
    log("\n── R. Reference anchor: " + args.ref_label + " (non-" + "/".join(SPECIES) + " species dropped) ──")
    Rfull, _ = load_g2m(args.ref_run, args.ref_label)
    KR = int(Rfull["mg"].max()) + 1
    ref_sp_all = sorted(Rfull["species"].unique())
    cntR, nspR, _ = tier_tables(Rfull, max(KR, 1), ref_sp_all)   # tier incl. all ref species (e.g. Spur)
    R = Rfull[Rfull["species"].isin(SPECIES)].set_index("key")
    common = sorted(set(R.index) & set(shared))
    pos = pd.Series(np.arange(N), index=shared).loc[common].values
    lr = R.loc[common, "mg"].values
    lra, lrb, spc = la[pos], lb[pos], spv[pos]
    log("Genes common to ref and both runs: " + str(len(common)))

    rr = []
    for scope in ["ALL"] + SPECIES:
        m = np.ones(len(common), bool) if scope == "ALL" else (spc == scope)
        rr.append({"scope": scope, "n_genes": int(m.sum()),
                   "ARI_ref_vs_" + LA: adjusted_rand_score(lr[m], lra[m]),
                   "ARI_ref_vs_" + LB: adjusted_rand_score(lr[m], lrb[m]),
                   "ARI_" + LA + "_vs_" + LB: adjusted_rand_score(lra[m], lrb[m]),
                   "AMI_ref_vs_" + LA: adjusted_mutual_info_score(lr[m], lra[m]),
                   "AMI_ref_vs_" + LB: adjusted_mutual_info_score(lr[m], lrb[m]),
                   "AMI_" + LA + "_vs_" + LB: adjusted_mutual_info_score(lra[m], lrb[m])})
    rr = pd.DataFrame(rr)
    save(rr, "R1_reference_agreement.csv")
    log(rr.round(3).to_string(index=False))

    KX = max(KR, K)
    mRA = match_table(contingency(lr, lra, KX, K), "ref", LA)
    mRB = match_table(contingency(lr, lrb, KX, K), "ref", LB)
    def sfx(lab):
        return lambda x: x if x in ("mg_ref", "n_genes", "best_mg_" + lab) else x + "_" + lab
    ra = mRA.rename(columns=sfx(LA))
    rb = mRB.drop(columns="n_genes").rename(columns=sfx(LB))
    refm = ra.merge(rb, on="mg_ref")
    refm.insert(1, "n_species_ref", nspR[refm["mg_ref"].values])
    save(refm, "R2_reference_mg_retention.csv")
    log("\nRetention of canonical macrogenes (median best-match Jaccard), by canonical tier:")
    summ = refm.groupby("n_species_ref").agg(
        n_mg=("mg_ref", "size"),
        **{"J_" + LA: ("jaccard_" + LA, "median"), "J_" + LB: ("jaccard_" + LB, "median"),
           "stable_" + LA: ("status_" + LA, lambda v: int((v == "stable").sum())),
           "stable_" + LB: ("status_" + LB, lambda v: int((v == "stable").sum()))})
    log(summ.round(3).to_string())

    focus = [int(x) for x in args.focus_ref_mgs.split(",") if x.strip()]
    ft = []
    for mg in focus:
        m = lr == mg
        if not m.any():
            log("Focus MG" + str(mg) + ": no genes in common set")
            continue
        bA = pd.Series(lra[m]).mode().iloc[0]
        bB = pd.Series(lrb[m]).mode().iloc[0]
        sub = pd.DataFrame({"ref_mg": mg, "species": spc[m],
                            "gene": [k.split("_", 1)[1] for k in np.array(common)[m]],
                            "mg_" + LA: lra[m], "mg_" + LB: lrb[m]})
        sub["with_majority_" + LA] = sub["mg_" + LA] == bA
        sub["with_majority_" + LB] = sub["mg_" + LB] == bB
        ft.append(sub)
        log("Focus MG" + str(mg) + " (" + str(m.sum()) + " genes, ref tier " + str(nspR[mg]) + "): " +
            LA + " keeps " + str(int(sub["with_majority_" + LA].sum())) + " together in MG" + str(bA) +
            "; " + LB + " keeps " + str(int(sub["with_majority_" + LB].sum())) + " together in MG" + str(bB))
    if ft:
        save(pd.concat(ft), "R3_focus_macrogenes_gene_trace.csv")

# ════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK,
                     "ytick.color": INK, "axes.titleweight": "bold", "axes.titlesize": 11})

# Fig 1: tier distribution
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
x = np.arange(nsp_max)
for ax, pref, ttl in [(axes[0], "n_mg_", "Species per macrogene (>=1 gene)"),
                      (axes[1], "n_mg_min3_", "Species per macrogene (>=3 genes)")]:
    ax.bar(x - 0.2, tiers[pref + LA], 0.38, color=COL_A, label=LA)
    ax.bar(x + 0.2, tiers[pref + LB], 0.38, color=COL_B, label=LB)
    ax.set_xticks(x, [str(t) + "/" + str(nsp_max) for t in tiers["n_species"]])
    ax.set_title(ttl); ax.set_ylabel("Macrogenes")
axes[0].legend(frameon=False)
fig.tight_layout(); fig.savefig(OUT / "figures" / "Fig1_tier_distribution.png", dpi=200); plt.close(fig)

# Fig 2: best-match Jaccard distributions
fig, ax = plt.subplots(figsize=(5.5, 3.6))
bins = np.linspace(0, 1, 26)
ax.hist(mAB["jaccard"], bins=bins, histtype="step", lw=2, color=COL_A, label=LA + " MG -> best " + LB)
ax.hist(mBA["jaccard"], bins=bins, histtype="step", lw=2, color=COL_B, label=LB + " MG -> best " + LA)
for v in (0.2, 0.5):
    ax.axvline(v, color=MUTED, lw=1, ls=":")
ax.set_xlabel("Best-match Jaccard"); ax.set_ylabel("Macrogenes")
ax.set_title("Macrogene best-match overlap"); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "figures" / "Fig2_bestmatch_jaccard.png", dpi=200); plt.close(fig)

# Fig 3: per-gene stability by species
metrics = [("comember_jaccard", "Co-member Jaccard (all)"),
           ("xsp_comember_jaccard", "Co-member Jaccard (cross-species)")]
if not args.skip_knn:
    metrics += [("knn_overlap_all", "Soft-weight kNN overlap (all)"),
                ("knn_overlap_xsp", "Soft-weight kNN overlap (cross-species)")]
fig, axes = plt.subplots(1, len(metrics), figsize=(3.2 * len(metrics), 3.6), sharey=True)
for ax, (col, ttl) in zip(np.atleast_1d(axes), metrics):
    data = [gene_tab.loc[gene_tab["species"] == s, col].dropna().values for s in SPECIES]
    bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                    medianprops=dict(color=INK, lw=2))
    for b in bp["boxes"]:
        b.set_facecolor("#cde2fb"); b.set_edgecolor(COL_A)
    ax.set_xticks(range(1, len(SPECIES) + 1), SPECIES)
    ax.set_title(ttl, fontsize=9); ax.set_ylim(-0.02, 1.02)
np.atleast_1d(axes)[0].set_ylabel(LA + " vs " + LB + " agreement")
fig.tight_layout(); fig.savefig(OUT / "figures" / "Fig3_per_gene_stability.png", dpi=200); plt.close(fig)

# Fig 4: reference retention
if ref_ok:
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    ax.hist(refm["jaccard_" + LA], bins=bins, histtype="step", lw=2, color=COL_A, label=LA + " 4sp")
    ax.hist(refm["jaccard_" + LB], bins=bins, histtype="step", lw=2, color=COL_B, label=LB + " 4sp")
    ax.set_xlabel("Best-match Jaccard to canonical macrogene"); ax.set_ylabel("Canonical macrogenes")
    ax.set_title("Retention of canonical 5sp macrogenes"); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "figures" / "Fig4_reference_retention.png", dpi=200); plt.close(fig)

log("\nDone in " + str(round(time.time() - t0)) + " s. Outputs: " + str(OUT))
(OUT / "summary.txt").write_text("\n".join(SUMMARY) + "\n")
