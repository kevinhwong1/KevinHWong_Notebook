#!/usr/bin/env python3
"""
05d_embedding_species_structure.py

Pre-SATURN check of the raw protein embeddings (ESM-1b vs ESMC-600M) that fed
the two 4-species macrogene-only runs. Asks whether ESMC's stronger species
separation in the macrogenes (esp. Drer) is already present in the raw
embeddings, and whether any species looks like a conversion artifact.

Gene set: HV genes present in BOTH runs' genes_to_macrogenes.pkl (the same
31,808 genes compared in 05c). Embedding files are read from each run's own
in_data.csv, so these are exactly the vectors SATURN used.

Two views of each model's embeddings:
  raw      : cosine on the vectors as stored
  centered : cosine after subtracting the global mean across ALL genes and
             species (removes the shared ESM direction; NOT species-centered,
             so species signal is kept)

Checks (per model, per species)
  E0  integrity   : missing genes, exact-duplicate vectors, L2-norm distribution
  E1  species sep : between-species share of variance (R^2 of species label),
                    mean within- vs between-species cosine, centroid cosine matrix
  E2  kNN mixing  : fraction of each gene's k nearest neighbours from OTHER
                    species (random expectation printed alongside)
  E3  cross-model : within-species kNN overlap ESM-1b vs ESMC for the same gene.
                    A species whose gene->protein bridge is wrong in one model
                    shows much lower agreement than the others.
Outputs: --out_dir (CSV tables, PNG figures, summary.txt)
"""
import argparse
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUNS = Path("/scratch/dark_genes/SATURN_Mnemi/04_saturn_runs")

ap = argparse.ArgumentParser()
ap.add_argument("--run_a", default=str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESM1b_macrogene"))
ap.add_argument("--run_b", default=str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESMC600M_macrogene"))
ap.add_argument("--label_a", default="ESM1b")
ap.add_argument("--label_b", default="ESMC600M")
ap.add_argument("--species", default="Mlei,Cgig,Crob,Drer")
ap.add_argument("--knn_k", type=int, default=10)
ap.add_argument("--out_dir", default=str(RUNS / "20261005_embedding_species_structure_ESM1b_vs_ESMC600M"))
args = ap.parse_args()

SPECIES = args.species.split(",")
LA, LB = args.label_a, args.label_b
LABELS = [LA, LB]
OUT = Path(args.out_dir)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
K = args.knn_k

SP_COL = {"Mlei": "#2a78d6", "Cgig": "#eb6834", "Crob": "#1baf7a", "Drer": "#eda100"}
COL = {LA: "#2a78d6", LB: "#eb6834"}
INK, MUTED = "#222222", "#8a8a85"
SUMMARY = []


def log(msg=""):
    print(msg, flush=True)
    SUMMARY.append(msg)


def save(df, name, index=False):
    df.to_csv(OUT / name, index=index)
    print("  wrote " + name, flush=True)


# ════════════════════════════════════════════════════════════
# Gene set + embedding paths
# ════════════════════════════════════════════════════════════
def pkl_keys(run_dir):
    hits = sorted((Path(run_dir) / "saturn_results").glob("*genes_to_macrogenes*.pkl"))
    if len(hits) != 1:
        raise RuntimeError("Expected one genes_to_macrogenes pkl in " + str(run_dir) + ", found " + str(hits))
    return set(pickle.load(open(hits[0], "rb")).keys())


def emb_paths(run_dir):
    d = pd.read_csv(Path(run_dir) / "in_data.csv")
    return dict(zip(d["species"], d["embedding_path"]))


def load_pt(path):
    import torch
    obj = torch.load(path, map_location="cpu")
    return obj


def vec(v):
    if hasattr(v, "detach"):
        v = v.detach().cpu().numpy()
    return np.asarray(v, dtype=np.float32).ravel()


t0 = time.time()
log("=" * 70)
log("Raw embedding species structure: " + LA + " vs " + LB)
log("=" * 70)
shared = sorted(pkl_keys(args.run_a) & pkl_keys(args.run_b))
shared = [k for k in shared if k.split("_", 1)[0] in SPECIES]
genes = pd.DataFrame({"key": shared,
                      "species": [k.split("_", 1)[0] for k in shared],
                      "gene": [k.split("_", 1)[1] for k in shared]})
log("Shared HV genes: " + str(len(genes)))
log(genes["species"].value_counts().reindex(SPECIES).to_string())

# ── Load embedding dicts and resolve every HV gene in both models ──
# Exact key match first, then a case-insensitive fallback (Drer mixes upper-
# and lower-case symbols). Genes unresolved in EITHER model are dropped from
# all checks so both models are compared on the identical gene set.
DICTS, RESOLVED, unresolved = {}, {}, []
for lab, run in [(LA, args.run_a), (LB, args.run_b)]:
    paths = emb_paths(run)
    for s in SPECIES:
        log(lab + " " + s + ": " + paths[s])
        d = load_pt(paths[s])
        DICTS[(lab, s)] = d
        lower = {}
        for k in d:
            lower.setdefault(k.lower(), []).append(k)
        n_exact = n_ci = 0
        for g in genes.loc[genes["species"] == s, "gene"]:
            if g in d:
                RESOLVED[(lab, s, g)] = g; n_exact += 1
            elif len(lower.get(g.lower(), [])) == 1:
                RESOLVED[(lab, s, g)] = lower[g.lower()][0]; n_ci += 1
            else:
                cands = lower.get(g.lower(), [])
                unresolved.append({"model": lab, "species": s, "gene": g,
                                   "reason": "ambiguous_case" if cands else "absent",
                                   "candidates": ";".join(cands)})
        log("  resolved exact=" + str(n_exact) + " case-insensitive=" + str(n_ci))

unres = pd.DataFrame(unresolved, columns=["model", "species", "gene", "reason", "candidates"])
save(unres, "E0_unresolved_genes.csv")
if len(unres):
    log("Unresolved HV genes (dropped from all checks, both models): " + str(len(unres)))
    log(unres.to_string(index=False))
    drop = set(zip(unres["species"], unres["gene"]))
    genes = genes[[(s, g) not in drop for s, g in zip(genes["species"], genes["gene"])]]
log("Genes analysed: " + str(len(genes)))

E = {}
integ = []
for lab in [LA, LB]:
    blocks = []
    for s in SPECIES:
        d = DICTS[(lab, s)]
        g = genes.loc[genes["species"] == s, "gene"].tolist()
        M = np.vstack([vec(d[RESOLVED[(lab, s, x)]]) for x in g])
        blocks.append(M)
        # exact duplicate vectors among this species' HV genes
        _, inv, cnt = np.unique(M.round(6), axis=0, return_inverse=True, return_counts=True)
        dup_genes = int((cnt[inv.ravel()] > 1).sum())
        nrm = np.linalg.norm(M, axis=1)
        integ.append({"model": lab, "species": s, "n_genes": len(g), "dim": M.shape[1],
                      "n_genes_in_file": len(d), "n_genes_sharing_identical_vector": dup_genes,
                      "n_distinct_vectors": int(len(cnt)),
                      "norm_mean": float(nrm.mean()), "norm_sd": float(nrm.std()),
                      "norm_min": float(nrm.min()), "norm_max": float(nrm.max())})
    E[lab] = np.vstack(blocks)
del DICTS

# row order of E[*] is species-blocked in SPECIES order; align the gene table
genes = pd.concat([genes[genes["species"] == s] for s in SPECIES]).reset_index(drop=True)
spv = genes["species"].values
spc = pd.Categorical(spv, categories=SPECIES).codes
N = len(genes)

log("\n── E0. Integrity ──")
integ = pd.DataFrame(integ)
save(integ, "E0_integrity.csv")
log(integ.round(3).to_string(index=False))


# ════════════════════════════════════════════════════════════
# Helpers
# ════════════════════════════════════════════════════════════
def unit(X):
    return X / np.clip(np.linalg.norm(X, axis=1, keepdims=True), 1e-12, None)


def views(X):
    return {"raw": unit(X), "centered": unit(X - X.mean(0, keepdims=True))}


def species_r2(X):
    """Between-species sum of squares / total sum of squares."""
    mu = X.mean(0)
    tot = ((X - mu) ** 2).sum()
    btw = sum((spc == i).sum() * ((X[spc == i].mean(0) - mu) ** 2).sum() for i in range(len(SPECIES)))
    return btw / tot


def knn_all(U, k, chunk=1024):
    n = U.shape[0]
    nn = np.empty((n, k), np.int32)
    for st in range(0, n, chunk):
        en = min(st + chunk, n)
        S = U[st:en] @ U.T
        S[np.arange(en - st), np.arange(st, en)] = -np.inf
        nn[st:en] = np.argpartition(-S, k, axis=1)[:, :k]
    return nn


def knn_within(U, k, chunk=1024):
    """Nearest neighbours restricted to the gene's own species (global indices)."""
    nn = np.empty((U.shape[0], k), np.int64)
    for i in range(len(SPECIES)):
        idx = np.where(spc == i)[0]
        sub = knn_all(U[idx], k, chunk)
        nn[idx] = idx[sub]
    return nn


def mean_cos_blocks(U):
    """Mean pairwise cosine within each species pair, via block centroids."""
    out = np.zeros((len(SPECIES), len(SPECIES)))
    sums = [U[spc == i].sum(0) for i in range(len(SPECIES))]
    ns = [(spc == i).sum() for i in range(len(SPECIES))]
    for i in range(len(SPECIES)):
        for j in range(len(SPECIES)):
            if i == j:   # exclude self-pairs (cos = 1)
                out[i, j] = (sums[i] @ sums[i] - ns[i]) / (ns[i] * (ns[i] - 1))
            else:
                out[i, j] = sums[i] @ sums[j] / (ns[i] * ns[j])
    return out


# ════════════════════════════════════════════════════════════
# E1 / E2 per model and view
# ════════════════════════════════════════════════════════════
log("\n── E1. Species separation ──")
sep_rows, cos_rows, knn_rows = [], [], []
nn_store, nnw_store = {}, {}
gene_tab = genes[["species", "gene"]].copy()
for lab in LABELS:
    for vname, U in views(E[lab]).items():
        r2 = species_r2(U)
        C = mean_cos_blocks(U)
        within = np.diag(C)
        betw = np.array([np.mean(np.delete(C[i], i)) for i in range(len(SPECIES))])
        for i, s in enumerate(SPECIES):
            sep_rows.append({"model": lab, "view": vname, "species": s,
                             "mean_cos_within": within[i], "mean_cos_to_other_species": betw[i],
                             "within_minus_between": within[i] - betw[i]})
        sep_rows.append({"model": lab, "view": vname, "species": "ALL_R2_species", "mean_cos_within": np.nan,
                         "mean_cos_to_other_species": np.nan, "within_minus_between": np.nan, "R2_species": r2})
        cm = pd.DataFrame(C, index=SPECIES, columns=SPECIES)
        cos_rows.append(cm.assign(model=lab, view=vname).reset_index().rename(columns={"index": "species"}))
        log(lab + " [" + vname + "]  R2(species) = " + str(round(r2, 3)))
        log("  mean pairwise cosine (rows/cols = species):\n" + cm.round(3).to_string())

        t1 = time.time()
        nn = knn_all(U, K)
        nn_store[(lab, vname)] = nn
        xfrac = (spc[nn] != spc[:, None]).mean(1)
        gene_tab["xsp_nn_frac_" + lab + "_" + vname] = xfrac
        for i, s in enumerate(SPECIES):
            m = spc == i
            knn_rows.append({"model": lab, "view": vname, "species": s,
                             "xsp_nn_frac_mean": xfrac[m].mean(),
                             "frac_genes_with_any_xsp_nn": (xfrac[m] > 0).mean(),
                             "random_expectation": (N - m.sum()) / (N - 1)})
        if vname == "centered":
            nnw_store[lab] = knn_within(U, K)
        print("  kNN " + lab + " " + vname + " " + str(round(time.time() - t1)) + " s", flush=True)

save(pd.DataFrame(sep_rows), "E1_species_separation.csv")
save(pd.concat(cos_rows), "E1_mean_cosine_species_matrix.csv")

log("\n── E2. Cross-species kNN mixing (k=" + str(K) + ") ──")
knn_df = pd.DataFrame(knn_rows)
save(knn_df, "E2_knn_cross_species_mixing.csv")
log(knn_df.round(3).to_string(index=False))

# Which species do each species' cross-species neighbours come from?
nbr_rows = []
for lab in LABELS:
    nn = nn_store[(lab, "centered")]
    for i, s in enumerate(SPECIES):
        nb = spc[nn[spc == i]].ravel()
        cnt = np.bincount(nb, minlength=len(SPECIES)) / nb.size
        nbr_rows.append({"model": lab, "query_species": s, **{"nn_" + t: cnt[j] for j, t in enumerate(SPECIES)}})
nbr = pd.DataFrame(nbr_rows)
save(nbr, "E2_neighbour_species_composition_centered.csv")
log("\nNeighbour species composition (centered view; rows sum to 1):")
log(nbr.round(3).to_string(index=False))

# ════════════════════════════════════════════════════════════
# E3. Cross-model agreement within species
# ════════════════════════════════════════════════════════════
log("\n── E3. Within-species kNN agreement, " + LA + " vs " + LB + " (centered, k=" + str(K) + ") ──")
a, b = nnw_store[LA], nnw_store[LB]
ov = (a[:, :, None] == b[:, None, :]).any(2).mean(1)
gene_tab["within_sp_knn_overlap"] = ov
e3 = []
for i, s in enumerate(SPECIES):
    m = spc == i
    e3.append({"species": s, "n_genes": int(m.sum()), "overlap_mean": ov[m].mean(),
               "overlap_median": np.median(ov[m]), "frac_zero_overlap": (ov[m] == 0).mean(),
               "random_expectation": K / (m.sum() - 1)})
e3 = pd.DataFrame(e3)
save(e3, "E3_cross_model_within_species_agreement.csv")
log(e3.round(3).to_string(index=False))
log("If one species is far below the others, suspect its gene->protein bridge in one model.")
save(gene_tab, "E_per_gene.csv")

# ════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.titleweight": "bold", "axes.titlesize": 11})

# Fig 1: cross-species kNN fraction by species, both models, both views
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
x = np.arange(len(SPECIES))
for ax, vname in zip(axes, ["raw", "centered"]):
    for off, lab in [(-0.2, LA), (0.2, LB)]:
        d = knn_df[(knn_df["model"] == lab) & (knn_df["view"] == vname)].set_index("species").loc[SPECIES]
        ax.bar(x + off, d["xsp_nn_frac_mean"], 0.38, color=COL[lab], label=lab)
    rnd = knn_df[(knn_df["model"] == LA) & (knn_df["view"] == vname)].set_index("species").loc[SPECIES, "random_expectation"]
    ax.scatter(x, rnd, marker="_", s=400, color=INK, lw=2, label="random", zorder=3)
    ax.set_xticks(x, SPECIES); ax.set_title("Cross-species neighbours (" + vname + ")")
axes[0].set_ylabel("Mean fraction of " + str(K) + "-NN from other species")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "figures" / "FigE1_xsp_knn_fraction.png", dpi=200); plt.close(fig)

# Fig 2: within-species cross-model agreement
fig, ax = plt.subplots(figsize=(4.5, 3.6))
data = [ov[spc == i] for i in range(len(SPECIES))]
bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=2))
for bx in bp["boxes"]:
    bx.set_facecolor("#cde2fb"); bx.set_edgecolor(COL[LA])
ax.set_xticks(range(1, len(SPECIES) + 1), SPECIES); ax.set_ylim(-0.02, 1.02)
ax.set_ylabel("Within-species " + str(K) + "-NN overlap")
ax.set_title(LA + " vs " + LB + " neighbour agreement", fontsize=10)
fig.tight_layout(); fig.savefig(OUT / "figures" / "FigE2_within_species_model_agreement.png", dpi=200); plt.close(fig)

# Fig 3: PCA (centered) coloured by species, one panel per model
fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))
rng = np.random.default_rng(42)
order = rng.permutation(N)
for ax, lab in zip(axes, LABELS):
    X = E[lab] - E[lab].mean(0)
    U_, S_, Vt = np.linalg.svd(X[rng.choice(N, min(N, 8000), replace=False)], full_matrices=False)
    P = X @ Vt[:2].T
    ve = (S_[:2] ** 2) / (S_ ** 2).sum()
    ax.scatter(P[order, 0], P[order, 1], s=2, c=[SP_COL.get(s, MUTED) for s in spv[order]], alpha=0.5, lw=0)
    ax.set_xlabel("PC1 (" + str(round(100 * ve[0], 1)) + "%)"); ax.set_ylabel("PC2 (" + str(round(100 * ve[1], 1)) + "%)")
    ax.set_title(lab)
handles = [plt.Line2D([], [], marker="o", ls="", color=SP_COL.get(s, MUTED), label=s) for s in SPECIES]
axes[1].legend(handles=handles, frameon=False, fontsize=8, loc="best")
fig.tight_layout(); fig.savefig(OUT / "figures" / "FigE3_PCA_by_species.png", dpi=200); plt.close(fig)

log("\nDone in " + str(round(time.time() - t0)) + " s. Outputs: " + str(OUT))
(OUT / "summary.txt").write_text("\n".join(SUMMARY) + "\n")
