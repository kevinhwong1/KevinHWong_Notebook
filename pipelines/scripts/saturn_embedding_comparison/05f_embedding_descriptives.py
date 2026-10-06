#!/usr/bin/env python3
"""
05f_embedding_descriptives.py

Descriptive statistics and figures for the protein-embedding sets fed to SATURN,
per model x species, written for sharing with collaborators.

Models (each skipped with a warning if its files are missing):
  ESM1b              03_embeddings/{SP}/{SP}_gene_embeddings.pt            (also used by the canonical 5sp run)
  ESMC600M_raw       03_embeddings_ESMC600M/{F}/{F}_gene_embeddings.pt     (FANTASIA v4, as converted by 05)
  ESMC600M_rescaled  03_embeddings_ESMC600M_rescaled/{F}/...               (05e: unit length x mean ESM-1b norm)
  ({F} = "Mnemi" for Mlei, otherwise the species code)

Gene sets
  all_in_file : every gene in that model's embedding file (differs by model)
  HV_shared   : the HV genes present in both 4sp macrogene runs (ESM1b + ESMC);
                all geometry / neighbour / ortholog checks use this set so the
                models are compared on identical genes

Sections
  D1 vector length (L2 norm) distribution
  D2 vector length vs protein length (Spearman) -- tests sum- vs mean-pooling
  D3 global geometry: anisotropy (mean pairwise cosine), effective dimensionality,
     share of variance explained by species (Euclidean = what SATURN's KMeans
     sees; unit = direction only)
  D4 cross-species neighbour mixing (k-NN; Euclidean and cosine)
  D5 ortholog recovery in embedding space: for a gene in species s, is its
     nearest neighbour in species t in the same eggNOG Metazoa OG? (hit@1, hit@10)

Outputs: --out_dir (CSV tables, PNG figures, summary.txt)
"""
import argparse
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path("/scratch/dark_genes/SATURN_Mnemi")
RUNS = BASE / "04_saturn_runs"

ap = argparse.ArgumentParser()
ap.add_argument("--hv_run_a", default=str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESM1b_macrogene"))
ap.add_argument("--hv_run_b", default=str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESMC600M_macrogene"))
ap.add_argument("--eggnog", default=str(BASE / "01_proteomes/eggnog_merged/ALL4_gene_protein_metazoaOG_slim.csv"))
ap.add_argument("--species", default="Mlei,Cgig,Crob,Drer")
ap.add_argument("--knn_k", type=int, default=10)
ap.add_argument("--out_dir", default=str(RUNS / "20261005_embedding_descriptives"))
args = ap.parse_args()

SPECIES = args.species.split(",")
FOLDER = {"Mlei": "Mnemi"}
FASTA = {"Mlei": BASE / "01_proteomes/Mlei/Mlei.protein.faa",
         "Cgig": BASE / "01_proteomes/Cgig/Cgig.protein.faa",
         "Crob": BASE / "01_proteomes/Crob/Crob_KY21_ghost.faa",
         "Drer": BASE / "01_proteomes/Drer/Drer.protein.faa"}
MODELS = {
    "ESM1b": lambda sp: BASE / "03_embeddings" / sp / f"{sp}_gene_embeddings.pt",
    "ESMC600M_raw": lambda sp: BASE / "03_embeddings_ESMC600M" / FOLDER.get(sp, sp) / f"{FOLDER.get(sp, sp)}_gene_embeddings.pt",
    "ESMC600M_rescaled": lambda sp: BASE / "03_embeddings_ESMC600M_rescaled" / FOLDER.get(sp, sp) / f"{FOLDER.get(sp, sp)}_gene_embeddings.pt",
}
# categorical palette: models keep the colours used for their SATURN runs in 05g
MCOL = {"ESM1b": "#eb6834", "ESMC600M_raw": "#1baf7a", "ESMC600M_rescaled": "#eda100"}
SCOL = dict(zip(SPECIES, ["#e87ba4", "#008300", "#4a3aa7", "#e34948"]))
INK, MUTED = "#222222", "#8a8a85"
K = args.knn_k
OUT = Path(args.out_dir)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
SUMMARY = []


def log(msg=""):
    print(msg, flush=True)
    SUMMARY.append(str(msg))


def save(df, name, index=False):
    df.to_csv(OUT / name, index=index)
    print("  wrote " + name, flush=True)


def vec(v):
    if hasattr(v, "detach"):
        v = v.detach().cpu().numpy()
    return np.asarray(v, dtype=np.float32).ravel()


def load_pt(p):
    import torch
    return torch.load(p, map_location="cpu")


def pkl_keys(run_dir):
    hits = sorted((Path(run_dir) / "saturn_results").glob("*genes_to_macrogenes*.pkl"))
    if len(hits) != 1:
        raise RuntimeError("Expected one genes_to_macrogenes pkl in " + str(run_dir))
    return set(pickle.load(open(hits[0], "rb")).keys())


def resolve(d, genes):
    """Exact key, then unique case-insensitive match (SATURN-equivalent)."""
    lower = {}
    for k in d:
        lower.setdefault(k.lower(), []).append(k)
    out = {}
    for g in genes:
        if g in d:
            out[g] = g
        elif len(lower.get(g.lower(), [])) == 1:
            out[g] = lower[g.lower()][0]
    return out


def fasta_lengths(path):
    lens, name, n = {}, None, 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens


t0 = time.time()
log("=" * 70)
log("Embedding descriptives: " + ", ".join(MODELS))
log("=" * 70)

# ── HV shared gene set ──
shared = sorted(pkl_keys(args.hv_run_a) & pkl_keys(args.hv_run_b))
hv = pd.DataFrame({"key": shared, "species": [k.split("_", 1)[0] for k in shared],
                   "gene": [k.split("_", 1)[1] for k in shared]})
hv = hv[hv["species"].isin(SPECIES)]
hv = pd.concat([hv[hv["species"] == s] for s in SPECIES]).reset_index(drop=True)
log("HV_shared genes: " + str(len(hv)) + "  " + str(hv["species"].value_counts().reindex(SPECIES).to_dict()))

# ── eggNOG: gene -> Metazoa OG(s), gene -> protein ids ──
egg = None
if Path(args.eggnog).exists():
    egg = pd.read_csv(args.eggnog)
    egg = egg[egg["species"].isin(SPECIES)]
    log("eggNOG slim table: " + args.eggnog + " (" + str(len(egg)) + " rows)")
else:
    log("WARNING: eggNOG table not found (" + args.eggnog + "); D2 protein lengths and D5 skipped")

# ── Load embeddings ──
EMB_HV, ALLNORM, avail = {}, [], []
for m, pathf in MODELS.items():
    paths = {s: pathf(s) for s in SPECIES}
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        log("SKIP model " + m + " (missing: " + "; ".join(missing) + ")")
        continue
    avail.append(m)
    blocks, keep = [], np.ones(len(hv), bool)
    for s in SPECIES:
        d = load_pt(paths[s])
        A = np.vstack([vec(v) for v in d.values()])
        ALLNORM.append(pd.DataFrame({"model": m, "species": s, "gene": list(d.keys()),
                                     "norm": np.linalg.norm(A, axis=1), "dim": A.shape[1]}))
        g = hv.loc[hv["species"] == s, "gene"].tolist()
        r = resolve(d, g)
        blocks.append(np.vstack([vec(d[r[x]]) if x in r else np.full(A.shape[1], np.nan, np.float32) for x in g]))
        log(m + " " + s + ": " + str(len(d)) + " genes in file, dim " + str(A.shape[1]) +
            ", HV resolved " + str(len(r)) + "/" + str(len(g)))
        del d, A
    EMB_HV[m] = np.vstack(blocks)

# genes missing in any model are dropped everywhere
ok = np.ones(len(hv), bool)
for m in avail:
    ok &= ~np.isnan(EMB_HV[m]).any(1)
if (~ok).sum():
    log("Dropping " + str(int((~ok).sum())) + " HV genes unresolved in at least one model")
hv = hv[ok].reset_index(drop=True)
for m in avail:
    EMB_HV[m] = EMB_HV[m][ok]
spv = hv["species"].values
spc = pd.Categorical(spv, categories=SPECIES).codes
N = len(hv)
allnorm = pd.concat(ALLNORM, ignore_index=True)

# ════════════════════════════════════════════════════════════
# D1. Vector length
# ════════════════════════════════════════════════════════════
log("\n── D1. Vector length (L2 norm) ──")


def norm_stats(x):
    return {"n": len(x), "mean": x.mean(), "sd": x.std(), "CV": x.std() / x.mean(),
            "median": np.median(x), "q05": np.quantile(x, 0.05), "q95": np.quantile(x, 0.95),
            "min": x.min(), "max": x.max(), "fold_q95_q05": np.quantile(x, 0.95) / np.quantile(x, 0.05)}


rows = []
hv_norm = {}
for m in avail:
    nh = np.linalg.norm(EMB_HV[m], axis=1)
    hv_norm[m] = nh
    for s in SPECIES:
        a = allnorm[(allnorm["model"] == m) & (allnorm["species"] == s)]
        rows.append({"model": m, "species": s, "gene_set": "all_in_file",
                     "dim": int(a["dim"].iloc[0]), **norm_stats(a["norm"].values)})
        rows.append({"model": m, "species": s, "gene_set": "HV_shared",
                     "dim": int(a["dim"].iloc[0]), **norm_stats(nh[spc == SPECIES.index(s)])})
d1 = pd.DataFrame(rows)
save(d1, "D1_vector_length_stats.csv")
log(d1[d1["gene_set"] == "HV_shared"].drop(columns="gene_set").round(3).to_string(index=False))

# ════════════════════════════════════════════════════════════
# D2. Vector length vs protein length
# ════════════════════════════════════════════════════════════
d2 = None
if egg is not None:
    log("\n── D2. Vector length vs protein length (gene length = mean of its proteins in the eggNOG table) ──")
    glen = []
    for s in SPECIES:
        if not FASTA[s].exists():
            log("  " + s + ": FASTA missing (" + str(FASTA[s]) + "), skipped")
            continue
        L = fasta_lengths(FASTA[s])
        e = egg[egg["species"] == s][["gene_id", "protein_id"]].drop_duplicates()
        e["len"] = e["protein_id"].map(L)
        log("  " + s + ": " + str(round(100 * e["len"].notna().mean(), 1)) + "% of eggNOG protein IDs found in FASTA")
        glen.append(e.dropna().groupby("gene_id")["len"].mean().rename("prot_len").reset_index().assign(species=s))
    if glen:
        glen = pd.concat(glen)
        hvl = hv.merge(glen, left_on=["species", "gene"], right_on=["species", "gene_id"], how="left")["prot_len"].values
        rows = []
        for m in avail:
            for s in SPECIES:
                msk = (spc == SPECIES.index(s)) & ~np.isnan(hvl)
                if msk.sum() > 10:
                    nm = hv_norm[m][msk]
                    # constant-length vectors (e.g. rescaled) -> correlation undefined
                    rho = spearmanr(hvl[msk], nm).correlation if nm.std() / nm.mean() > 1e-4 else np.nan
                    rows.append({"model": m, "species": s, "n_genes_with_length": int(msk.sum()),
                                 "spearman_rho_norm_vs_length": rho,
                                 "median_protein_length": float(np.median(hvl[msk]))})
        d2 = pd.DataFrame(rows)
        save(d2, "D2_norm_vs_protein_length.csv")
        log(d2.round(3).to_string(index=False))
        log("rho near 1 => length grows with protein length (consistent with sum-pooling or a length-scaled layer)")

# ════════════════════════════════════════════════════════════
# D3. Global geometry
# ════════════════════════════════════════════════════════════
log("\n── D3. Global geometry (HV_shared) ──")


def unit(X):
    return X / np.clip(np.linalg.norm(X, axis=1, keepdims=True), 1e-12, None)


def species_r2(X):
    mu = X.mean(0)
    tot = ((X - mu) ** 2).sum()
    btw = sum((spc == i).sum() * ((X[spc == i].mean(0) - mu) ** 2).sum() for i in range(len(SPECIES)))
    return btw / tot


def eff_dim(X, n_sub=8000, seed=0):
    rng = np.random.default_rng(seed)
    Xs = X[rng.choice(len(X), min(len(X), n_sub), replace=False)]
    Xs = Xs - Xs.mean(0)
    _, sv, Vt = np.linalg.svd(Xs, full_matrices=False)
    ev = sv ** 2
    pr = ev.sum() ** 2 / (ev ** 2).sum()
    n90 = int(np.searchsorted(np.cumsum(ev) / ev.sum(), 0.9) + 1)
    return pr, n90, ev[0] / ev.sum(), Vt[0]


rows = []
for m in avail:
    X = EMB_HV[m]
    U = unit(X)
    s_ = U.sum(0)
    aniso = (s_ @ s_ - N) / (N * (N - 1))
    pr, n90, pc1, v1 = eff_dim(X)
    pc1_scores = (X - X.mean(0)) @ v1
    rows.append({"model": m, "dim": X.shape[1],
                 "mean_pairwise_cosine_all": aniso,
                 "participation_ratio": pr, "n_PCs_90pct_var": n90, "PC1_var_frac": pc1,
                 "R2_species_euclidean": species_r2(X), "R2_species_unit": species_r2(U),
                 "abs_corr_PC1_vs_norm": abs(np.corrcoef(pc1_scores, hv_norm[m])[0, 1])})
d3 = pd.DataFrame(rows)
save(d3, "D3_global_geometry.csv")
log(d3.round(3).to_string(index=False))
log("R2_species_euclidean = what SATURN's Euclidean KMeans sees; R2_species_unit = direction only")

# ════════════════════════════════════════════════════════════
# D4 / D5. Neighbours
# ════════════════════════════════════════════════════════════
og_sets = None
if egg is not None:
    e = egg.dropna(subset=["metazoa_OG"]).drop_duplicates(["species", "gene_id", "metazoa_OG"])
    hv_idx = pd.Series(np.arange(N), index=pd.MultiIndex.from_arrays([hv["species"], hv["gene"]]))
    e = e[pd.MultiIndex.from_arrays([e["species"], e["gene_id"]]).isin(hv_idx.index)]
    e["i"] = hv_idx.loc[pd.MultiIndex.from_arrays([e["species"], e["gene_id"]])].values
    og_of = e.groupby("i")["metazoa_OG"].apply(set)
    og_sets = [og_of.get(i, set()) for i in range(N)]
    # per target species: OG -> set of HV gene indices in that species
    OG_MEMBERS = []
    for t in range(len(SPECIES)):
        mm = {}
        for j in np.where(spc == t)[0]:
            for og in og_sets[j]:
                mm.setdefault(og, set()).add(j)
        OG_MEMBERS.append(mm)
    log("\nHV genes with a Metazoa OG: " + str(int(sum(1 for x in og_sets if x))) + " / " + str(N))


def neighbours(X, metric, k, chunk=1024):
    """Return all-species kNN and, per target species, the k nearest within that species."""
    if metric == "cosine":
        X = unit(X)
    sq = (X ** 2).sum(1)
    idx_t = [np.where(spc == t)[0] for t in range(len(SPECIES))]
    nn_all = np.empty((N, k), np.int64)
    nn_t = {t: np.empty((N, k), np.int64) for t in range(len(SPECIES))}
    for st in range(0, N, chunk):
        en = min(st + chunk, N)
        G = X[st:en] @ X.T
        D = -G if metric == "cosine" else sq[st:en, None] + sq[None, :] - 2 * G
        D[np.arange(en - st), np.arange(st, en)] = np.inf
        part = np.argpartition(D, k, axis=1)[:, :k]
        nn_all[st:en] = part
        for t, it in enumerate(idx_t):
            Dt = D[:, it]
            p = np.argpartition(Dt, k, axis=1)[:, :k]
            # order the k by distance so column 0 is the nearest
            o = np.argsort(np.take_along_axis(Dt, p, 1), axis=1)
            nn_t[t][st:en] = it[np.take_along_axis(p, o, 1)]
    return nn_all, nn_t


log("\n── D4. Cross-species neighbour mixing (k=" + str(K) + ") ──")
mix_rows, orth_rows = [], []
for m in avail:
    for metric in ["euclidean", "cosine"]:
        t1 = time.time()
        nn_all, nn_t = neighbours(EMB_HV[m], metric, K)
        xfrac = (spc[nn_all] != spc[:, None]).mean(1)
        for i, s in enumerate(SPECIES):
            msk = spc == i
            comp = np.bincount(spc[nn_all[msk]].ravel(), minlength=len(SPECIES)) / (msk.sum() * K)
            mix_rows.append({"model": m, "metric": metric, "species": s,
                             "xsp_nn_frac": xfrac[msk].mean(),
                             "random_expectation": (N - msk.sum()) / (N - 1),
                             **{"nn_from_" + t: comp[j] for j, t in enumerate(SPECIES)}})
        if og_sets is not None:
            for si, s in enumerate(SPECIES):
                for ti, t in enumerate(SPECIES):
                    if si == ti:
                        continue
                    tgt = np.where(spc == ti)[0]
                    members_t = OG_MEMBERS[ti]
                    # genes in s with >=1 OG partner among HV genes of t, and how many partners
                    q, n_part = [], []
                    for i in np.where(spc == si)[0]:
                        part = set()
                        for og in og_sets[i]:
                            part |= members_t.get(og, set())
                        if part:
                            q.append(i); n_part.append(len(part))
                    if not q:
                        continue
                    q, n_part = np.array(q), np.array(n_part)
                    hit1 = np.array([bool(og_sets[i] & og_sets[nn_t[ti][i, 0]]) for i in q])
                    hitk = np.array([any(og_sets[i] & og_sets[j] for j in nn_t[ti][i]) for i in q])
                    orth_rows.append({"model": m, "metric": metric, "query": s, "target": t,
                                      "n_query_genes": len(q), "hit_at_1": hit1.mean(),
                                      "hit_at_" + str(K): hitk.mean(),
                                      "random_hit_at_1": (n_part / len(tgt)).mean()})
        print("  " + m + " " + metric + " " + str(round(time.time() - t1)) + " s", flush=True)

d4 = pd.DataFrame(mix_rows)
save(d4, "D4_cross_species_mixing.csv")
log(d4[["model", "metric", "species", "xsp_nn_frac", "random_expectation"]].round(3).to_string(index=False))

d5 = None
if orth_rows:
    log("\n── D5. Ortholog recovery: nearest neighbour in target species shares a Metazoa OG ──")
    d5 = pd.DataFrame(orth_rows)
    save(d5, "D5_ortholog_nn_recovery.csv")
    piv = d5.pivot_table(index=["query", "target"], columns=["metric", "model"], values="hit_at_1")
    log(piv.round(3).to_string())
    summ = d5.groupby(["model", "metric"])[["hit_at_1", "hit_at_" + str(K), "random_hit_at_1"]].mean()
    log("\nMean over species pairs:\n" + summ.round(3).to_string())
    save(summ.reset_index(), "D5_ortholog_nn_recovery_summary.csv")

# ════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.titleweight": "bold", "axes.titlesize": 11})
FIG = OUT / "figures"

def pad_flat_ylim(ax, v):
    """Near-constant lengths: show +-10% instead of zooming into float noise."""
    v = np.asarray(v); mu = float(np.mean(v))
    if mu > 0 and (v.max() - v.min()) / mu < 0.01:
        ax.set_ylim(0.9 * mu, 1.1 * mu)


# FigD1: norm distributions, one panel per model (own y-axis: scales differ ~200x)
fig, axes = plt.subplots(1, len(avail), figsize=(3.6 * len(avail), 3.8))
for ax, m in zip(np.atleast_1d(axes), avail):
    data = [hv_norm[m][spc == i] for i in range(len(SPECIES))]
    vp = ax.violinplot(data, showextrema=False, showmedians=True)
    for b, s in zip(vp["bodies"], SPECIES):
        b.set_facecolor(SCOL[s]); b.set_edgecolor(SCOL[s]); b.set_alpha(0.6)
    vp["cmedians"].set_color(INK)
    ax.set_xticks(range(1, len(SPECIES) + 1), SPECIES)
    ax.set_title(m, fontsize=10); ax.set_ylabel("Vector length (L2 norm)")
    ax.ticklabel_format(axis="y", useOffset=False, style="plain")
    pad_flat_ylim(ax, hv_norm[m])
fig.suptitle("Gene embedding vector length, HV genes", fontweight="bold")
fig.tight_layout(); fig.savefig(FIG / "FigD1_vector_length.png", dpi=200); plt.close(fig)

# FigD2: norm vs protein length
if d2 is not None:
    fig, axes = plt.subplots(1, len(avail), figsize=(3.8 * len(avail), 3.6))
    for ax, m in zip(np.atleast_1d(axes), avail):
        for i, s in enumerate(SPECIES):
            msk = (spc == i) & ~np.isnan(hvl)
            ax.scatter(hvl[msk], hv_norm[m][msk], s=2, alpha=0.25, color=SCOL[s], lw=0, label=s)
        ax.set_xscale("log"); ax.set_xlabel("Protein length (aa, log)"); ax.set_ylabel("Vector length")
        ax.ticklabel_format(axis="y", useOffset=False, style="plain")
        pad_flat_ylim(ax, hv_norm[m])
        r = d2[d2["model"] == m]
        ax.set_title(m + "\nSpearman rho " + ", ".join(s + " " + ("n/a" if np.isnan(v) else str(round(v, 2))) for s, v in zip(r["species"], r["spearman_rho_norm_vs_length"])), fontsize=8)
    handles = [plt.Line2D([], [], marker="o", ls="", color=SCOL[s], label=s) for s in SPECIES]
    np.atleast_1d(axes)[0].legend(handles=handles, frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "FigD2_length_vs_protein_length.png", dpi=200); plt.close(fig)

# FigD3: PCA as SATURN's KMeans sees it (Euclidean, uncentred scale kept)
fig, axes = plt.subplots(1, len(avail), figsize=(3.8 * len(avail), 3.8))
rng = np.random.default_rng(42)
order = rng.permutation(N)
for ax, m in zip(np.atleast_1d(axes), avail):
    X = EMB_HV[m] - EMB_HV[m].mean(0)
    _, S_, Vt = np.linalg.svd(X[rng.choice(N, min(N, 8000), replace=False)], full_matrices=False)
    P = X @ Vt[:2].T
    ve = S_[:2] ** 2 / (S_ ** 2).sum()
    ax.scatter(P[order, 0], P[order, 1], s=2, c=[SCOL[s] for s in spv[order]], alpha=0.5, lw=0)
    ax.set_xlabel("PC1 (" + str(round(100 * ve[0], 1)) + "%)"); ax.set_ylabel("PC2 (" + str(round(100 * ve[1], 1)) + "%)")
    ax.set_title(m, fontsize=10)
handles = [plt.Line2D([], [], marker="o", ls="", color=SCOL[s], label=s) for s in SPECIES]
np.atleast_1d(axes)[-1].legend(handles=handles, frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "FigD3_PCA_by_species.png", dpi=200); plt.close(fig)

# FigD4: cross-species neighbour mixing
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
x = np.arange(len(SPECIES))
w = 0.8 / len(avail)
for ax, metric in zip(axes, ["euclidean", "cosine"]):
    for j, m in enumerate(avail):
        d = d4[(d4["model"] == m) & (d4["metric"] == metric)].set_index("species").loc[SPECIES]
        ax.bar(x - 0.4 + w / 2 + j * w, d["xsp_nn_frac"], w * 0.95, color=MCOL[m], label=m)
    rnd = d4[(d4["model"] == avail[0]) & (d4["metric"] == metric)].set_index("species").loc[SPECIES, "random_expectation"]
    ax.scatter(x, rnd, marker="_", s=600, color=INK, lw=2, label="random mixing", zorder=3)
    ax.set_xticks(x, SPECIES); ax.set_title("Cross-species neighbours (" + metric + ")")
axes[0].set_ylabel("Fraction of " + str(K) + " nearest neighbours\nfrom other species")
axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(FIG / "FigD4_cross_species_mixing.png", dpi=200); plt.close(fig)

# FigD5: ortholog hit@1 per species pair (pairs averaged over both directions)
if d5 is not None:
    d5p = d5.assign(pair=[" - ".join(sorted([a, b], key=SPECIES.index)) for a, b in zip(d5["query"], d5["target"])])
    agg = d5p.groupby(["pair", "model", "metric"], sort=False)[["hit_at_1", "random_hit_at_1"]].mean().reset_index()
    pairs = list(dict.fromkeys(agg["pair"]))
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    x = np.arange(len(pairs))
    for ax, metric in zip(axes, ["euclidean", "cosine"]):
        for j, m in enumerate(avail):
            d = agg[(agg["model"] == m) & (agg["metric"] == metric)].set_index("pair").reindex(pairs)
            ax.bar(x - 0.4 + w / 2 + j * w, d["hit_at_1"], w * 0.95, color=MCOL[m], label=m)
        r = agg[(agg["model"] == avail[0]) & (agg["metric"] == metric)].set_index("pair").reindex(pairs)["random_hit_at_1"]
        ax.scatter(x, r, marker="_", s=500, color=INK, lw=2, label="random", zorder=3)
        ax.set_xticks(x, pairs, rotation=30, ha="right"); ax.set_title("Nearest neighbour is an OG partner (" + metric + ")")
    axes[0].set_ylabel("Hit@1"); axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "FigD5_ortholog_nn_recovery.png", dpi=200); plt.close(fig)

log("\nDone in " + str(round(time.time() - t0)) + " s. Outputs: " + str(OUT))
(OUT / "summary.txt").write_text("\n".join(SUMMARY) + "\n")
