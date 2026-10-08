#!/usr/bin/env python3
"""
05i_embedding_standardization_check.py

Can simple post-processing rescue the ESMC-600M embeddings? Cheap test, no SATURN run:
apply each transformation to the gene embeddings, then measure how often a gene's nearest
neighbour in another species shares its eggNOG Metazoa OG (same metric as 05f D5).

Transformations (applied to every gene of all 4 species; statistics from ALL genes in the files):
  raw         as delivered
  unit        each vector scaled to length 1 (= the 05e rescaling, up to a constant)
  center      subtract the mean vector
  zscore      per dimension: subtract mean, divide by SD (fix for "rogue"/outlier dimensions)
  zscore_unit zscore, then length 1
  abtt1/2/3   "all-but-the-top": center, then remove the top 1/2/3 principal components
Run on ESM-1b too, as a control (a fix that helps ESMC should not need to help ESM-1b).

Metrics per model x transform x metric (euclidean = what SATURN's KMeans sees; cosine):
  hit@1, hit@10 (mean over the 12 directed species pairs; also Mlei pairs only),
  cross-species neighbour fraction, effective dimensionality (participation ratio).

--write_pt TRANSFORM writes ESMC embeddings transformed that way (all genes) to
03_embeddings_ESMC600M_<TRANSFORM>/ so a SATURN macrogene run can use them if it wins.
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

BASE = Path("/scratch/dark_genes/SATURN_Mnemi")
RUNS = BASE / "04_saturn_runs"
ap = argparse.ArgumentParser()
ap.add_argument("--hv_run_a", default=str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESM1b_macrogene"))
ap.add_argument("--hv_run_b", default=str(RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESMC600M_macrogene"))
ap.add_argument("--eggnog", default=str(BASE / "01_proteomes/eggnog_merged/ALL4_gene_protein_metazoaOG_slim.csv"))
ap.add_argument("--species", default="Mlei,Cgig,Crob,Drer")
ap.add_argument("--knn_k", type=int, default=10)
ap.add_argument("--write_pt", default="", help="ESMC transform to write as SATURN .pt files (e.g. zscore)")
ap.add_argument("--out_dir", default=str(RUNS / "20261008_embedding_standardization_check"))
args = ap.parse_args()

SPECIES = args.species.split(",")
FOLDER = {"Mlei": "Mnemi"}
MODELS = {
    "ESM1b": lambda sp: BASE / "03_embeddings" / sp / f"{sp}_gene_embeddings.pt",
    "ESMC600M": lambda sp: BASE / "03_embeddings_ESMC600M" / FOLDER.get(sp, sp) / f"{FOLDER.get(sp, sp)}_gene_embeddings.pt",
}
TRANSFORMS = ["raw", "unit", "center", "zscore", "zscore_unit", "abtt1", "abtt2", "abtt3"]
K = args.knn_k
OUT = Path(args.out_dir)
(OUT / "figures").mkdir(parents=True, exist_ok=True)
LOG = []


def log(m=""):
    print(m, flush=True); LOG.append(str(m))


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


def resolve(keys, genes):
    lower = {}
    for k in keys:
        lower.setdefault(k.lower(), []).append(k)
    return {g: (g if g in keys else lower[g.lower()][0]) for g in genes
            if g in keys or len(lower.get(g.lower(), [])) == 1}


t0 = time.time()
shared = sorted(pkl_keys(args.hv_run_a) & pkl_keys(args.hv_run_b))
hv = pd.DataFrame({"species": [k.split("_", 1)[0] for k in shared], "gene": [k.split("_", 1)[1] for k in shared]})
hv = pd.concat([hv[hv.species == s] for s in SPECIES]).reset_index(drop=True)

# ── load all genes per model; remember HV rows ──
ALL, HVIDX, KEYS = {}, {}, {}
ok = np.ones(len(hv), bool)
for m, pf in MODELS.items():
    blocks, keys, hv_rows, off = [], [], np.full(len(hv), -1), 0
    for s in SPECIES:
        d = load_pt(pf(s))
        names = list(d.keys())
        blocks.append(np.vstack([vec(d[k]) for k in names]))
        keys += [(s, k) for k in names]
        r = resolve(set(names), hv.loc[hv.species == s, "gene"])
        pos = {k: i for i, k in enumerate(names)}
        for i in np.where(hv.species.values == s)[0]:
            gname = hv.at[i, "gene"]
            if gname in r:
                hv_rows[i] = off + pos[r[gname]]
        off += len(names)
        del d
    ALL[m], KEYS[m] = np.vstack(blocks), keys
    HVIDX[m] = hv_rows
    ok &= hv_rows >= 0
    log(f"{m}: {ALL[m].shape[0]} genes in files, dim {ALL[m].shape[1]}")
hv = hv[ok].reset_index(drop=True)
for m in MODELS:
    HVIDX[m] = HVIDX[m][ok]
spc = pd.Categorical(hv.species, categories=SPECIES).codes
N = len(hv)
log(f"HV genes analysed: {N}")


def transform(X, name, hv_rows):
    if name == "raw":
        return X
    if name == "unit":
        return X / np.clip(np.linalg.norm(X, axis=1, keepdims=True), 1e-12, None)
    mu = X.mean(0, keepdims=True)
    if name == "center":
        return X - mu
    if name.startswith("zscore"):
        Z = (X - mu) / np.clip(X.std(0, keepdims=True), 1e-8, None)
        return Z / np.clip(np.linalg.norm(Z, axis=1, keepdims=True), 1e-12, None) if name == "zscore_unit" else Z
    if name.startswith("abtt"):
        k = int(name[4:])
        C = X - mu
        rng = np.random.default_rng(0)
        sub = C[rng.choice(len(C), min(len(C), 20000), replace=False)]
        _, _, Vt = np.linalg.svd(sub, full_matrices=False)
        U = Vt[:k]
        return C - (C @ U.T) @ U
    raise ValueError(name)


def part_ratio(X):
    rng = np.random.default_rng(1)
    Xs = X[rng.choice(len(X), min(len(X), 8000), replace=False)]
    ev = np.linalg.svd(Xs - Xs.mean(0), compute_uv=False) ** 2
    return ev.sum() ** 2 / (ev ** 2).sum()


# ── eggNOG OG sets ──
egg = pd.read_csv(args.eggnog).dropna(subset=["metazoa_OG"])
egg = egg[egg.species.isin(SPECIES)].drop_duplicates(["species", "gene_id", "metazoa_OG"])
og_of = egg.groupby(["species", "gene_id"])["metazoa_OG"].apply(set)
og_sets = [og_of.get((s, gname), set()) for s, gname in zip(hv.species, hv.gene)]
OGM = []
for t in range(len(SPECIES)):
    mm = {}
    for j in np.where(spc == t)[0]:
        for og in og_sets[j]:
            mm.setdefault(og, set()).add(j)
    OGM.append(mm)
QUERY = {}
for si in range(len(SPECIES)):
    for ti in range(len(SPECIES)):
        if si == ti:
            continue
        q = [i for i in np.where(spc == si)[0] if any(og in OGM[ti] for og in og_sets[i])]
        QUERY[(si, ti)] = np.array(q, dtype=int)


def evaluate(H, metric):
    if metric == "cosine":
        H = H / np.clip(np.linalg.norm(H, axis=1, keepdims=True), 1e-12, None)
    sq = (H ** 2).sum(1)
    idx_t = [np.where(spc == t)[0] for t in range(len(SPECIES))]
    nn_t = {t: np.empty((N, K), np.int64) for t in range(len(SPECIES))}
    xsp = np.empty(N)
    for st in range(0, N, 1024):
        en = min(st + 1024, N)
        G = H[st:en] @ H.T
        D = -G if metric == "cosine" else sq[st:en, None] + sq[None, :] - 2 * G
        D[np.arange(en - st), np.arange(st, en)] = np.inf
        allk = np.argpartition(D, K, axis=1)[:, :K]
        xsp[st:en] = (spc[allk] != spc[st:en, None]).mean(1)
        for t, it in enumerate(idx_t):
            Dt = D[:, it]
            p = np.argpartition(Dt, K, axis=1)[:, :K]
            o = np.argsort(np.take_along_axis(Dt, p, 1), axis=1)
            nn_t[t][st:en] = it[np.take_along_axis(p, o, 1)]
    res = {}
    for (si, ti), q in QUERY.items():
        h1 = np.mean([bool(og_sets[i] & og_sets[nn_t[ti][i, 0]]) for i in q])
        hk = np.mean([any(og_sets[i] & og_sets[j] for j in nn_t[ti][i]) for i in q])
        res[(si, ti)] = (h1, hk)
    return res, xsp


rows, pair_rows = [], []
for m in MODELS:
    for tf in TRANSFORMS:
        X = transform(ALL[m], tf, HVIDX[m])
        H = X[HVIDX[m]]
        pr = part_ratio(H)
        for metric in ["euclidean", "cosine"]:
            t1 = time.time()
            res, xsp = evaluate(H, metric)
            h1 = np.mean([v[0] for v in res.values()]); hk = np.mean([v[1] for v in res.values()])
            mlei = [v[0] for (a, b), v in res.items() if SPECIES[a] == "Mlei" or SPECIES[b] == "Mlei"]
            rows.append({"model": m, "transform": tf, "metric": metric, "hit_at_1": h1, f"hit_at_{K}": hk,
                         "hit_at_1_Mlei_pairs": np.mean(mlei), "xsp_nn_frac": xsp.mean(),
                         "participation_ratio": pr})
            for (a, b), v in res.items():
                pair_rows.append({"model": m, "transform": tf, "metric": metric, "query": SPECIES[a],
                                  "target": SPECIES[b], "hit_at_1": v[0], f"hit_at_{K}": v[1]})
            print(f"  {m} {tf} {metric}: hit@1 {h1:.3f} ({time.time() - t1:.0f}s)", flush=True)
res = pd.DataFrame(rows)
res.to_csv(OUT / "I1_standardization_summary.csv", index=False)
pd.DataFrame(pair_rows).to_csv(OUT / "I2_standardization_by_species_pair.csv", index=False)
log("\nNearest-neighbour ortholog recovery (hit@1, mean of 12 species pairs):")
log(res.pivot_table(index="transform", columns=["model", "metric"], values="hit_at_1").loc[TRANSFORMS].round(3).to_string())
log("\nSame, Mlei pairs only:")
log(res.pivot_table(index="transform", columns=["model", "metric"], values="hit_at_1_Mlei_pairs").loc[TRANSFORMS].round(3).to_string())
log("\nCross-species neighbour fraction (euclidean; random ~0.75):")
log(res[res.metric == "euclidean"].pivot(index="transform", columns="model", values="xsp_nn_frac").loc[TRANSFORMS].round(3).to_string())
log("\nEffective dimensionality (participation ratio):")
log(res[res.metric == "euclidean"].pivot(index="transform", columns="model", values="participation_ratio").loc[TRANSFORMS].round(1).to_string())

# ── figure ──
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#8a8a85", "axes.titleweight": "bold", "axes.titlesize": 10})
MC = {"ESM1b": "#eb6834", "ESMC600M": "#1baf7a"}
fig, axes = plt.subplots(1, 2, figsize=(12, 3.9), sharey=True)
x = np.arange(len(TRANSFORMS))
for ax, metric in zip(axes, ["euclidean", "cosine"]):
    for k, m in enumerate(MODELS):
        v = res[(res.model == m) & (res.metric == metric)].set_index("transform").loc[TRANSFORMS, "hit_at_1"]
        ax.bar(x - 0.2 + k * 0.4, v, 0.38, color=MC[m], label=m)
    ax.set_xticks(x, TRANSFORMS, rotation=30, ha="right")
    ax.set_title(f"Nearest neighbour in other species shares an OG ({metric})")
axes[0].set_ylabel("Hit@1 (mean of species pairs)"); axes[0].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(OUT / "figures" / "FigI1_standardization_hit1.png", dpi=200); plt.close(fig)

# ── optional: write transformed ESMC for SATURN ──
if args.write_pt:
    import torch
    X = transform(ALL["ESMC600M"], args.write_pt, HVIDX["ESMC600M"])
    keys = KEYS["ESMC600M"]
    for s in SPECIES:
        f = FOLDER.get(s, s)
        dst = BASE / f"03_embeddings_ESMC600M_{args.write_pt}" / f / f"{f}_gene_embeddings.pt"
        dst.parent.mkdir(parents=True, exist_ok=True)
        d = {k: torch.tensor(X[i]) for i, (sp_, k) in enumerate(keys) if sp_ == s}
        torch.save(d, dst)
        log(f"wrote {len(d)} {s} genes -> {dst}")

log(f"\nDone in {time.time() - t0:.0f} s. Outputs: {OUT}")
(OUT / "summary.txt").write_text("\n".join(LOG) + "\n")
