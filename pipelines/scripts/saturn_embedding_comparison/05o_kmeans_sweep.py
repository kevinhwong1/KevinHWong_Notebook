#!/usr/bin/env python3
"""
05o_kmeans_sweep.py  (runs on Pegasus)

How does the number of macrogenes (K) change the ortholog benchmark, for ESM-1b vs ESMC-600M?
SATURN seeds its macrogenes with KMeans on the stacked gene embeddings (train-saturn.py ->
model/saturn_model.py::make_centroids). This script repeats exactly that KMeans call,
    KMeans(n_clusters=K, random_state=seed).fit(X)          (no other arguments, as in SATURN)
for a range of K, without running SATURN. Only one thing differs between models: the embeddings.

Gene set: genes in BOTH macrogene-only runs (ESM1b_4sp and ESMC600M_resc_4sp genes_to_macrogenes.pkl keys),
in the ESM-1b run's key order, so both models cluster the identical genes.
Embeddings are used as SATURN used them (ESM-1b: 03_embeddings; ESMC: 03_embeddings_ESMC600M_rescaled,
i.e. the 05e length-normalized vectors of the 20261005 run). No other transformation.

Modes
  --mode sweep    --model {ESM1b,ESMC600M_rescaled,ESMC600M_raw} --k K
                  writes labels/<model>_K<K>_seed<seed>.csv.gz (key, label) + a one-line timing/metadata json
  --mode validate checks that KMeans is a fair stand-in for SATURN's final macrogenes:
                  (a) for each K3000 run: argmax of centroids_init.pkl vs argmax of the final
                      genes_to_macrogenes.pkl (how much does pretraining move genes?)
                  (b) reproduce centroids_init.pkl: refit KMeans(3000, random_state=42) on the run's own
                      genes in the pkl's own order and compare with the init argmax (ARI 1 = identical)
                  writes validate.csv
Scoring is done locally by 05p_score_kmeans_sweep.py.
"""
import argparse
import json
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd

import sys
BASE = Path(next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--base=")), "/scratch/dark_genes/SATURN_Mnemi"))
RUNS = BASE / "04_saturn_runs"
RUN_DIRS = {
    "ESM1b": RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESM1b_macrogene",
    "ESMC600M_rescaled": RUNS / "20261005_Mlei_Cgig_Crob_Drer_K3000_ESMC600M_rescaled_macrogene",
    "ESMC600M_raw": RUNS / "20261002_Mlei_Cgig_Crob_Drer_K3000_ESMC600M_macrogene",
}
FOLDER = {"Mlei": "Mnemi"}
EMB = {
    "ESM1b": lambda sp: BASE / "03_embeddings" / sp / f"{sp}_gene_embeddings.pt",
    "ESMC600M_rescaled": lambda sp: BASE / "03_embeddings_ESMC600M_rescaled" / FOLDER.get(sp, sp) / f"{FOLDER.get(sp, sp)}_gene_embeddings.pt",
    "ESMC600M_raw": lambda sp: BASE / "03_embeddings_ESMC600M" / FOLDER.get(sp, sp) / f"{FOLDER.get(sp, sp)}_gene_embeddings.pt",
}

ap = argparse.ArgumentParser()
ap.add_argument("--base", default=str(BASE), help="project root (use --base=PATH form; for testing)")
ap.add_argument("--mode", choices=["sweep", "validate"], required=True)
ap.add_argument("--model", choices=list(EMB), default="ESM1b")
ap.add_argument("--k", type=int, default=3000)
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--gene_set_runs", default="ESM1b,ESMC600M_rescaled",
                help="runs whose genes_to_macrogenes keys are intersected to define the gene set")
ap.add_argument("--out_dir", default=str(RUNS / "20261009_kmeans_sweep"))
args = ap.parse_args()

OUT = Path(args.out_dir)
(OUT / "labels").mkdir(parents=True, exist_ok=True)


def log(m=""):
    print(m, flush=True)


def vec(v):
    if hasattr(v, "detach"):
        v = v.detach().cpu().numpy()
    return np.asarray(v, dtype=np.float32).ravel()


def load_pkl(run, pattern):
    hits = sorted((RUN_DIRS[run] / "saturn_results").glob(pattern)) if pattern != "centroids_init.pkl" \
        else [RUN_DIRS[run] / "centroids_init.pkl"]
    if len(hits) != 1 or not hits[0].exists():
        raise RuntimeError(f"Expected one {pattern} for {run}: {hits}")
    return pickle.load(open(hits[0], "rb"))


def embed_matrix(model, keys):
    """stack embeddings for 'Species_gene' keys, in the given order (exact name, then case-insensitive)"""
    import torch
    sp_of = [k.split("_", 1)[0] for k in keys]
    rows = [None] * len(keys)
    for s in dict.fromkeys(sp_of):
        d = torch.load(EMB[model](s), map_location="cpu")
        lower = {}
        for n in d:
            lower.setdefault(n.lower(), []).append(n)
        for i, k in enumerate(keys):
            if sp_of[i] != s:
                continue
            gname = k.split("_", 1)[1]
            n = gname if gname in d else (lower.get(gname.lower(), [None])[0] if len(lower.get(gname.lower(), [])) == 1 else None)
            if n is None:
                raise KeyError(f"{model}: no embedding for {k}")
            rows[i] = vec(d[n])
        del d
    return np.vstack(rows)


def kmeans(X, k, seed):
    import sklearn
    from sklearn.cluster import KMeans
    t = time.time()
    km = KMeans(n_clusters=k, random_state=seed).fit(X)   # identical call to SATURN's make_centroids
    meta = {"sklearn": sklearn.__version__, "n_init_used": getattr(km, "_n_init", getattr(km, "n_init", None)),
            "n_iter": int(km.n_iter_), "inertia": float(km.inertia_), "seconds": round(time.time() - t, 1)}
    return km.labels_, meta


if args.mode == "sweep":
    runs = args.gene_set_runs.split(",")
    keysets = [list(load_pkl(r, "*genes_to_macrogenes*.pkl").keys()) for r in runs]
    common = set(keysets[0]).intersection(*map(set, keysets[1:]))
    keys = [k for k in keysets[0] if k in common]
    log(f"Gene set: {len(keys)} genes common to {runs}  " +
        str(pd.Series([k.split('_', 1)[0] for k in keys]).value_counts().to_dict()))
    X = embed_matrix(args.model, keys)
    log(f"{args.model}: X {X.shape}, mean L2 norm {np.linalg.norm(X, axis=1).mean():.2f}")
    labels, meta = kmeans(X, args.k, args.seed)
    meta.update({"model": args.model, "k": args.k, "seed": args.seed, "n_genes": len(keys),
                 "n_nonempty_clusters": int(len(np.unique(labels)))})
    log(json.dumps(meta))
    stem = f"{args.model}_K{args.k}_seed{args.seed}"
    pd.DataFrame({"key": keys, "label": labels}).to_csv(OUT / "labels" / f"{stem}.csv.gz", index=False)
    (OUT / "labels" / f"{stem}.json").write_text(json.dumps(meta) + "\n")
    log("Done: " + str(OUT / "labels" / f"{stem}.csv.gz"))

else:
    from sklearn.metrics import adjusted_rand_score as ari
    rows = []
    for run in ["ESM1b", "ESMC600M_rescaled"]:
        init = load_pkl(run, "centroids_init.pkl")
        final = load_pkl(run, "*genes_to_macrogenes*.pkl")
        keys = [k for k in init if k in final]
        a0 = np.array([int(np.argmax(np.asarray(init[k]))) for k in keys])
        a1 = np.array([int(np.argmax(np.asarray(final[k]))) for k in keys])
        row = {"run": run, "n_genes": len(keys), "K": len(np.asarray(init[keys[0]])),
               "pct_same_mg_init_vs_final": 100 * float((a0 == a1).mean()), "ARI_init_vs_final": ari(a0, a1)}
        log(f"{run}: init vs final argmax: {row['pct_same_mg_init_vs_final']:.1f}% same, ARI {row['ARI_init_vs_final']:.3f}")
        # (b) reproduce SATURN's KMeans on this run's own genes, in the pkl order (= SATURN's row order)
        X = embed_matrix(run, list(init.keys()))
        lab, meta = kmeans(X, row["K"], args.seed)
        pos = {k: i for i, k in enumerate(init.keys())}
        lab_k = lab[[pos[k] for k in keys]]
        row.update({"ARI_refit_vs_init": ari(lab_k, a0), "ARI_refit_vs_final": ari(lab_k, a1),
                    **{"refit_" + k: v for k, v in meta.items()}})
        log(f"{run}: refit KMeans vs SATURN init: ARI {row['ARI_refit_vs_init']:.3f}; vs final: {row['ARI_refit_vs_final']:.3f}  {meta}")
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "validate.csv", index=False)
    log("Done: " + str(OUT / "validate.csv"))
