#!/usr/bin/env python3
"""
05t_esm1b_longest_isoform.py   (2026-10-27)

Isoform test for the 6-species ESM-1b benchmark (RESULTS_6sp_benchmark_ESM1b_vs_ESMCpostnorm_20261027):
human pairs score lowest with every partner. Hypothesis: averaging all Ensembl pep.all isoforms per gene
(18.6% of human pep.all proteins are < 100 aa) blurs the human gene vectors.

Mode 1 (default): rebuild ESM-1b gene vectors from ONE protein per gene (the longest in {SP}_filtered.faa),
using the per-protein .pt files that 02_embed.sh / 02b_resume already wrote. No new embedding.
  in : 03_embeddings/{SP}/{SP}_gene_embeddings.pt              (current: all isoforms averaged)
       03_embeddings/{SP}/{SP}_gene_to_protein.filtered.json    (gene -> proteins with a .pt)
       03_embeddings/{SP}/{SP}_pt_embeddings/{protein}.pt       (fair-esm extract.py, mean_representations[33])
       03_embeddings/{SP}/{SP}_filtered.faa                     (protein lengths)
  out: 03_embeddings_ESM1b_longest/{F}/{F}_gene_embeddings.pt   (F = Mnemi for Mlei, the layout 05i expects)
       03_embeddings_ESM1b_longest/longest_isoform_summary.tsv, longest_isoform_choice_{SP}.tsv
  Same gene set as the averaged file, so SATURN HV genes and the 05i gene universe are unchanged.
  Single-protein genes are copied unchanged. Ties in length -> alphabetically first protein ID.
  Safety check: for --n_check random multi-isoform genes, the mean of all isoform vectors is recomputed from
  the .pt files and must match the existing averaged vector (confirms layer / file format). Exits 1 if not.

Mode 2: --report_pairs <I2_standardization_by_species_pair.csv>  prints per-pair hit@1 (raw, cosine) for each
  model, and the human rows side by side.

Run on Pegasus via 05u_isoform_longest_test.sh (normal queue).
"""
import argparse
import json
import os
import random
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

BASE = Path(os.environ.get("SATURN_BASE", "/scratch/dark_genes/SATURN_Mnemi"))
FOLDER = {"Mlei": "Mnemi"}

ap = argparse.ArgumentParser()
ap.add_argument("--species", default="Mlei,Cgig,Crob,Drer,Dmel,Hsap")
ap.add_argument("--out_dir", default=str(BASE / "03_embeddings_ESM1b_longest"))
ap.add_argument("--layer", type=int, default=33)
ap.add_argument("--n_check", type=int, default=200)
ap.add_argument("--report_pairs", default="", help="I2_standardization_by_species_pair.csv from 05i (mode 2)")
args = ap.parse_args()


def report_pairs(path):
    d = pd.read_csv(path)
    d = d[(d["transform"] == "raw") & (d["metric"] == "cosine")].copy()
    for m in d.model.unique():
        p = d[d.model == m].pivot(index="query", columns="target", values="hit_at_1").round(3)
        print(f"\n== hit@1, raw, cosine: {m} (rows = query species, columns = target species) ==")
        print(p.to_string())
    h = d[np.logical_or(d["query"] == "Hsap", d["target"] == "Hsap")]
    w = h.pivot_table(index=["query", "target"], columns="model", values="hit_at_1").round(3)
    print("\n== pairs involving Hsap, raw, cosine ==")
    print(w.to_string())
    nh = d[(d["query"] != "Hsap") & (d["target"] != "Hsap")].groupby("model").hit_at_1.mean().round(3)
    print("\nmean hit@1 over pairs WITHOUT Hsap:")
    print(nh.to_string())
    print("mean hit@1 over pairs WITH Hsap:")
    print(h.groupby("model").hit_at_1.mean().round(3).to_string())


if args.report_pairs:
    report_pairs(args.report_pairs)
    sys.exit(0)


def read_lengths(faa):
    lens, cur, n = {}, None, 0
    with open(faa) as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    lens[cur] = n
                cur, n = line[1:].split()[0].strip(), 0
            else:
                n += len(line.strip().rstrip("*"))
    if cur is not None:
        lens[cur] = n
    return lens


def load_protein(pt_dir, pid):
    d = torch.load(pt_dir / f"{pid}.pt", map_location="cpu")
    return d["mean_representations"][args.layer].float()


def cos(a, b):
    a = a.float().ravel(); b = b.float().ravel()
    return float(torch.dot(a, b) / (a.norm() * b.norm()).clamp_min(1e-12))


OUT = Path(args.out_dir)
OUT.mkdir(parents=True, exist_ok=True)
random.seed(0)
summary = []
failed = False
t0 = time.time()

for sp in [s.strip() for s in args.species.split(",") if s.strip()]:
    t1 = time.time()
    ed = BASE / "03_embeddings" / sp
    orig_p = ed / f"{sp}_gene_embeddings.pt"
    map_p = ed / f"{sp}_gene_to_protein.filtered.json"
    pt_dir = ed / f"{sp}_pt_embeddings"
    faa = ed / f"{sp}_filtered.faa"
    print(f"\n=== {sp} ===", flush=True)

    orig = torch.load(orig_p, map_location="cpu")
    gmap = json.load(open(map_p))
    genes = list(orig.keys())
    no_map = [g for g in genes if g not in gmap]
    print(f"  genes in averaged file: {len(genes):,}; in filtered map: {len(gmap):,}; "
          f"averaged-file genes missing from map: {len(no_map)}")
    if no_map:
        print("  ERROR: genes without a protein list; cannot rebuild them", no_map[:10])
        failed = True
        continue

    multi = [g for g in genes if len(gmap[g]) > 1]
    if multi and not pt_dir.is_dir():
        print(f"  ERROR: {len(multi):,} multi-isoform genes but no per-protein folder {pt_dir}")
        failed = True
        continue
    lens = read_lengths(faa)
    prots_used = sorted({p for g in genes for p in gmap[g]})
    no_len = [p for p in prots_used if p not in lens]
    if no_len:
        print(f"  WARNING: {len(no_len)} proteins missing from {faa.name} (length -1): {no_len[:5]}")
    L = np.array([lens.get(p, -1) for p in prots_used])

    # safety check: recomputed isoform mean == existing averaged vector
    chk = random.sample(multi, min(args.n_check, len(multi)))
    maxdiff = 0.0
    for g in chk:
        m = torch.stack([load_protein(pt_dir, p) for p in gmap[g]]).mean(0)
        o = orig[g].float().ravel()
        maxdiff = max(maxdiff, float((m.ravel() - o).abs().max() / o.abs().max()))
    print(f"  check: {len(chk)} multi-isoform genes, max relative diff recomputed-mean vs file = {maxdiff:.2e}")
    if chk and maxdiff > 1e-3:
        print("  ERROR: recomputed mean does not match the averaged file (layer/format?). Not writing.")
        failed = True
        continue

    new, rows, cs = {}, [], []
    n_trunc = 0
    for g in genes:
        ps = gmap[g]
        if len(ps) == 1:
            new[g] = orig[g].clone()
            continue
        best = sorted(ps, key=lambda p: (-lens.get(p, -1), p))[0]
        v = load_protein(pt_dir, best).to(orig[g].dtype).reshape(orig[g].shape)
        new[g] = v
        c = cos(v, orig[g])
        cs.append(c)
        bl = lens.get(best, -1)
        n_trunc += bl > 1022
        short = sum(lens.get(p, -1) < 100 for p in ps)
        rows.append({"gene": g, "n_isoforms": len(ps), "chosen": best, "chosen_len": bl,
                     "min_len": min(lens.get(p, -1) for p in ps), "n_lt100": short, "cos_longest_vs_avg": round(c, 4)})

    f = FOLDER.get(sp, sp)
    (OUT / f).mkdir(parents=True, exist_ok=True)
    torch.save(new, OUT / f / f"{f}_gene_embeddings.pt")
    pd.DataFrame(rows).to_csv(OUT / f"longest_isoform_choice_{sp}.tsv", sep="\t", index=False)

    cs = np.array(cs) if cs else np.array([np.nan])
    s = {"species": sp, "genes": len(genes), "multi_isoform_genes": len(multi),
         "pct_multi": round(100 * len(multi) / len(genes), 1),
         "mean_isoforms_multi": round(float(np.mean([len(gmap[g]) for g in multi])), 2) if multi else 0,
         "proteins": len(prots_used), "pct_proteins_lt100aa": round(100 * float(np.mean(L < 100)), 1),
         "cos_longest_vs_avg_median": round(float(np.nanmedian(cs)), 4),
         "cos_p10": round(float(np.nanpercentile(cs, 10)), 4),
         "pct_multi_cos_lt_0.95": round(100 * float(np.nanmean(cs < 0.95)), 1),
         "longest_over_1022aa": int(n_trunc), "check_maxdiff": f"{maxdiff:.1e}",
         "seconds": round(time.time() - t1)}
    summary.append(s)
    print("  " + "  ".join(f"{k}={v}" for k, v in s.items() if k != "species"), flush=True)
    del orig, new

if summary:
    sm = pd.DataFrame(summary)
    sm.to_csv(OUT / "longest_isoform_summary.tsv", sep="\t", index=False)
    print("\n=== Summary (cos = longest-isoform vector vs current averaged vector, multi-isoform genes) ===")
    print(sm.drop(columns=["seconds"]).to_string(index=False))
print(f"\nDone in {time.time() - t0:.0f} s. Output: {OUT}")
sys.exit(1 if failed else 0)
