#!/usr/bin/env python3
"""
05y_check_gene_sets.py   (2026-10-29)

Pre-flight check for the 2x2 full-SATURN comparison (model x protein-per-gene). SATURN keeps the genes that
are both in the h5ad and in the embedding file, then picks HV genes from those, so identical gene sets across
the four embedding sets => identical HV genes (same cells, same h5ad). This script confirms that, per species:

  ESM1b_avg      03_embeddings/{SP}/{SP}_gene_embeddings.pt                         (Mlei: 03_embeddings/Mlei/)
  ESMC_avg       03_embeddings_ESMC600M_postnorm_avgmatched/{F}/{F}_gene_embeddings.pt   (05q --match_esm1b)
  ESM1b_longest  03_embeddings_ESM1b_longest/{F}/{F}_gene_embeddings.pt             (05t)
  ESMC_longest   03_embeddings_ESMC600M_postnorm_longest/{F}/{F}_gene_embeddings.pt  (05q --longest + restriction)
  (F = Mnemi for Mlei)

Reports per species and set: n genes, n genes also in the h5ad (= what SATURN can use), dim, non-finite
values; whether all four "usable" gene sets are identical (FAIL otherwise, exit code 1).
Sanity checks:
  * median / p10 cosine(avg, longest) per model (ESM-1b 10-27: ~1.000 / 0.98-0.999; Hsap lowest)
  * ESMC_avg vs the old unrestricted averaged ESMC (03_embeddings_ESMC600M_postnorm): genes in both, and
    max |difference| for genes whose isoform set did not change (Dmel/Hsap should be 0: already restricted)
  * md5 of each .pt (recorded for provenance)

Usage (saturn_protein env; a few minutes; Hsap h5ad var_names read in backed mode):
  python 05y_check_gene_sets.py [--base /scratch/dark_genes/SATURN_Mnemi] [--out_tsv path]
"""
import argparse
import hashlib
import sys
from pathlib import Path

import anndata as ad
import torch

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--base", default="/scratch/dark_genes/SATURN_Mnemi")
ap.add_argument("--esmc_avg_dir", default="03_embeddings_ESMC600M_postnorm_avgmatched")
ap.add_argument("--out_tsv", default="")
args = ap.parse_args()
BASE = Path(args.base)
OUT_TSV = Path(args.out_tsv) if args.out_tsv else BASE / args.esmc_avg_dir / "05y_gene_set_check.tsv"
SPECIES = ["Mlei", "Cgig", "Crob", "Drer", "Dmel", "Hsap"]
FOLDER = {"Mlei": "Mnemi"}


def paths(sp):
    f = FOLDER.get(sp, sp)
    return {
        "ESM1b_avg": BASE / "03_embeddings" / sp / f"{sp}_gene_embeddings.pt",
        "ESMC_avg": BASE / args.esmc_avg_dir / f / f"{f}_gene_embeddings.pt",
        "ESM1b_longest": BASE / "03_embeddings_ESM1b_longest" / f / f"{f}_gene_embeddings.pt",
        "ESMC_longest": BASE / "03_embeddings_ESMC600M_postnorm_longest" / f / f"{f}_gene_embeddings.pt",
    }


def md5(p, chunk=1 << 22):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def cos_summary(a, b):
    genes = sorted(a.keys() & b.keys())
    if not genes:
        return float("nan"), float("nan")
    A = torch.stack([a[g].float() for g in genes])
    B = torch.stack([b[g].float() for g in genes])
    c = torch.nn.functional.cosine_similarity(A, B, dim=1)
    return c.median().item(), torch.quantile(c, 0.10).item()


fail = False
rows = []
print(f"=== 05y gene-set check  base={BASE} ===")
for sp in SPECIES:
    print(f"\n--- {sp} ---")
    h5ad = BASE / "02_processed_h5ad" / f"{sp}_saturn.h5ad"
    a = ad.read_h5ad(h5ad, backed="r")
    var = set(map(str, a.var_names))
    a.file.close()
    emb, usable = {}, {}
    for name, p in paths(sp).items():
        if not p.exists():
            print(f"  {name:14s} MISSING {p}")
            fail = True
            continue
        d = torch.load(p, map_location="cpu")
        emb[name] = d
        usable[name] = set(d.keys()) & var
        X = torch.stack([v.float().flatten() for v in d.values()])
        nonfinite = int((~torch.isfinite(X)).any(dim=1).sum())
        dig = md5(p)
        print(f"  {name:14s} genes {len(d):>7,}  in h5ad {len(usable[name]):>7,}  dim {X.shape[1]:>5}  "
              f"non-finite genes {nonfinite}  md5 {dig}")
        if nonfinite:
            fail = True
        rows.append((sp, name, len(d), len(usable[name]), X.shape[1], nonfinite, dig, str(p)))
    if len(usable) == 4:
        ref = usable["ESM1b_avg"]
        same = all(u == ref for u in usable.values())
        print(f"  usable gene sets identical across all 4: {'YES' if same else 'NO  <-- FAIL'}")
        if not same:
            fail = True
            for name, u in usable.items():
                if u != ref:
                    print(f"    {name}: only here {len(u - ref):,} (e.g. {sorted(u - ref)[:5]}); "
                          f"missing {len(ref - u):,} (e.g. {sorted(ref - u)[:5]})")
        for m in ("ESM1b", "ESMC"):
            med, p10 = cos_summary(emb[f"{m}_avg"], emb[f"{m}_longest"])
            print(f"  cos({m} avg, longest) median {med:.4f}  p10 {p10:.4f}")
    old = BASE / "03_embeddings_ESMC600M_postnorm" / FOLDER.get(sp, sp) / f"{FOLDER.get(sp, sp)}_gene_embeddings.pt"
    if old.exists() and "ESMC_avg" in emb:
        o = torch.load(old, map_location="cpu")
        both = sorted(o.keys() & emb["ESMC_avg"].keys())
        diffs = torch.stack([(o[g].float() - emb["ESMC_avg"][g].float()).abs().max() for g in both]) if both else None
        n_same = int((diffs < 1e-6).sum()) if diffs is not None else 0
        print(f"  vs old averaged ESMC ({len(o):,} genes): in both {len(both):,}; identical vectors {n_same:,}; "
              f"changed {len(both) - n_same:,} (isoform set restricted)")

OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_TSV, "w") as fh:
    fh.write("species\tset\tn_genes\tn_in_h5ad\tdim\tnonfinite_genes\tmd5\tpath\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
print(f"\nwritten -> {OUT_TSV}")
print("RESULT:", "FAIL" if fail else "PASS (gene sets identical, all finite)")
sys.exit(1 if fail else 0)
