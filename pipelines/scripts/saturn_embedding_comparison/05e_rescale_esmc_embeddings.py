#!/usr/bin/env python3
"""
05e_rescale_esmc_embeddings.py

Writes a length-normalized copy of the ESMC-600M gene embeddings for SATURN.

Why: SATURN's make_centroids() runs Euclidean KMeans on the raw stacked gene
embeddings (normalize=False; no rescaling in data/gene_embeddings.py). The
FANTASIA ESMC vectors have L2 norms of ~1,700-8,900 (Drer systematically
lower), whereas ESM-1b norms are ~16-20. With ESMC, vector length therefore
drives macrogene initialisation as much as direction.

What: every ESMC gene vector v -> v / ||v|| * TARGET, where TARGET is the mean
L2 norm of all ESM-1b gene vectors across the same species (computed here, not
hard-coded). Directions are unchanged; only lengths are equalised. The
original ESMC files in 03_embeddings_ESMC600M/ are not modified.

Output: 03_embeddings_ESMC600M_rescaled/<folder>/<folder>_gene_embeddings.pt
(same folder names as the source, incl. "Mnemi" for Mlei)
"""
from pathlib import Path
import torch

BASE = Path("/scratch/dark_genes/SATURN_Mnemi")
SRC = BASE / "03_embeddings_ESMC600M"
DST = BASE / "03_embeddings_ESMC600M_rescaled"
ESM1B = BASE / "03_embeddings"
# species label -> ESMC folder name (FANTASIA used "Mnemi" for Mlei)
ESMC_FOLDERS = {"Mlei": "Mnemi", "Cgig": "Cgig", "Crob": "Crob", "Drer": "Drer"}


def stack(d):
    return torch.stack([v.float().flatten() for v in d.values()])


# Target norm = mean ESM-1b gene-vector norm over the same four species
norms = []
for sp in ESMC_FOLDERS:
    d = torch.load(ESM1B / sp / f"{sp}_gene_embeddings.pt", map_location="cpu")
    norms.append(stack(d).norm(dim=1))
norms = torch.cat(norms)
TARGET = norms.mean().item()
print(f"ESM-1b gene-vector norm: mean {TARGET:.4f}  sd {norms.std().item():.4f}  "
      f"(n={norms.numel():,}) -> TARGET = {TARGET:.4f}")

for sp, folder in ESMC_FOLDERS.items():
    src = SRC / folder / f"{folder}_gene_embeddings.pt"
    d = torch.load(src, map_location="cpu")
    out, n_zero = {}, 0
    before = stack(d).norm(dim=1)
    for g, v in d.items():
        v = v.float()
        n = v.norm()
        if n == 0:
            n_zero += 1
            out[g] = v.clone()
            continue
        out[g] = v / n * TARGET
    after = stack(out).norm(dim=1)
    dst = DST / folder / f"{folder}_gene_embeddings.pt"
    dst.parent.mkdir(parents=True, exist_ok=True)
    torch.save(out, dst)
    print(f"{sp:5s} ({folder}): {len(out):,} genes | norm before mean {before.mean():.1f} "
          f"[{before.min():.1f}-{before.max():.1f}] -> after mean {after.mean():.4f} "
          f"[{after.min():.4f}-{after.max():.4f}] | zero-norm kept as-is: {n_zero}")
    print(f"      -> {dst}")
print("Done.")
