#!/usr/bin/env python3
"""
05_convert_esmc_embeddings.py

Converts FANTASIA v4 ESMC-600M per-protein embeddings (HDF5 + mapping CSV)
into SATURN-ready {gene_symbol: torch.Tensor} .pt files, one per species —
using the SAME gene-ID bridging strategies already validated in
check_proteome_overlap.py for the original ESM-1b run:

    Mlei  -> direct               (g_XXXX -> g-XXXX)
    Cgig  -> gtf_loc_to_xp         (reversed: XP/NP accession -> LOC gene_id)
    Crob  -> direct_strip_suffix   (KY21_ChrN_NNNN_... -> KY21.ChrN.NNNN)
    Drer  -> gtf_symbol_to_ensdarg (reversed: ENSDARP -> ENSDARG -> gene symbol)

Species covered here: Cgig, Crob, Drer, Mnemi (longest-isoform set only).
Spur is intentionally excluded (poor SATURN integration -> dropped from
the project; see project notes).

Usage:
    python 05_convert_esmc_embeddings.py                       # all 4 species
    python 05_convert_esmc_embeddings.py --species Crob Drer   # just these
    python 05_convert_esmc_embeddings.py --check-coverage      # also diff
                                                                 # against each
                                                                 # species' h5ad
                                                                 # var_names

Run inside the saturn_protein conda env (needs h5py, torch, and scanpy
only if --check-coverage is used).
"""

import re
import csv
import argparse
import warnings
from pathlib import Path
from collections import defaultdict

import h5py
import torch

# ============================================================
# Paths — ADJUST THESE if your raw_embeddings / proteome layout
# differs from what we set up together. Everything else in the
# script should not need touching.
# ============================================================
BASE = Path("/scratch/dark_genes/SATURN_Mnemi")
RAW = BASE / "raw_embeddings" / "ESMC600M_FANTASIA"
OUT_DIR = BASE / "03_embeddings_ESMC600M"

MAPPING_CSV = RAW / "proteome_h5_id_mapping.csv"

H5_FILES = {
    "Cgig": RAW / "Cgig" / "Cgig_v4_embeddings_ESMC-600M.h5",
    "Crob": RAW / "Crob" / "Crob_v4_embeddings_ESMC-600M.h5",
    "Drer": RAW / "Drer" / "Drer_v4_embeddings_ESM3c.h5",
    "Mnemi": RAW / "Mnemi" / "Mnemi-longestisoform_v4_embeddings_ESM3c.h5",
}

PROTEOME_NAMES = {
    "Cgig": "Cgig_proteome.fasta",
    "Crob": "Crob_proteome.fasta",
    "Drer": "Drer_proteome.fasta",
    "Mnemi": "Mnemi-longestisoform_proteome.fasta",
}

PROTEOMES_DIR = BASE / "01_proteomes"

GTF_PATHS = {
    "Cgig": PROTEOMES_DIR / "Cgig" / "Cgig.genomic.gtf",
    "Drer": PROTEOMES_DIR / "Drer" / "Drer.genomic.gtf",
}
FASTA_PATHS = {
    "Cgig": PROTEOMES_DIR / "Cgig" / "Cgig.protein.faa",
    "Drer": PROTEOMES_DIR / "Drer" / "Drer.protein.faa",
}

# Only used with --check-coverage, to diff against the real var_names.
# Adjust paths if these h5ad files live somewhere else now.
H5AD_PATHS = {
    "Cgig": BASE / "02_processed_h5ad" / "Cgig_saturn.h5ad",
    "Crob": BASE / "02_processed_h5ad" / "Crob_saturn.h5ad",
    "Drer": BASE / "02_processed_h5ad" / "Drer_saturn.h5ad",
    "Mnemi": BASE / "02_processed_h5ad" / "Mlei_saturn.h5ad",
}

# Confirmed identical across all 5 FANTASIA .h5 files during inspection:
# {group}/type_5/layer_0/embedding, shape (1152,), float64.
EMBED_DATASET_SUFFIX = "type_5/layer_0/embedding"


# ============================================================
# Per-species fasta_id -> gene_symbol bridges
# (mirrors the strategies validated in check_proteome_overlap.py,
#  but run in reverse: we start from a protein accession and need
#  the gene symbol used in the h5ad var_names, not the other way)
# ============================================================

def build_bridge_mlei(fasta_ids):
    """g16414_anno2.3794_t / g10000_anno1.g10373.t1 -> g-16414 / g-10000"""
    bridge = {}
    pat = re.compile(r'^g(\d+)_')
    for fid in fasta_ids:
        m = pat.match(fid)
        if m:
            bridge[fid] = f"g-{m.group(1)}"
    return bridge


def build_bridge_crob(fasta_ids):
    """KY21_Chr10_1000_v1_SL1-1 -> KY21.Chr10.1000

    This is the same direct_strip_suffix strategy used for the ESM-1b run
    (regex ^(KY21\\.[^.]+\\.\\d+) against dot-separated Ghost DB headers),
    adapted because FANTASIA's HDF5 group names use underscores instead
    of dots.
    """
    bridge = {}
    pat = re.compile(r'^(KY21_[^_]+_\d+)')
    for fid in fasta_ids:
        m = pat.match(fid)
        if m:
            bridge[fid] = m.group(1).replace("_", ".")
    return bridge


def build_bridge_cgig(fasta_ids):
    """NP_/XP_ protein accession -> LOC gene_id, via GTF CDS attributes.

    Same source data as the gtf_loc_to_xp strategy (gene_id/protein_id on
    CDS lines), just inverted: there we needed gene_id -> protein_id to
    check FASTA coverage; here we need protein_id -> gene_id directly.
    """
    xp_to_loc = {}
    with open(GTF_PATHS["Cgig"]) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attr = parts[8]
            m_loc = re.search(r'gene_id "([^"]+)"', attr)
            m_xp = re.search(r'protein_id "([^"]+)"', attr)
            if m_loc and m_xp:
                xp_to_loc[m_xp.group(1)] = m_loc.group(1)
    return {fid: xp_to_loc[fid] for fid in fasta_ids if fid in xp_to_loc}


def build_bridge_drer(fasta_ids):
    """ENSDARP accession -> gene symbol.

    Mirrors gtf_symbol_to_ensdarg, run in reverse:
      1. GTF 'gene' lines give gene_id (ENSDARG) -> gene_name (symbol)
      2. protein.faa headers carry 'gene:ENSDARG...' per ENSDARP accession
      3. chain protein -> ENSDARG -> symbol
    """
    ensdarg_to_symbol = {}
    with open(GTF_PATHS["Drer"]) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attr = parts[8]
            m_id = re.search(r'gene_id "([^"]+)"', attr)
            m_sym = re.search(r'gene_name "([^"]+)"', attr)
            if m_id and m_sym:
                ensdarg_to_symbol[m_id.group(1)] = m_sym.group(1)

    ensdarp_to_ensdarg = {}
    gene_re = re.compile(r'gene:(ENSDARG\d+)')
    with open(FASTA_PATHS["Drer"]) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            acc = line[1:].split()[0].strip()
            m = gene_re.search(line)
            if m:
                ensdarp_to_ensdarg[acc] = m.group(1)

    bridge = {}
    for fid in fasta_ids:
        ensdarg = ensdarp_to_ensdarg.get(fid)
        if ensdarg and ensdarg in ensdarg_to_symbol:
            bridge[fid] = ensdarg_to_symbol[ensdarg]
    return bridge


BRIDGE_BUILDERS = {
    "Mnemi": build_bridge_mlei,
    "Crob": build_bridge_crob,
    "Cgig": build_bridge_cgig,
    "Drer": build_bridge_drer,
}


# ============================================================
# Core conversion
# ============================================================

def load_mapping_rows(species):
    """fasta_id -> hdf5_group_path, filtered to this species' proteome rows."""
    proteome = PROTEOME_NAMES[species]
    rows = {}
    with open(MAPPING_CSV, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["proteome"] == proteome:
                rows[row["fasta_id"]] = row["hdf5_group_path"]
    return rows


def convert_species(species, check_coverage=False):
    print(f"\n{'=' * 60}")
    print(f"  {species}")
    print(f"{'=' * 60}")

    rows = load_mapping_rows(species)
    print(f"  Mapping CSV rows      : {len(rows):,}")

    bridge = BRIDGE_BUILDERS[species](rows.keys())
    mapped = {fid: gsym for fid, gsym in bridge.items() if fid in rows}
    pct = 100 * len(mapped) / len(rows) if rows else 0
    print(f"  Bridged to gene symbol: {len(mapped):,}  ({pct:.1f}%)")

    # Multiple fasta_ids can map to the same gene symbol (isoform variants,
    # e.g. Crob's nonSL-tagged splice variants, or multiple RefSeq XP_
    # transcripts per Cgig LOC). SATURN's own
    # convert_protein_embeddings_to_gene_embeddings.py handles this by
    # averaging every isoform's embedding into one gene-level vector
    # (torch.mean over torch.stack of the per-protein embeddings) — that's
    # also what your original ESM-1b run did for every species whose
    # strategy routes through that script, so we match it here rather than
    # arbitrarily keeping one isoform.
    gsym_to_fids = defaultdict(list)
    for fid, gsym in mapped.items():
        gsym_to_fids[gsym].append(fid)
    collisions = {g: f for g, f in gsym_to_fids.items() if len(f) > 1}
    if collisions:
        print(f"  {len(collisions)} gene symbols matched by >1 fasta_id "
              f"-- averaging their embeddings (matches SATURN's own "
              f"convert_protein_embeddings_to_gene_embeddings.py)")
        for g, fids in list(collisions.items())[:3]:
            print(f"    {g}: {fids}")

    h5_path = H5_FILES[species]
    result = {}
    skipped_missing_dataset = 0
    with h5py.File(h5_path, "r") as hf:
        for gsym, fids in gsym_to_fids.items():
            vecs = []
            for fid in fids:
                group_path = rows[fid]
                ds_path = f"{group_path}/{EMBED_DATASET_SUFFIX}"
                if ds_path not in hf:
                    skipped_missing_dataset += 1
                    continue
                vecs.append(torch.tensor(hf[ds_path][()], dtype=torch.float32))
            if not vecs:
                continue
            result[gsym] = torch.mean(torch.stack(vecs), dim=0) if len(vecs) > 1 else vecs[0]

    if skipped_missing_dataset:
        print(f"  WARNING: {skipped_missing_dataset} mapped IDs had no "
              f"embedding dataset in the .h5 file")

    print(f"  Final gene embeddings : {len(result):,}")

    out_path = OUT_DIR / species / f"{species}_gene_embeddings.pt"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(result, out_path)
    print(f"  Written -> {out_path}")

    if check_coverage:
        h5ad_path = H5AD_PATHS.get(species)
        if h5ad_path and h5ad_path.exists():
            warnings.filterwarnings("ignore")
            import scanpy as sc
            ad = sc.read_h5ad(str(h5ad_path), backed="r")
            var_genes = set(ad.var_names.astype(str))
            covered = var_genes & result.keys()
            cov_pct = 100 * len(covered) / len(var_genes) if var_genes else 0
            print(f"  h5ad var_names        : {len(var_genes):,}")
            print(f"  Covered by embeddings : {len(covered):,}  ({cov_pct:.1f}%)")
            missing = sorted(var_genes - result.keys())[:10]
            if missing:
                print(f"  Sample missing genes  : {missing}")
        else:
            print(f"  (--check-coverage requested but h5ad not found at "
                  f"{h5ad_path} -- skipping)")

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--species", nargs="+", default=list(BRIDGE_BUILDERS.keys()),
        choices=list(BRIDGE_BUILDERS.keys()),
        help="Which species to convert (default: all 4)",
    )
    parser.add_argument(
        "--check-coverage", action="store_true",
        help="Also compare the resulting gene set against each species' "
             "h5ad var_names and report coverage %%.",
    )
    args = parser.parse_args()

    print("=== ESMC-600M -> SATURN .pt conversion ===")
    for sp in args.species:
        convert_species(sp, check_coverage=args.check_coverage)

    print(f"\n{'=' * 60}")
    print("Done.")


if __name__ == "__main__":
    main()
