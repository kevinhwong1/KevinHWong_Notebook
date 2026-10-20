#!/usr/bin/env python3
"""
05q_convert_esmc_postnorm.py   (2026-10-19)

Converts Ingrid's second ESMC-600M export (FANTASIA v4, "postnorm_maxlen2000", delivered 2026-10-16/19)
into SATURN-ready {gene: torch.Tensor} .pt files, one per species. Companion to
05_convert_esmc_embeddings.py, which converted the first (pre-norm, v4) export; that script is unchanged
so the Oct 2-9 runs stay reproducible.

What is the same as 05 (so the two ESMC exports are compared like for like):
  * protein -> gene bridges for Mlei/Cgig/Crob/Drer are imported from 05_convert_esmc_embeddings.py
    (one source of truth); Mlei uses the longest-isoform proteome, as in 05
  * several proteins per gene -> their embeddings are averaged (torch.mean), as SATURN's own
    convert_protein_embeddings_to_gene_embeddings.py and the ESM-1b run do
  * output layout 03_embeddings_ESMC600M_postnorm/{F}/{F}_gene_embeddings.pt, F = "Mnemi" for Mlei
    (the FOLDER convention 05f / 05i already use)

New:
  * Dmel and Hsap (Ensembl r113 pep.all, all isoforms). Ingrid's mapping CSVs (Dmel_to_pep_all.csv,
    Hsap_to_pep_all.csv: id, full_header, embedding_name) carry the Ensembl gene ID in full_header
    ("gene:FBgn..." / "gene:ENSG....<version>"); the version is stripped, exactly as 02_embed.sh
    step1_ensembl_gene_id.py does, so keys match the h5ad var_names and the ESM-1b .pt files.
  * Dmel/Hsap are RESTRICTED (default) to the proteins in 03_embeddings/{SP}/{SP}_filtered.faa, i.e. the
    exact protein set embedded with ESM-1b (mapped to an h5ad gene; internal-stop proteins removed:
    21 Hsap proteins / 18 genes, 0 Dmel). Ingrid embedded all of pep.all (Hsap 123,845, Dmel 30,799), so
    without this the two models would average different isoform sets. --no_restrict turns it off.
  * --compare_esm1b reports, per species, how the gene set compares with the ESM-1b .pt file.

Truncation differs between models and is NOT harmonised here: ESMC postnorm = first 2,000 aa,
ESM-1b = first 1,022 aa (1-4% of proteins are > 2,000 aa; more are > 1,022). Note it in Methods.

Inputs (copy Ingrid's folder to Pegasus first; see the plan doc):
  RAW = /scratch/dark_genes/SATURN_Mnemi/raw_embeddings/ESMC600M_postnorm/
        {Cgig,Crob,Drer,Mnemi,Dmel,Hsap}/*_postnorm_maxlen2000.h5, proteome_h5_id_mapping.csv,
        Dmel_to_pep_all.csv, Hsap_to_pep_all.csv

Usage (saturn_protein env; minutes, login node is fine for the 4 species, use bigmem for Hsap):
  python 05q_convert_esmc_postnorm.py                              # Mlei Cgig Crob Drer (4-species check)
  python 05q_convert_esmc_postnorm.py --species Dmel Hsap --compare_esm1b
  python 05q_convert_esmc_postnorm.py --species all --compare_esm1b
"""
import argparse
import csv
import importlib.util
import re
from collections import defaultdict
from pathlib import Path

import h5py
import torch

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--base", default="/scratch/dark_genes/SATURN_Mnemi")
ap.add_argument("--raw", default="", help="default: <base>/raw_embeddings/ESMC600M_postnorm")
ap.add_argument("--out", default="", help="default: <base>/03_embeddings_ESMC600M_postnorm")
ap.add_argument("--species", nargs="+", default=["Mlei", "Cgig", "Crob", "Drer"],
                help="any of Mlei Cgig Crob Drer Dmel Hsap, or 'all'")
ap.add_argument("--no_restrict", action="store_true",
                help="Dmel/Hsap: use every pep.all protein instead of only those in {SP}_filtered.faa")
ap.add_argument("--compare_esm1b", action="store_true",
                help="report overlap with 03_embeddings/{SP}/{SP}_gene_embeddings.pt")
args = ap.parse_args()

BASE = Path(args.base)
RAW = Path(args.raw) if args.raw else BASE / "raw_embeddings" / "ESMC600M_postnorm"
OUT = Path(args.out) if args.out else BASE / "03_embeddings_ESMC600M_postnorm"
ALL_SP = ["Mlei", "Cgig", "Crob", "Drer", "Dmel", "Hsap"]
SPECIES = ALL_SP if args.species == ["all"] else args.species
bad = [s for s in SPECIES if s not in ALL_SP]
if bad:
    raise SystemExit(f"unknown species {bad}")

FOLDER = {"Mlei": "Mnemi"}                       # embedding-folder name (05f/05i convention)
EMBED_SUFFIX = "type_5/layer_0/embedding"        # same in both FANTASIA exports
H5 = {
    "Mlei": RAW / "Mnemi" / "Mnemi-longestisoform_esmc-600m_postnorm_maxlen2000.h5",
    "Cgig": RAW / "Cgig" / "Cgig_esmc-600m_postnorm_maxlen2000.h5",
    "Crob": RAW / "Crob" / "Crob_esmc-600m_postnorm_maxlen2000.h5",
    "Drer": RAW / "Drer" / "Drer_esmc-600m_postnorm_maxlen2000.h5",
    "Dmel": RAW / "Dmel" / "Dmel_esmc_postnorm_maxlen2000.h5",
    "Hsap": RAW / "Hsap" / "Hsap_esmc_postnorm_maxlen2000.h5",
}
PROTEOME_NAME = {  # 'proteome' column of proteome_h5_id_mapping.csv (same file as for the v4 export)
    "Mlei": "Mnemi-longestisoform_proteome.fasta", "Cgig": "Cgig_proteome.fasta",
    "Crob": "Crob_proteome.fasta", "Drer": "Drer_proteome.fasta",
}


def load_bridges():
    """Import the Mlei/Cgig/Crob/Drer protein->gene bridges from 05 (it only defines functions/constants)."""
    p = Path(__file__).resolve().parent / "05_convert_esmc_embeddings.py"
    spec = importlib.util.spec_from_file_location("conv05", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return {"Mlei": m.build_bridge_mlei, "Cgig": m.build_bridge_cgig,
            "Crob": m.build_bridge_crob, "Drer": m.build_bridge_drer}


def rows_fantasia(sp):
    """protein id -> (h5 group path, gene) for the 4 original species."""
    rows = {}
    with open(RAW / "proteome_h5_id_mapping.csv", newline="") as f:
        for r in csv.DictReader(f):
            if r["proteome"] == PROTEOME_NAME[sp]:
                rows[r["fasta_id"]] = r["hdf5_group_path"]
    print(f"  mapping rows          : {len(rows):,}")
    bridge = load_bridges()[sp](rows.keys())
    return {pid: (rows[pid], g) for pid, g in bridge.items() if pid in rows}, len(rows)


def filtered_ids(sp):
    fa = BASE / "03_embeddings" / sp / f"{sp}_filtered.faa"
    if not fa.exists():
        raise SystemExit(f"{fa} not found (needed to match the ESM-1b protein set; or pass --no_restrict)")
    with open(fa) as f:
        return {line[1:].split()[0] for line in f if line.startswith(">")}


def rows_ensembl(sp):
    """protein id -> (h5 group path, gene) for Dmel/Hsap via Ingrid's *_to_pep_all.csv."""
    gre = re.compile(r"gene:([^ .]+)")           # identical to 02_embed.sh step1_ensembl_gene_id.py
    out, n = {}, 0
    with open(RAW / f"{sp}_to_pep_all.csv", newline="") as f:
        for r in csv.DictReader(f):
            n += 1
            m = gre.search(r["full_header"])
            if m:
                out[r["id"]] = ("/" + r["embedding_name"].lstrip("/"), m.group(1))
    print(f"  mapping rows          : {n:,}  ({len(out):,} with a gene: field)")
    if not args.no_restrict:
        keep = filtered_ids(sp)
        missing = keep - out.keys()
        out = {pid: v for pid, v in out.items() if pid in keep}
        print(f"  restricted to ESM-1b protein set ({sp}_filtered.faa: {len(keep):,} proteins): "
              f"{len(out):,} kept; {len(missing):,} filtered proteins absent from Ingrid's mapping")
        if missing:
            print(f"    e.g. {sorted(missing)[:5]}")
    return out, n


def convert(sp):
    print(f"\n{'=' * 60}\n  {sp}  ->  {H5[sp].name}\n{'=' * 60}")
    prot, n_rows = rows_ensembl(sp) if sp in ("Dmel", "Hsap") else rows_fantasia(sp)
    by_gene = defaultdict(list)
    for pid, (grp, g) in prot.items():
        by_gene[g].append(grp)
    print(f"  proteins bridged      : {len(prot):,} -> {len(by_gene):,} genes "
          f"({sum(len(v) > 1 for v in by_gene.values()):,} genes average >1 protein)")
    res, n_missing = {}, 0
    with h5py.File(H5[sp], "r") as h:
        for g, grps in by_gene.items():
            vecs = []
            for grp in grps:
                ds = f"{grp}/{EMBED_SUFFIX}"
                if ds in h:
                    vecs.append(torch.tensor(h[ds][()], dtype=torch.float32))
                else:
                    n_missing += 1
            if vecs:
                res[g] = torch.mean(torch.stack(vecs), dim=0) if len(vecs) > 1 else vecs[0]
    if n_missing:
        print(f"  WARNING: {n_missing:,} mapped proteins have no dataset in the .h5")
    norms = torch.stack(list(res.values())).norm(dim=1)
    print(f"  gene embeddings       : {len(res):,}  dim {next(iter(res.values())).numel()}  "
          f"norm min/median/max {norms.min():.3f} / {norms.median():.3f} / {norms.max():.3f}")
    f = FOLDER.get(sp, sp)
    dst = OUT / f / f"{f}_gene_embeddings.pt"
    dst.parent.mkdir(parents=True, exist_ok=True)
    torch.save(res, dst)
    print(f"  written -> {dst}")
    if args.compare_esm1b:
        p = BASE / "03_embeddings" / sp / f"{sp}_gene_embeddings.pt"
        if p.exists():
            e = set(torch.load(p, map_location="cpu").keys())
            el = {k.lower() for k in e}
            both = e & res.keys()
            ci = {k for k in res if k not in e and k.lower() in el}
            print(f"  vs ESM-1b ({len(e):,} genes): in both {len(both):,}; case-only differences {len(ci):,}; "
                  f"ESM-1b only {len(e) - len(both) - len(ci):,}; ESMC only {len(res) - len(both) - len(ci):,}")
            only_e = sorted(e - res.keys() - {k for k in e if k.lower() in {x.lower() for x in ci}})
            if only_e:
                print(f"    ESM-1b-only examples: {only_e[:8]}")
        else:
            print(f"  (no ESM-1b file at {p})")


print(f"=== ESMC-600M postnorm -> SATURN .pt ===\nRAW {RAW}\nOUT {OUT}\nrestrict Dmel/Hsap to ESM-1b proteins: "
      f"{not args.no_restrict}")
for sp in SPECIES:
    convert(sp)
print("\nDone.")
