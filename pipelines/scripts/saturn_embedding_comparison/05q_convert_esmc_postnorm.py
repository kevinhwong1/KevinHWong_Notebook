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
  * --longest (added 2026-10-27): ONE protein per gene instead of averaging isoforms. For each gene the
    protein is the one ESM-1b's longest-isoform rebuild used (05t_esm1b_longest_isoform.py: longest in
    03_embeddings/{SP}/{SP}_filtered.faa, ties -> alphabetically first ID), matched to Ingrid's protein IDs
    (exact, version-stripped, or '.'/'_'-insensitive). Genes without that match: longest ESMC protein by
    length from {SP}_filtered.faa / 01_proteomes FASTA; if no length is known, isoforms are averaged and
    counted. Default output 03_embeddings_ESMC600M_postnorm_longest/. Writes esmc_longest_choice_{SP}.tsv.
  * --match_esm1b (added 2026-10-29): isoforms AVERAGED, but over exactly the proteins ESM-1b averaged
    (03_embeddings/{SP}/{SP}_gene_to_protein.filtered.json), in ALL six species, and only for genes in the
    ESM-1b .pt (keys renamed to ESM-1b's spelling when they differ only by case). Without it, the averaged
    ESMC files differ from ESM-1b in gene set (Drer 19,937 vs 14,179 genes; Cgig extra Cg_lysoz3) and, for
    Mlei/Cgig/Crob/Drer, possibly in isoform set. Matching per protein uses id_keys (exact, version-stripped,
    '.'/'_'-insensitive). Default output 03_embeddings_ESMC600M_postnorm_avgmatched/. Writes
    esmc_avgmatch_{SP}.tsv (per gene: ESM-1b proteins, ESMC proteins, matched, how). Not combinable with
    --longest (which already uses ESM-1b's protein).

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
  python 05q_convert_esmc_postnorm.py --species all --longest --compare_esm1b   # one protein per gene
  python 05q_convert_esmc_postnorm.py --species all --match_esm1b --compare_esm1b  # averaged, ESM-1b's isoforms
"""
import argparse
import csv
import importlib.util
import json
import re
from collections import Counter, defaultdict
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
ap.add_argument("--longest", action="store_true",
                help="one protein per gene (the one ESM-1b's longest-isoform rebuild used) instead of averaging")
ap.add_argument("--match_esm1b", action="store_true",
                help="average over exactly ESM-1b's isoforms per gene, ESM-1b gene set only (all species)")
args = ap.parse_args()
if args.longest and args.match_esm1b:
    raise SystemExit("--longest and --match_esm1b are exclusive (--longest already uses ESM-1b's protein)")

BASE = Path(args.base)
RAW = Path(args.raw) if args.raw else BASE / "raw_embeddings" / "ESMC600M_postnorm"
OUT = Path(args.out) if args.out else BASE / ("03_embeddings_ESMC600M_postnorm_longest" if args.longest
                                              else "03_embeddings_ESMC600M_postnorm_avgmatched" if args.match_esm1b
                                              else "03_embeddings_ESMC600M_postnorm")
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


PROTEOME_FASTA = {  # extra length sources for proteins outside {SP}_filtered.faa (same files as run_config.sh)
    "Cgig": BASE / "01_proteomes" / "Cgig" / "Cgig.protein.faa",
    "Crob": BASE / "01_proteomes" / "Crob" / "Crob_KY21_ghost.faa",
    "Drer": BASE / "01_proteomes" / "Drer" / "Drer.protein.faa",
}


def id_keys(pid):
    """Match keys for one protein ID, strictest first: exact; version stripped; '.'/'_' removed."""
    k1 = pid
    k2 = re.sub(r"\.\d+$", "", pid)
    k3 = re.sub(r"[._]", "", k2).lower()
    return [k1, k2, k3]


def fasta_lengths(path, lens):
    if not Path(path).exists():
        return
    cur, n = None, 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    for k in id_keys(cur):
                        lens.setdefault(k, n)
                cur, n = line[1:].split()[0].strip(), 0
            else:
                n += len(line.strip().rstrip("*"))
    if cur is not None:
        for k in id_keys(cur):
            lens.setdefault(k, n)


def esm1b_choice(sp, lens_esm):
    """gene -> protein ESM-1b's longest-isoform rebuild used (same rule as 05t)."""
    ed = BASE / "03_embeddings" / sp
    gmap = json.load(open(ed / f"{sp}_gene_to_protein.filtered.json"))
    return {g: sorted(ps, key=lambda p: (-lens_esm.get(p, -1), p))[0] for g, ps in gmap.items() if ps}


def pick_longest(sp, by_gene):
    """by_gene: gene -> [(pid, grp)]  ->  gene -> [(pid, grp)] with one entry (or all, if averaging is the fallback)."""
    lens_esm = {}
    with open(BASE / "03_embeddings" / sp / f"{sp}_filtered.faa") as f:
        cur, n = None, 0
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    lens_esm[cur] = n
                cur, n = line[1:].split()[0].strip(), 0
            else:
                n += len(line.strip().rstrip("*"))
        if cur is not None:
            lens_esm[cur] = n
    choice = esm1b_choice(sp, lens_esm)
    lens = {}
    for pid, n in lens_esm.items():
        for k in id_keys(pid):
            lens.setdefault(k, n)
    if sp in PROTEOME_FASTA:
        fasta_lengths(PROTEOME_FASTA[sp], lens)
    out, rows = {}, []
    cnt = {"single": 0, "matched_esm1b": 0, "by_length": 0, "avg_fallback": 0}
    for g, items in by_gene.items():
        if len(items) == 1:
            out[g] = items; cnt["single"] += 1
            rows.append((g, len(items), items[0][0], "single"))
            continue
        hit = None
        if g in choice:
            want = id_keys(choice[g])
            for level in range(3):
                cand = [it for it in items if id_keys(it[0])[level] == want[level]]
                if len(cand) == 1:
                    hit = cand[0]; break
        if hit is not None:
            out[g] = [hit]; cnt["matched_esm1b"] += 1
            rows.append((g, len(items), hit[0], "matched_esm1b"))
            continue
        L = []
        for it in items:
            n = next((lens[k] for k in id_keys(it[0]) if k in lens), None)
            L.append(n)
        if all(n is not None for n in L):
            best = sorted(zip(items, L), key=lambda x: (-x[1], x[0][0]))[0][0]
            out[g] = [best]; cnt["by_length"] += 1
            rows.append((g, len(items), best[0], "by_length"))
        else:
            out[g] = items; cnt["avg_fallback"] += 1
            rows.append((g, len(items), "|".join(i[0] for i in items), "avg_fallback"))
    multi = len(by_gene) - cnt["single"]
    print(f"  longest-isoform choice: {cnt}  (multi-protein genes: {multi:,}; "
          f"matched to ESM-1b's protein: {100 * cnt['matched_esm1b'] / max(multi, 1):.1f}%)")
    f = FOLDER.get(sp, sp)
    (OUT / f).mkdir(parents=True, exist_ok=True)
    with open(OUT / f"esmc_longest_choice_{sp}.tsv", "w") as fh:
        fh.write("gene\tn_esmc_proteins\tchosen\thow\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    return out


def match_esm1b_isoforms(sp, by_gene):
    """gene -> [(pid, grp)]  ->  same, restricted to ESM-1b's gene set and, per gene, to ESM-1b's isoforms."""
    ed = BASE / "03_embeddings" / sp
    esm_genes = list(torch.load(ed / f"{sp}_gene_embeddings.pt", map_location="cpu").keys())
    gmap = json.load(open(ed / f"{sp}_gene_to_protein.filtered.json"))
    esm_set = set(esm_genes)
    lower = defaultdict(list)
    for g in esm_genes:
        lower[g.lower()].append(g)
    # 1. gene keys: exact, else a unique case-insensitive ESM-1b gene that ESMC lacks under its exact name
    renamed, n_case, dropped = {}, 0, []
    for g, items in by_gene.items():
        if g in esm_set:
            key = g
        elif len(lower.get(g.lower(), [])) == 1 and lower[g.lower()][0] not in by_gene:
            key = lower[g.lower()][0]
            n_case += 1
        else:
            dropped.append(g)
            continue
        renamed.setdefault(key, []).extend(items)
    # 2. isoforms: keep the ESMC proteins that match ESM-1b's list for that gene
    out, rows, cnt = {}, [], Counter()
    for g, items in renamed.items():
        want = gmap.get(g, [])
        got = []
        for w in want:
            wk = id_keys(w)
            for level in range(3):
                cand = [it for it in items if id_keys(it[0])[level] == wk[level]]
                if len(cand) == 1:
                    if cand[0] not in got:
                        got.append(cand[0])
                    break
        if not want:
            got, how = items, "no_esm1b_list_kept_all"
        elif len(got) == len(want):
            how = "all_matched"
        elif got:
            how = "partial"
        else:
            got, how = items, "none_matched_kept_all"
        cnt[how] += 1
        out[g] = got
        rows.append((g, len(want), len(items), len(got), how))
    missing = sorted(esm_set - out.keys())
    print(f"  match ESM-1b: genes kept {len(out):,} of ESM-1b {len(esm_set):,}; case-renamed {n_case:,}; "
          f"ESMC genes dropped (not in ESM-1b) {len(dropped):,}; ESM-1b genes missing in ESMC {len(missing):,}")
    if dropped:
        print(f"    dropped e.g. {sorted(dropped)[:8]}")
    if missing:
        print(f"    MISSING e.g. {missing[:8]}")
    print(f"  isoform matching      : {dict(cnt)}  (ESMC isoforms removed: "
          f"{sum(r[2] - r[3] for r in rows):,})")
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / f"esmc_avgmatch_{sp}.tsv", "w") as fh:
        fh.write("gene\tn_esm1b_proteins\tn_esmc_proteins\tn_used\thow\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    return out


def convert(sp):
    print(f"\n{'=' * 60}\n  {sp}  ->  {H5[sp].name}\n{'=' * 60}")
    prot, n_rows = rows_ensembl(sp) if sp in ("Dmel", "Hsap") else rows_fantasia(sp)
    by_gene = defaultdict(list)
    for pid, (grp, g) in prot.items():
        by_gene[g].append((pid, grp))
    print(f"  proteins bridged      : {len(prot):,} -> {len(by_gene):,} genes "
          f"({sum(len(v) > 1 for v in by_gene.values()):,} genes have >1 protein)")
    if args.longest:
        by_gene = pick_longest(sp, by_gene)
    if args.match_esm1b:
        by_gene = match_esm1b_isoforms(sp, by_gene)
    by_gene = {g: [grp for _, grp in items] for g, items in by_gene.items()}
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
      f"{not args.no_restrict}\none protein per gene (--longest): {args.longest}\n"
      f"average over ESM-1b's isoforms, ESM-1b gene set (--match_esm1b): {args.match_esm1b}")
for sp in SPECIES:
    convert(sp)
print("\nDone.")
