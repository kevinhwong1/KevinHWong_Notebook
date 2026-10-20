#!/usr/bin/env python3
"""
06a_orthofinder_prep.py  (runs on Pegasus; a few minutes, login node or bigmem)

Prepares OrthoFinder input as a second, de novo orthology reference for the SATURN benchmarks
(the first is eggNOG Metazoa OGs). One protein per gene, the LONGEST (zebrafish: primary-assembly proteins preferred over ALT_CTG copies), using the same proteome FASTAs and
the same protein -> gene bridges as the SATURN embeddings (check_proteome_overlap.py / 05_convert):
    Mlei  direct                 FASTA header = gene id (g-XXXX); fallback g12345_... -> g-12345
    Cgig  gtf_loc_to_xp          XP_/NP_ protein_id -> LOC gene_id   (GTF CDS lines)
    Crob  direct_strip_suffix    KY21.Chr1.123.v1.SL1-1 -> KY21.Chr1.123
    Drer  gtf_symbol_to_ensdarg  ENSDARP -> ENSDARG (FASTA header 'gene:') -> gene_name (GTF);
                                 genes with no gene_name are kept under their ENSDARG id (not SATURN genes,
                                 but part of the proteome OrthoFinder should see), EXCEPT unnamed genes on GRCz11
                                 alternate-haplotype contigs (ALT_CTG), which duplicate primary loci and are dropped
                                 (named ALT genes collapse with their primary copy by symbol)
    Dmel  ensembl_gene_id        Ensembl pep.all header 'gene:FBgn...' (version stripped) = h5ad var_names
    Hsap  ensembl_gene_id        Ensembl pep.all header 'gene:ENSG...'  (version stripped) = h5ad var_names;
                                 proteins on non-primary scaffolds are dropped: alternate haplotypes (HSCHR*) and
                                 fix/novel patches (HG*). In Ensembl r113 pep.all that is 11,512 proteins from
                                 3,337 genes, 98% of whose proteins share a symbol with a primary-assembly gene (they would
                                 otherwise look like paralogs). Unplaced/unlocalised primary scaffolds (KI*, GL*) are kept.
    Sros  gff_cds_locus          outgroup (choanoflagellate Salpingoeca rosetta, NCBI; 06a0_download_outgroup.sh):
                                 protein_id -> locus_tag (or gene=) from GFF3 CDS lines; not a SATURN species, it
                                 only roots the tree so OrthoFinder writes root-level HOGs
The bridged-fraction check counts dropped ALT proteins as intentionally excluded.
All genes in the proteome are kept (not only HV genes): orthogroups are better inferred from whole proteomes;
the benchmarks restrict to the SATURN genes afterwards.

Sequence clean-up: trailing '*' removed; internal '*' -> 'X' (OrthoFinder/DIAMOND do not accept stops).
FASTA headers are short neutral IDs (Mlei_000001 ...); id_map.tsv maps them back.

Species sets (--species; default = the original 4, so the 2026-10-15 run is reproducible):
    Mlei,Cgig,Crob,Drer                      tree (Mlei,(Cgig,(Crob,Drer)));
    Mlei,Cgig,Crob,Drer,Dmel,Hsap,Sros       tree (Sros,(Mlei,((Cgig,Dmel),(Crob,(Drer,Hsap)))));
        With the outgroup, N0 = Choanozoa (root) and the Metazoa-level HOGs are the node joining Mlei and the
        bilaterians (check the label in SpeciesTree_rooted_node_labels.txt; expected N1).
        Ctenophore-first is assumed; with one ctenophore and no sponge, sponge-first gives the same Metazoa node.

Outputs (--out_dir, default 06_orthofinder/):
  input/<sp>.faa            one longest protein per gene
  id_map.tsv                species, of_id, gene_id, protein_id, length, n_proteins_for_gene
  species_tree.nwk          fixed tree for OrthoFinder -s (see above)
  prep_summary.txt          proteins read, bridged, genes written, per species
"""
import argparse
import re
from collections import defaultdict
from pathlib import Path

BASE = Path("/scratch/dark_genes/SATURN_Mnemi")
ap = argparse.ArgumentParser()
ap.add_argument("--base", default=str(BASE))
ap.add_argument("--out_dir", default="")
ap.add_argument("--species", default="Mlei,Cgig,Crob,Drer",
                help="comma-separated; 6-species run: Mlei,Cgig,Crob,Drer,Dmel,Hsap,Sros")
ap.add_argument("--min_bridged_frac", type=float, default=0.8,
                help="stop with an error if fewer proteins than this are bridged to a gene")
args = ap.parse_args()
BASE = Path(args.base)
P = BASE / "01_proteomes"
OUT = Path(args.out_dir) if args.out_dir else BASE / "06_orthofinder"
(OUT / "input").mkdir(parents=True, exist_ok=True)

ALL_CFG = {
    "Mlei": {"fasta": P / "Mlei/Mlei.protein.faa", "gtf": None, "strategy": "direct"},
    "Cgig": {"fasta": P / "Cgig/Cgig.protein.faa", "gtf": P / "Cgig/Cgig.genomic.gtf", "strategy": "gtf_loc_to_xp"},
    "Crob": {"fasta": P / "Crob/Crob_KY21_ghost.faa", "gtf": None, "strategy": "direct_strip_suffix"},
    "Drer": {"fasta": P / "Drer/Drer.protein.faa", "gtf": P / "Drer/Drer.genomic.gtf", "strategy": "gtf_symbol_to_ensdarg"},
    "Dmel": {"fasta": P / "Dmel/Drosophila_melanogaster.BDGP6.46.pep.all.fa", "gtf": None, "strategy": "ensembl_gene_id"},
    "Hsap": {"fasta": P / "Hsap/Homo_sapiens.GRCh38.pep.all.fa", "gtf": None, "strategy": "ensembl_gene_id"},
    "Sros": {"fasta": P / "Sros/Sros.protein.faa", "gtf": P / "Sros/Sros.genomic.gff", "strategy": "gff_cds_locus"},
}
TREES = {
    ("Cgig", "Crob", "Drer", "Mlei"): "(Mlei,(Cgig,(Crob,Drer)));",
    ("Cgig", "Crob", "Dmel", "Drer", "Hsap", "Mlei", "Sros"): "(Sros,(Mlei,((Cgig,Dmel),(Crob,(Drer,Hsap)))));",
}
SEL = [x.strip() for x in args.species.split(",") if x.strip()]
unknown = [x for x in SEL if x not in ALL_CFG]
if unknown:
    raise SystemExit(f"Unknown species {unknown}; known: {list(ALL_CFG)}")
if tuple(sorted(SEL)) not in TREES:
    raise SystemExit(f"No species tree defined for {sorted(SEL)}; add one to TREES")
CFG = {sp: ALL_CFG[sp] for sp in SEL}
for sp, c in CFG.items():
    for k in ("fasta", "gtf"):
        if c[k] is not None and not Path(c[k]).exists():
            raise SystemExit(f"{sp}: {k} not found: {c[k]}")
LOG = []
ALT_PIDS = set()  # Drer proteins on ALT_CTG contigs (a primary-assembly protein is preferred when a gene has one)
EXCLUDED = {}   # species -> proteins intentionally left out (Drer unnamed ALT_CTG duplicates; Hsap non-primary scaffolds)


def log(m=""):
    print(m, flush=True); LOG.append(str(m))


def read_fasta(path):
    """returns {id: (header_line, seq)}"""
    out, name, head, buf = {}, None, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out[name] = (head, "".join(buf))
                head = line[1:]
                name, buf = head.split()[0], []
            elif line:
                buf.append(line.strip())
    if name is not None:
        out[name] = (head, "".join(buf))
    return out


def gtf_attr_pairs(gtf, feature, key_a, key_b):
    pat_a, pat_b = re.compile(key_a + r' "([^"]+)"'), re.compile(key_b + r' "([^"]+)"')
    out = {}
    with open(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != feature:
                continue
            a, b = pat_a.search(p[8]), pat_b.search(p[8])
            if a and b:
                out[a.group(1)] = b.group(1)
    return out


def bridge(sp, recs):
    s = CFG[sp]["strategy"]
    if s == "direct":
        pat = re.compile(r"^g(\d+)_")
        out = {}
        for pid in recs:
            m = pat.match(pid)
            out[pid] = f"g-{m.group(1)}" if m else pid
        return out
    if s == "direct_strip_suffix":
        pat = re.compile(r"^(KY21[._][^._]+[._]\d+)")
        return {pid: m.group(1).replace("_", ".") for pid in recs if (m := pat.match(pid))}
    if s == "gtf_loc_to_xp":
        xp2loc = gtf_attr_pairs(CFG[sp]["gtf"], "CDS", "protein_id", "gene_id")
        return {pid: xp2loc[pid] for pid in recs if pid in xp2loc}
    if s == "gtf_symbol_to_ensdarg":
        g2sym = gtf_attr_pairs(CFG[sp]["gtf"], "gene", "gene_id", "gene_name")
        gre = re.compile(r"gene:(ENSDARG\d+)")
        out = {}
        alt = re.compile(r"[:\s](CHR_)?ALT_CTG")
        n_fallback = n_alt = n_alt_named = n_alt_dropped = 0
        for pid, (head, _) in recs.items():
            m = gre.search(head)
            is_alt = bool(alt.search(head))
            n_alt += is_alt
            if is_alt:
                ALT_PIDS.add(pid)
            if m and m.group(1) in g2sym:
                # named gene: ALT-haplotype copies share the symbol, so they collapse into the same gene
                out[pid] = g2sym[m.group(1)]
                n_alt_named += is_alt
            elif m and is_alt:
                # unnamed gene on an alternate-haplotype contig (GRCz11 ALT_CTG): a duplicate of a primary-assembly
                # locus under a separate ENSDARG id -> dropped, as is standard for orthology inference
                n_alt_dropped += 1
            elif m:
                # unnamed gene on the primary assembly: not among SATURN's genes (the h5ad uses symbols), but kept
                # under its ENSDARG id so OrthoFinder sees the whole proteome
                out[pid] = m.group(1)
                n_fallback += 1
        log(f"  Drer: {n_alt:,} proteins on ALT_CTG contigs: {n_alt_named:,} from named genes (collapse with their "
            f"primary copy by symbol), {n_alt_dropped:,} from unnamed genes (dropped as alternate-haplotype duplicates)")
        EXCLUDED[sp] = n_alt_dropped
        log(f"  Drer: {n_fallback:,} proteins from unnamed primary-assembly genes kept under their ENSDARG id")
        return out
    if s == "ensembl_gene_id":
        gre = re.compile(r"gene:([^ .]+)")          # same regex as 02_embed.sh step1_ensembl_gene_id.py
        loc = re.compile(r"^\S+\s+pep\s+(\S+)")   # e.g. chromosome:GRCh38:14:..., scaffold:GRCh38:HSCHR7_2_CTG6:...
        out, n_nonprim, nonprim_genes = {}, 0, set()
        for pid, (head, _) in recs.items():
            m = gre.search(head)
            if not m:
                continue
            lm = loc.search(head)
            parts = lm.group(1).split(":") if lm else []
            if len(parts) > 2 and parts[0] == "scaffold" and (parts[2].startswith("HSCHR") or parts[2].startswith("HG")):
                n_nonprim += 1; nonprim_genes.add(m.group(1))
                continue
            out[pid] = m.group(1)
        EXCLUDED[sp] = n_nonprim
        log(f"  {sp}: {n_nonprim:,} proteins ({len(nonprim_genes):,} genes) on alternate-haplotype / patch scaffolds "
            f"(HSCHR*, HG*) dropped")
        return out
    if s == "gff_cds_locus":
        pid_re = re.compile(r"protein_id=([^;]+)")
        key_res = [re.compile(r"locus_tag=([^;]+)"), re.compile(r"gene=([^;]+)")]
        p2g = {}
        with open(CFG[sp]["gtf"]) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[2] != "CDS":
                    continue
                m = pid_re.search(p[8])
                if not m:
                    continue
                for kr in key_res:
                    k = kr.search(p[8])
                    if k:
                        p2g[m.group(1)] = k.group(1); break
        out = {pid: p2g.get(pid, pid) for pid in recs}
        n_self = sum(1 for pid in recs if pid not in p2g)
        log(f"  {sp}: {len(recs) - n_self:,} proteins bridged to a locus via GFF CDS lines; "
            f"{n_self:,} kept as their own gene (no CDS record)")
        return out
    raise ValueError(s)


def clean(seq):
    seq = seq.rstrip("*")
    return seq.replace("*", "X")


rows = ["species\tof_id\tgene_id\tprotein_id\tlength\tn_proteins_for_gene"]
for sp, cfg in CFG.items():
    recs = read_fasta(cfg["fasta"])
    b = bridge(sp, recs)
    frac = len(b) / max(len(recs) - EXCLUDED.get(sp, 0), 1)
    by_gene = defaultdict(list)
    for pid, gene in b.items():
        by_gene[gene].append(pid)
    n_stop = sum(1 for pid in b if "*" in recs[pid][1].rstrip("*"))
    log(f"{sp}: {len(recs):,} proteins in FASTA, {len(b):,} bridged to a gene ({100 * frac:.1f}%), "
        f"{len(by_gene):,} genes; {EXCLUDED.get(sp, 0):,} excluded on purpose; {n_stop} proteins with internal stops (-> X)")
    if frac < args.min_bridged_frac:
        raise SystemExit(f"{sp}: only {100 * frac:.1f}% of proteins bridged -- check the FASTA headers / strategy")
    with open(OUT / "input" / f"{sp}.faa", "w") as fo:
        for i, gene in enumerate(sorted(by_gene), 1):
            pids = by_gene[gene]
            best = max(pids, key=lambda p: (p not in ALT_PIDS, len(clean(recs[p][1])), p))  # primary assembly first, then longest
            seq = clean(recs[best][1])
            of_id = f"{sp}_{i:06d}"
            fo.write(f">{of_id}\n")
            for j in range(0, len(seq), 80):
                fo.write(seq[j:j + 80] + "\n")
            rows.append(f"{sp}\t{of_id}\t{gene}\t{best}\t{len(seq)}\t{len(pids)}")
(OUT / "id_map.tsv").write_text("\n".join(rows) + "\n")
(OUT / "species_tree.nwk").write_text(TREES[tuple(sorted(SEL))] + "\n")
(OUT / "prep_summary.txt").write_text("\n".join(LOG) + "\n")
log("Wrote " + str(OUT / "input") + ", id_map.tsv, species_tree.nwk")
