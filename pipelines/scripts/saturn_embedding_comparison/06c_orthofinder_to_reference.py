#!/usr/bin/env python3
"""
06c_orthofinder_to_reference.py  (runs on Pegasus after 06b / 06b6; a minute)

Turns OrthoFinder output into gene-level reference tables in the same shape the benchmark scripts use
(species, gene_id, protein_id, group), so 05g / 05i / 05n / 05p can score against OrthoFinder.

  <PREFIX>_gene_orthofinder_OG.csv           Orthogroups.tsv (+ Orthogroups_UnassignedGenes.tsv as singletons)
  <PREFIX>_gene_orthofinder_<NODE>HOG.csv    Phylogenetic_Hierarchical_Orthogroups/<NODE>.tsv, NODE = the
                                             Metazoa-level node (genes absent from it become singletons)
  Columns: species, gene_id, protein_id, group, singleton (True if the gene's group has only that gene)
  coverage_vs_saturn.csv                     per species: SATURN genes found in the tables, and in groups with
                                             >= 2 species (sanity check of the gene bridges); needs --saturn_run
  06c_summary.txt

2026-10-19 (6-species panel; defaults reproduce the 2026-10-15 4-species run):
  --of_dir      06_orthofinder (4 sp, default) or 06_orthofinder_6sp
  species       read from <of_dir>/id_map.tsv; --outgroup species (default Sros) are dropped from the tables
                (they only root the tree); PREFIX = ALL<number of species kept>, e.g. ALL4, ALL6
  --hog_node    'auto' (default): the node of SpeciesTree_rooted_node_labels.txt whose leaves are exactly the
                kept species (= Metazoa). 4 sp without outgroup -> N0 (root; file name unchanged: ALL4_..._N0HOG.csv);
                6 sp with Sros -> expected N1. Or give a node name, e.g. N1.
  --results     OrthoFinder results folder; default: the newest Orthogroups.tsv under <of_dir>
  --saturn_run  run dir with saturn_results/*genes_to_macrogenes*.pkl for the coverage check ('' = skip)
"""
import argparse
import pickle
import re
from pathlib import Path

import pandas as pd

BASE = Path("/scratch/dark_genes/SATURN_Mnemi")
ap = argparse.ArgumentParser()
ap.add_argument("--of_dir", default=str(BASE / "06_orthofinder"))
ap.add_argument("--results", default="", help="OrthoFinder results folder (default: newest Orthogroups.tsv under --of_dir)")
ap.add_argument("--saturn_run", default=str(BASE / "04_saturn_runs/20261002_Mlei_Cgig_Crob_Drer_K3000_ESM1b_macrogene"))
ap.add_argument("--outgroup", default="Sros", help="comma-separated species to drop from the tables")
ap.add_argument("--hog_node", default="auto", help="'auto' = Metazoa node from the species tree, or e.g. N1")
args = ap.parse_args()
OF = Path(args.of_dir)
idmap = pd.read_csv(OF / "id_map.tsv", sep="\t")
OUTGROUP = [x for x in args.outgroup.split(",") if x]
ALL_SP = list(dict.fromkeys(idmap["species"]))
SPECIES = [s for s in ALL_SP if s not in OUTGROUP]
NSP = len(SPECIES)
PREFIX = f"ALL{NSP}"
idmap = idmap[idmap["species"].isin(SPECIES)]


def find_one(root, name):
    hits = sorted(Path(root).rglob(name), key=lambda p: p.stat().st_mtime)
    if not hits:
        raise SystemExit(f"{name} not found under {root}")
    return hits[-1]


root = Path(args.results) if args.results else OF
og_tsv = find_one(root, "Orthogroups.tsv")
res = og_tsv.parent.parent  # .../Results_<date>/
log = [f"OrthoFinder results: {res}", f"species kept: {SPECIES}; dropped (outgroup): {[s for s in ALL_SP if s in OUTGROUP]}"]


def leaves_by_node(newick):
    """{internal node label: set of leaf names} from a rooted, node-labelled Newick string."""
    s = re.sub(r":[0-9.eE+-]+", "", newick.strip().rstrip(";"))   # drop branch lengths
    out, stack, i, cur = {}, [], 0, set()
    while i < len(s):
        c = s[i]
        if c == "(":
            stack.append(cur); cur = set(); i += 1
        elif c in ",)":
            if c == ")":
                j = i + 1
                while j < len(s) and s[j] not in ",()":
                    j += 1
                lab = s[i + 1:j]
                if lab:
                    out[lab] = set(cur)
                parent = stack.pop(); parent |= cur; cur = parent
                i = j
            else:
                i += 1
        else:
            j = i
            while j < len(s) and s[j] not in ",()":
                j += 1
            cur.add(s[i:j].strip()); i = j
    return out


def metazoa_node():
    if args.hog_node != "auto":
        return args.hog_node
    trees = sorted(res.rglob("SpeciesTree_rooted_node_labels.txt"))
    if not trees:
        log.append("SpeciesTree_rooted_node_labels.txt not found; falling back to N0")
        return "N0"
    tree = trees[-1].read_text()
    log.append(f"Species tree with node labels: {tree.strip()}")
    nodes = leaves_by_node(tree)
    # leaf names may carry a file suffix (e.g. Mlei.faa) depending on the OrthoFinder version
    norm = {k: {re.sub(r"\.(faa|fa|fasta)$", "", x) for x in v} for k, v in nodes.items()}
    hit = [k for k, v in norm.items() if v == set(SPECIES)]
    if not hit:
        raise SystemExit(f"No node whose leaves are exactly {SPECIES}: {norm}")
    log.append(f"Metazoa-level node (leaves = kept species): {hit[0]}")
    return hit[0]


def melt(tsv, id_col):
    d = pd.read_csv(tsv, sep="\t", dtype=str)
    spcols = [c for c in d.columns if re.sub(r"\.(faa|fa|fasta)$", "", c) in SPECIES]
    rows = []
    for _, r in d.iterrows():
        for s in spcols:
            v = r[s]
            if isinstance(v, str) and v.strip():
                for x in v.split(","):
                    rows.append((r[id_col], x.strip()))
    return pd.DataFrame(rows, columns=["group", "of_id"])


def finish(df, label):
    df = idmap.merge(df, on="of_id", how="left")
    miss = df["group"].isna()
    df.loc[miss, "group"] = label + "_single_" + df.loc[miss, "of_id"]
    size = df.groupby("group")["of_id"].transform("size")
    df["singleton"] = size == 1
    return df[["species", "gene_id", "protein_id", "group", "singleton"]]


og = melt(og_tsv, "Orthogroup")
un = og_tsv.with_name("Orthogroups_UnassignedGenes.tsv")
if un.exists():
    og = pd.concat([og, melt(un, "Orthogroup")])
OG = finish(og, "OG")
OG.to_csv(OF / f"{PREFIX}_gene_orthofinder_OG.csv", index=False)
TABLES = [("OG", OG)]
log.append(f"wrote {PREFIX}_gene_orthofinder_OG.csv")

node = metazoa_node()
hog_tsv = sorted(Path(res).rglob(f"{node}.tsv"))
if hog_tsv:
    HOG = finish(melt(hog_tsv[-1], "HOG"), "HOG")
    HOG.to_csv(OF / f"{PREFIX}_gene_orthofinder_{node}HOG.csv", index=False)
    TABLES.append((f"{node}HOG", HOG))
    log.append(f"wrote {PREFIX}_gene_orthofinder_{node}HOG.csv  (from {hog_tsv[-1]})")
else:
    hog_files = sorted(p.name for p in (res / "Phylogenetic_Hierarchical_Orthogroups").glob("N*.tsv"))
    log.append(f"{node}.tsv not found (HOG files present: {', '.join(hog_files) or 'none'}): "
               f"skipping the HOG table; use the OG table only.")

# per-table group statistics (+ SATURN coverage if a run is given)
sat = None
if args.saturn_run:
    pk = sorted((Path(args.saturn_run) / "saturn_results").glob("*genes_to_macrogenes*.pkl"))
    if pk:
        keys = list(pickle.load(open(pk[0], "rb")).keys())
        sat = pd.DataFrame({"species": [k.split("_", 1)[0] for k in keys], "gene_id": [k.split("_", 1)[1] for k in keys]})
    else:
        log.append(f"no genes_to_macrogenes pkl under {args.saturn_run}; coverage check skipped")
cov = []
for lab, T in TABLES:
    nsp = T.groupby("group")["species"].nunique()
    log.append(f"{lab}: {T['group'].nunique():,} groups, {int((~T['singleton']).sum()):,} genes in groups of >= 2, "
               f"{int((T['group'].map(nsp) == NSP).sum()):,} genes in {NSP}-species groups")
    if sat is None:
        continue
    T = T.assign(multi=T["group"].map(nsp) >= 2, gl=T["gene_id"].astype(str).str.lower())
    for s in SPECIES:
        g = sat[sat.species == s]
        if g.empty:
            continue
        ts = T[T.species == s]
        cov.append({"table": lab, "species": s, "saturn_genes": len(g),
                    "pct_found_exact": 100 * g["gene_id"].isin(ts["gene_id"]).mean(),
                    "pct_found_case_insensitive": 100 * g["gene_id"].str.lower().isin(ts["gl"]).mean(),
                    "pct_in_multispecies_group": 100 * g["gene_id"].str.lower().isin(ts.loc[ts["multi"], "gl"]).mean()})
if cov:
    C = pd.DataFrame(cov)
    C.to_csv(OF / "coverage_vs_saturn.csv", index=False)
    log.append(C.round(1).to_string(index=False))
(OF / "06c_summary.txt").write_text("\n".join(log) + "\n")
print("\n".join(log))
