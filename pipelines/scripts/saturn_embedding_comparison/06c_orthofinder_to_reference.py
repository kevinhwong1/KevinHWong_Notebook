#!/usr/bin/env python3
"""
06c_orthofinder_to_reference.py  (runs on Pegasus after 06b; a minute)

Turns OrthoFinder output into gene-level reference tables in the same shape the benchmark scripts use
(species, gene_id, protein_id, group), so 05n / 05p can be rerun with --reference of_og or of_hog.

  ALL4_gene_orthofinder_OG.csv     Orthogroups.tsv (+ Orthogroups_UnassignedGenes.tsv as singletons)
  ALL4_gene_orthofinder_N0HOG.csv  Phylogenetic_Hierarchical_Orthogroups/N0.tsv (root = Metazoa-level HOGs;
                                   genes absent from N0.tsv become singletons)
  Columns: species, gene_id, protein_id, group, singleton (True if the gene's group has only that gene)
  coverage_vs_saturn.csv           per species: SATURN genes (ESM-1b macrogene-only run) found in the tables,
                                   and in groups with >= 2 species  (sanity check of the gene bridges)
"""
import argparse
import pickle
from pathlib import Path
import pandas as pd

BASE = Path("/scratch/dark_genes/SATURN_Mnemi")
ap = argparse.ArgumentParser()
ap.add_argument("--of_dir", default=str(BASE / "06_orthofinder"))
ap.add_argument("--results", default="", help="OrthoFinder results folder (default: newest under 06_orthofinder/Results_20261009)")
ap.add_argument("--saturn_run", default=str(BASE / "04_saturn_runs/20261002_Mlei_Cgig_Crob_Drer_K3000_ESM1b_macrogene"))
args = ap.parse_args()
OF = Path(args.of_dir)
SPECIES = ["Mlei", "Cgig", "Crob", "Drer"]


def find_one(root, name):
    hits = sorted(Path(root).rglob(name), key=lambda p: p.stat().st_mtime)
    if not hits:
        raise SystemExit(f"{name} not found under {root}")
    return hits[-1]


root = Path(args.results) if args.results else OF / "Results_20261009"
og_tsv = find_one(root, "Orthogroups.tsv")
res = og_tsv.parent.parent  # .../Results_<date>/
log = [f"OrthoFinder results: {res}"]
idmap = pd.read_csv(OF / "id_map.tsv", sep="\t")


def melt(tsv, id_col, extra_cols=()):
    d = pd.read_csv(tsv, sep="\t", dtype=str)
    spcols = [c for c in d.columns if c in SPECIES]
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
OG.to_csv(OF / "ALL4_gene_orthofinder_OG.csv", index=False)
n0 = find_one(res, "N0.tsv")
HOG = finish(melt(n0, "HOG"), "HOG")
HOG.to_csv(OF / "ALL4_gene_orthofinder_N0HOG.csv", index=False)

# coverage of SATURN genes
pk = sorted((Path(args.saturn_run) / "saturn_results").glob("*genes_to_macrogenes*.pkl"))
keys = list(pickle.load(open(pk[0], "rb")).keys())
sat = pd.DataFrame({"species": [k.split("_", 1)[0] for k in keys], "gene_id": [k.split("_", 1)[1] for k in keys]})
cov = []
for lab, T in [("OG", OG), ("N0HOG", HOG)]:
    nsp = T.groupby("group")["species"].nunique()
    T = T.assign(multi=T["group"].map(nsp) >= 2)
    t_low = T.assign(gl=T["gene_id"].str.lower())
    for s in SPECIES:
        g = sat[sat.species == s]
        ts = t_low[t_low.species == s]
        exact = g["gene_id"].isin(ts["gene_id"]).mean()
        ci = g["gene_id"].str.lower().isin(ts["gl"]).mean()
        multi = g["gene_id"].str.lower().isin(ts.loc[ts["multi"], "gl"]).mean()
        cov.append({"table": lab, "species": s, "saturn_genes": len(g), "pct_found_exact": 100 * exact,
                    "pct_found_case_insensitive": 100 * ci, "pct_in_multispecies_group": 100 * multi})
    log.append(f"{lab}: {T['group'].nunique():,} groups, {int((~T['singleton']).sum()):,} genes in groups of >= 2, "
               f"{int((T['group'].map(nsp) == 4).sum()):,} genes in 4-species groups")
C = pd.DataFrame(cov)
C.to_csv(OF / "coverage_vs_saturn.csv", index=False)
log.append(C.round(1).to_string(index=False))
(OF / "06c_summary.txt").write_text("\n".join(log) + "\n")
print("\n".join(log))
