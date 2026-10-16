#!/bin/bash
#BSUB -J qc_umap_Dmel_Hsap
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 4
#BSUB -W 4:00
#BSUB -R "rusage[mem=64000] span[hosts=1]"
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/qc_umap_Dmel_Hsap.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/qc_umap_Dmel_Hsap.%J.err
#BSUB -u kxw755@earth.miami.edu
#BSUB -N
# =============================================================================
# 07f_QC_UMAP_Dmel_Hsap.sh  --  DIAGNOSTICS ONLY. Reads the 07d SATURN inputs, writes
# QC tables, per-species UMAPs and marker dot plots. Does NOT change or filter the
# *_saturn.h5ad files (the ESM embedding job and SATURN keep reading those as they are).
#
# Inputs : 02_processed_h5ad/{Dmel,Hsap}_saturn.h5ad      (from 07d; raw counts in .X)
#          02_processed_h5ad/{Mlei,Cgig,Crob,Drer}_saturn.h5ad   (cross-species QC comparison only)
# Outputs: 07_new_species_QC/20261016_Dmel_Hsap_QC/
#            {SP}_qc_by_celltype.csv        n cells, median UMI / genes / %mt / %ribo, MAD-outlier counts
#            {SP}_cells_umap_qc.csv         per-cell obs + QC + UMAP coords + leiden (for R/ggplot locally)
#            {SP}_qc_umap.h5ad              log-normalised, HVG/PCA/UMAP (for interactive checks)
#            {SP}_marker_genes_found.csv    which canonical markers were present in var
#            crossspecies_qc_per_cell.csv, crossspecies_qc_summary.csv
#            figures/*.png + *.pdf
#            summary.txt
#
# Processing (per species, no integration, no filtering):
#   normalize_total 1e4 -> log1p -> HVG 2,000 (seurat_v3 on counts; falls back to 'seurat'
#   on log data if scikit-misc is missing) -> scale (max 10) -> PCA 50 -> neighbors
#   (k=15, 30 PCs) -> UMAP -> Leiden res 1.0 (reference only). Hsap: no donor correction,
#   so the donor UMAP doubles as a batch check for the 2 donors (D496, D503).
#   QC outlier FLAGS (not filters): >5 MADs on log1p total counts or log1p n_genes,
#   or %mt > median + 3 MADs (scverse best-practice thresholds).
#
# cat-then-python. NOTE: the Python below contains no pipe characters (Pegasus heredoc issue).
# Submit:  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/07f_QC_UMAP_Dmel_Hsap.sh
# If 07d still has to be rerun for the progenitor/erythroid revision, chain it:
#   bsub < 07d_prepare_Dmel_Hsap.sh
#   bsub -w "done(prep_Dmel_Hsap)" < 07f_QC_UMAP_Dmel_Hsap.sh
# =============================================================================
set -eo pipefail

BASE="/scratch/dark_genes/SATURN_Mnemi"
OUT="${BASE}/07_new_species_QC/20261016_Dmel_Hsap_QC"
mkdir -p "${OUT}/figures" "${BASE}/logs"

# Batch jobs don't read ~/.bashrc: same activation as 07c (conda, micromamba fallback)
set +eu
if [ -f ~/anaconda3/bin/activate ] && ~/anaconda3/bin/conda env list | grep -q "saturn_protein"; then
    source ~/anaconda3/bin/activate
    conda activate saturn_protein
else
    export MAMBA_ROOT_PREFIX=/nethome/kxw755/micromamba
    eval "$(/nethome/kxw755/.local/bin/micromamba shell hook --shell bash)"
    micromamba activate saturn_protein
fi
set -eo pipefail
echo "python: $(which python)"
python -c "import scanpy, anndata; print('scanpy', scanpy.__version__, 'anndata', anndata.__version__)"
echo "=== 07f_QC_UMAP_Dmel_Hsap.sh ==="; date; hostname

cat > "${OUT}/07f_QC_UMAP_Dmel_Hsap.py" << 'EOF'
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

SEED = 42
BASE = Path(os.environ.get("QC_BASE", "/scratch/dark_genes/SATURN_Mnemi"))
H5 = BASE / "02_processed_h5ad"
OUT = BASE / "07_new_species_QC" / "20261016_Dmel_Hsap_QC"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sc.settings.verbosity = 1
sc.settings.n_jobs = 4
plt.rcParams["savefig.dpi"] = 200
plt.rcParams["savefig.bbox"] = "tight"

LOG = []
def say(s=""):
    print(s, flush=True)
    LOG.append(str(s))

def savefig(fig, name):
    fig.savefig(FIG / (name + ".png"))
    fig.savefig(FIG / (name + ".pdf"))
    plt.close(fig)

# ---------------------------------------------------------------------------
# Species settings
# ---------------------------------------------------------------------------
# Marker panels (by gene symbol; any symbol not in var is dropped and reported)
HSAP_MARKERS = {
    "Classical mono":  ["CD14", "LYZ", "S100A8", "VCAN"],
    "Nonclass. mono":  ["FCGR3A", "CDKN1C", "LST1"],
    "Macrophage":      ["C1QA", "C1QB", "CD68", "APOE", "MARCO", "FABP4", "CD163", "HMOX1"],
    "DC":              ["CLEC9A", "XCR1", "CD1C", "FCER1A", "CLEC10A", "CCR7", "LAMP3"],
    "Mast":            ["KIT", "TPSAB1", "CPA3"],
    "Progenitor":      ["CD34", "PRSS57", "SPINK2"],
    "Erythroid":       ["HBB", "GYPA", "ALAS2"],
    "T":               ["CD3E", "CD4", "CD8A", "TRDC", "KLRB1"],
    "NK":              ["NKG7", "GNLY", "NCAM1"],
    "B / plasma":      ["MS4A1", "CD79A", "JCHAIN", "MZB1"],
    "Cycling":         ["MKI67", "TOP2A"],
}
DMEL_MARKERS = {
    "Pan-hemocyte":    ["srp", "Hml", "He"],
    "Plasmatocyte":    ["NimC1", "eater", "Pxn", "crq", "Col4a1"],
    "PL subsets":      ["Rel", "Inos", "vir-1", "ImpL2", "Lsp2", "ppl", "robo2", "Pcd"],
    "AMP":             ["Drs", "Mtk", "AttA", "DptA", "IM1"],
    "Crystal cell":    ["lz", "PPO1", "PPO2", "peb"],
    "Lamellocyte":     ["atilla", "PPO3", "mys", "alphaPS4", "Itgbn"],
    "Proliferating":   ["stg", "PCNA", "polo", "Mcm2"],
}
HSAP_ORDER = ["Progenitor", "Cycling", "Classical monocytes", "Nonclassical monocytes",
              "Alveolar macrophages", "Intermediate macrophages", "Intestinal macrophages",
              "Erythrophagocytic macrophages", "DC1", "DC2", "migDC", "Mast cells", "Erythroid"]

SPECIES = {
    "Dmel": dict(mt_prefix=("mt:",), ribo_prefix=("RpS", "RpL"), markers=DMEL_MARKERS,
                 extra_cols=[], order=None),
    "Hsap": dict(mt_prefix=("MT-",), ribo_prefix=("RPS", "RPL"), markers=HSAP_MARKERS,
                 extra_cols=["tissue", "donor_id", "sample_group"], order=HSAP_ORDER),
}

def mad(x):
    return np.median(np.abs(x - np.median(x)))

def outlier(x, nmads):
    m = mad(x)
    return np.abs(x - np.median(x)) > nmads * m if m > 0 else np.zeros(len(x), dtype=bool)

def basic_qc(a, mt_prefix, ribo_prefix):
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    sym = a.var["gene_symbol"].astype(str) if "gene_symbol" in a.var else pd.Series(a.var_names, index=a.var_names)
    tot = np.asarray(X.sum(1)).ravel()
    ngen = np.asarray((X > 0).sum(1)).ravel()
    mt = sym.str.startswith(mt_prefix).values
    ribo = sym.str.startswith(ribo_prefix).values
    den = np.maximum(tot, 1)
    a.obs["total_counts"] = tot
    a.obs["n_genes"] = ngen
    a.obs["pct_mt"] = 100 * np.asarray(X[:, mt].sum(1)).ravel() / den
    a.obs["pct_ribo"] = 100 * np.asarray(X[:, ribo].sum(1)).ravel() / den
    return int(mt.sum()), int(ribo.sum())

def order_celltypes(a, order):
    cts = a.obs["cell_type"].astype(str)
    counts = cts.value_counts()
    if order:
        lymph = sorted([c for c in counts.index if c not in order])
        cats = [c for c in order if c in counts.index] + lymph
    else:
        cats = list(counts.index)
    a.obs["cell_type"] = pd.Categorical(cts, categories=cats)

def pick_markers(a, panel):
    sym_set = set(a.var_names)
    found, rows = {}, []
    for grp, genes in panel.items():
        g = [x for x in genes if x in sym_set]
        for x in genes:
            rows.append(dict(group=grp, gene=x, present=x in sym_set))
        if g:
            found[grp] = g
    return found, pd.DataFrame(rows)

# ---------------------------------------------------------------------------
# Per-species QC + UMAP
# ---------------------------------------------------------------------------
summary_rows = []
for SP, cfg in SPECIES.items():
    say("\n" + "=" * 80 + "\n" + SP + "\n" + "=" * 80)
    path = H5 / (SP + "_saturn.h5ad")
    a = sc.read_h5ad(path)
    say("loaded %s: %d cells x %d genes; %d cell_type labels" % (path.name, a.n_obs, a.n_vars, a.obs["cell_type"].nunique()))
    say("file modified: %s" % pd.Timestamp(path.stat().st_mtime, unit="s"))
    if SP == "Hsap":
        have = set(a.obs["cell_type"].astype(str))
        missing = [x for x in ["Progenitor", "Erythroid"] if x not in have]
        if missing:
            say("WARNING: %s not in Hsap cell_type -> this is the PRE-revision 07d output; rerun 07d to include them" % missing)
        else:
            say("Hsap includes Progenitor + Erythroid (post-revision 07d output)")

    n_mt, n_ribo = basic_qc(a, cfg["mt_prefix"], cfg["ribo_prefix"])
    say("mito genes: %d   ribosomal genes: %d" % (n_mt, n_ribo))
    lt = np.log1p(a.obs["total_counts"].values)
    lg = np.log1p(a.obs["n_genes"].values)
    pm = a.obs["pct_mt"].values
    f_tot = outlier(lt, 5)
    f_gen = outlier(lg, 5)
    m = mad(pm)
    f_mt = pm > (np.median(pm) + 3 * m) if m > 0 else np.zeros(len(pm), dtype=bool)
    a.obs["flag_counts"] = f_tot
    a.obs["flag_genes"] = f_gen
    a.obs["flag_mt"] = f_mt
    a.obs["qc_outlier"] = np.logical_or(np.logical_or(f_tot, f_gen), f_mt)
    say("MAD outlier flags (not removed): counts %d, genes %d, mt %d, any %d (%.1f%%)"
        % (f_tot.sum(), f_gen.sum(), f_mt.sum(), a.obs["qc_outlier"].sum(), 100 * a.obs["qc_outlier"].mean()))
    say("mt cutoff used for flag: %.2f%%" % (np.median(pm) + 3 * m))
    say(a.obs[["total_counts", "n_genes", "pct_mt", "pct_ribo"]].describe(percentiles=[0.05, 0.5, 0.95]).round(1).to_string())

    order_celltypes(a, cfg["order"])
    g = a.obs.groupby("cell_type", observed=True)
    tab = pd.DataFrame({
        "n_cells": g.size(),
        "median_total_counts": g["total_counts"].median().round(0),
        "median_n_genes": g["n_genes"].median().round(0),
        "median_pct_mt": g["pct_mt"].median().round(2),
        "median_pct_ribo": g["pct_ribo"].median().round(2),
        "n_flag_counts": g["flag_counts"].sum(),
        "n_flag_genes": g["flag_genes"].sum(),
        "n_flag_mt": g["flag_mt"].sum(),
        "n_qc_outlier": g["qc_outlier"].sum(),
    })
    tab["pct_qc_outlier"] = (100 * tab["n_qc_outlier"] / tab["n_cells"]).round(1)
    tab.to_csv(OUT / (SP + "_qc_by_celltype.csv"))
    say("\nQC by cell type:\n" + tab.to_string())
    if SP == "Hsap":
        say("\ncell_type x tissue:\n" + pd.crosstab(a.obs["cell_type"], a.obs["tissue"]).to_string())
        say("\ncell_type x donor:\n" + pd.crosstab(a.obs["cell_type"], a.obs["donor_id"]).to_string())

    # ---- processing (symbols as var_names for plotting) ----
    a.var["gene_id"] = a.var_names.astype(str)
    a.var_names = pd.Index(a.var["gene_symbol"].astype(str).values)
    a.var_names_make_unique()
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    try:
        sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat_v3", layer="counts")
        hvg_flavor = "seurat_v3"
    except Exception as e:
        say("seurat_v3 HVG failed (%s); using flavor='seurat' on log data" % type(e).__name__)
        sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat")
        hvg_flavor = "seurat"
    say("HVG flavor: %s (%d genes)" % (hvg_flavor, a.var["highly_variable"].sum()))
    a.raw = a
    b = a[:, a.var["highly_variable"]].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=50, random_state=SEED)
    a.obsm["X_pca"] = b.obsm["X_pca"]
    a.uns["pca"] = b.uns["pca"]
    del b
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)
    try:
        sc.tl.leiden(a, resolution=1.0, random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
    except TypeError:
        sc.tl.leiden(a, resolution=1.0, random_state=SEED)
    say("Leiden res1.0 clusters: %d" % a.obs["leiden"].nunique())
    say("\nleiden x cell_type (rows = leiden):\n" + pd.crosstab(a.obs["leiden"], a.obs["cell_type"]).to_string())

    # ---- per-cell export for local plotting ----
    cols = ["cell_type"] + cfg["extra_cols"] + ["total_counts", "n_genes", "pct_mt", "pct_ribo",
            "flag_counts", "flag_genes", "flag_mt", "qc_outlier", "leiden"]
    cells = a.obs[cols].copy()
    cells["UMAP1"] = a.obsm["X_umap"][:, 0]
    cells["UMAP2"] = a.obsm["X_umap"][:, 1]
    cells.to_csv(OUT / (SP + "_cells_umap_qc.csv"))

    # ---- figures ----
    big = dict(frameon=False, show=False, size=max(4, 120000 / a.n_obs) if a.n_obs > 0 else 10)
    big["size"] = min(big["size"], 25)
    fig, ax = plt.subplots(figsize=(9, 7))
    sc.pl.umap(a, color="cell_type", ax=ax, title=SP + " - author cell type (n=%d)" % a.n_obs, **big)
    savefig(fig, SP + "_umap_celltype")
    fig, ax = plt.subplots(figsize=(9, 8))
    sc.pl.umap(a, color="cell_type", ax=ax, legend_loc="on data", legend_fontsize=7,
               legend_fontoutline=2, title=SP + " - author cell type", **big)
    savefig(fig, SP + "_umap_celltype_ondata")

    qc_cols = ["total_counts", "n_genes", "pct_mt", "pct_ribo"]
    a.obs["log10_total_counts"] = np.log10(a.obs["total_counts"] + 1)
    fig, axs = plt.subplots(1, 4, figsize=(22, 5))
    for ax, c in zip(axs, ["log10_total_counts", "n_genes", "pct_mt", "pct_ribo"]):
        sc.pl.umap(a, color=c, ax=ax, cmap="viridis", **big)
    savefig(fig, SP + "_umap_qc_metrics")

    a.obs["qc_outlier_str"] = pd.Categorical(np.where(a.obs["qc_outlier"], "flagged", "ok"), categories=["ok", "flagged"])
    a.uns["qc_outlier_str_colors"] = ["#d0d0d0", "#d62728"]
    meta = ["qc_outlier_str", "leiden"] + cfg["extra_cols"]
    fig, axs = plt.subplots(1, len(meta), figsize=(9 * len(meta), 6))
    fig.subplots_adjust(wspace=0.45)
    for ax, c in zip(np.atleast_1d(axs), meta):
        sc.pl.umap(a, color=c, ax=ax, title=c, **big)
    savefig(fig, SP + "_umap_meta")

    cats = list(a.obs["cell_type"].cat.categories)
    cats = [c for c in cats if (a.obs["cell_type"] == c).sum() > 0]
    fig, axs = plt.subplots(len(qc_cols), 1, figsize=(max(8, 0.55 * len(cats) + 3), 3.2 * len(qc_cols)), sharex=True)
    for ax, c in zip(axs, qc_cols):
        vals = [a.obs.loc[a.obs["cell_type"] == k, c].values.astype(float) for k in cats]
        if c == "total_counts":
            vals = [np.log10(v + 1) for v in vals]
        parts = ax.violinplot(vals, positions=np.arange(len(cats)), showmedians=True, widths=0.8)
        for body in parts["bodies"]:
            body.set_facecolor("#4c72b0"); body.set_alpha(0.6)
        ax.set_ylabel("log10(total_counts + 1)" if c == "total_counts" else c)
        ax.grid(axis="y", alpha=0.3)
    axs[-1].set_xticks(np.arange(len(cats)))
    axs[-1].set_xticklabels(cats, rotation=90)
    axs[0].set_title(SP + " QC by author cell type (medians marked)")
    savefig(fig, SP + "_violin_qc_by_celltype")

    fig, ax = plt.subplots(figsize=(max(7, 0.45 * len(tab) + 3), 4))
    tab["n_cells"].plot.bar(ax=ax, color="#4c72b0")
    ax.set_ylabel("cells"); ax.set_xlabel(""); ax.set_title(SP + " cells per author cell type")
    for i, v in enumerate(tab["n_cells"].values):
        ax.text(i, v, str(int(v)), ha="center", va="bottom", fontsize=7)
    savefig(fig, SP + "_bar_celltype_counts")
    if SP == "Hsap":
        ct = pd.crosstab(a.obs["cell_type"], a.obs["tissue"])
        fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(ct) + 4), 5))
        ct.plot.bar(stacked=True, ax=ax, colormap="tab20", width=0.85)
        ax.set_ylabel("cells"); ax.set_xlabel(""); ax.legend(title="tissue", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7)
        ax.set_title("Hsap cell type x tissue")
        savefig(fig, "Hsap_bar_celltype_by_tissue")

    found, mrk = pick_markers(a, cfg["markers"])
    mrk.to_csv(OUT / (SP + "_marker_genes_found.csv"), index=False)
    miss = mrk.loc[~mrk["present"], "gene"].tolist()
    say("markers missing from var: %s" % (miss if miss else "none"))
    if found:
        dp = sc.pl.dotplot(a, found, groupby="cell_type", standard_scale="var", show=False, return_fig=True)
        dp.savefig(FIG / (SP + "_dotplot_markers_by_celltype.png"))
        dp.savefig(FIG / (SP + "_dotplot_markers_by_celltype.pdf"))
        plt.close("all")
        allg = [x for v in found.values() for x in v]
        nc = 6
        nr = int(np.ceil(len(allg) / nc))
        fig, axs = plt.subplots(nr, nc, figsize=(3.6 * nc, 3.3 * nr))
        for ax, gname in zip(axs.ravel(), allg):
            sc.pl.umap(a, color=gname, ax=ax, cmap="Reds", use_raw=True, frameon=False, show=False,
                       size=big["size"], colorbar_loc=None)
        for ax in axs.ravel()[len(allg):]:
            ax.axis("off")
        savefig(fig, SP + "_umap_markers")

    a.obs = a.obs.drop(columns=["qc_outlier_str", "log10_total_counts"])
    a.write_h5ad(OUT / (SP + "_qc_umap.h5ad"), compression="gzip")
    summary_rows.append(dict(species=SP, n_cells=a.n_obs, n_genes=a.n_vars,
                             n_cell_types=a.obs["cell_type"].nunique(),
                             median_total_counts=float(np.median(a.obs["total_counts"])),
                             median_n_genes=float(np.median(a.obs["n_genes"])),
                             median_pct_mt=float(np.median(a.obs["pct_mt"])),
                             n_qc_outlier=int(a.obs["qc_outlier"].sum()),
                             n_leiden=a.obs["leiden"].nunique(), hvg_flavor=hvg_flavor))
    del a

# ---------------------------------------------------------------------------
# Cross-species QC comparison (raw counts as SATURN sees them)
# ---------------------------------------------------------------------------
say("\n" + "=" * 80 + "\nCross-species QC (SATURN inputs, raw counts)\n" + "=" * 80)
per_cell = []
for SP in ["Mlei", "Cgig", "Crob", "Drer", "Dmel", "Hsap"]:
    p = H5 / (SP + "_saturn.h5ad")
    if not p.exists():
        say("  %s: %s not found, skipped" % (SP, p.name))
        continue
    a = sc.read_h5ad(p)
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    per_cell.append(pd.DataFrame({
        "species": SP,
        "cell_type": a.obs["cell_type"].astype(str).values if "cell_type" in a.obs else "NA",
        "total_counts": np.asarray(X.sum(1)).ravel(),
        "n_genes": np.asarray((X > 0).sum(1)).ravel(),
        "n_genes_in_matrix": a.n_vars,
    }, index=a.obs_names))
    del a, X
pc = pd.concat(per_cell)
pc.to_csv(OUT / "crossspecies_qc_per_cell.csv")
cs = pc.groupby("species", sort=False).agg(
    n_cells=("total_counts", "size"),
    n_genes_in_matrix=("n_genes_in_matrix", "first"),
    n_cell_types=("cell_type", "nunique"),
    median_total_counts=("total_counts", "median"),
    q05_total_counts=("total_counts", lambda x: np.percentile(x, 5)),
    q95_total_counts=("total_counts", lambda x: np.percentile(x, 95)),
    median_n_genes=("n_genes", "median"),
    q05_n_genes=("n_genes", lambda x: np.percentile(x, 5)),
    q95_n_genes=("n_genes", lambda x: np.percentile(x, 95)),
).round(0)
cs.to_csv(OUT / "crossspecies_qc_summary.csv")
say(cs.to_string())

sps = list(cs.index)
fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
for ax, c, lab in zip(axs, ["total_counts", "n_genes"], ["UMIs per cell", "genes detected per cell"]):
    data = [np.log10(pc.loc[pc["species"] == s, c].values + 1) for s in sps]
    parts = ax.violinplot(data, showmedians=True, widths=0.8)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor("#dd8452" if sps[i] in ("Dmel", "Hsap") else "#4c72b0")
        body.set_alpha(0.7)
    ax.set_xticks(range(1, len(sps) + 1)); ax.set_xticklabels(sps)
    ax.set_ylabel("log10(" + lab + " + 1)"); ax.set_title(lab + " (new species in orange)")
savefig(fig, "crossspecies_qc_violin")

pd.DataFrame(summary_rows).to_csv(OUT / "new_species_summary.csv", index=False)
(OUT / "summary.txt").write_text("\n".join(LOG) + "\n")
say("\nwrote outputs to " + str(OUT))
EOF

python "${OUT}/07f_QC_UMAP_Dmel_Hsap.py"
echo; echo "=== 07f DONE ==="; date
