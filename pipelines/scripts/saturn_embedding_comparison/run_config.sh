#!/usr/bin/env bash
# ============================================================
# run_config.sh  —  MASTER CONFIGURATION for SATURN_Mnemi
#
# 5-species immune integration:
#   Mlei  Mnemiopsis leidyi      (ctenophore)
#   Cgig  Crassostrea gigas      (oyster)
#   Spur  Strongylocentrotus purpuratus (sea urchin)
#   Crob  Ciona robusta          (tunicate)
#   Drer  Danio rerio            (zebrafish)
#
# SOURCE at the top of every pipeline script:
#   source "$(dirname "$0")/run_config.sh"
# ============================================================

# ------------------------------------------------------------
# TOP-LEVEL PATHS
# ------------------------------------------------------------
export BASE="/scratch/dark_genes/SATURN_Mnemi"
export SCRIPTS_DIR="${BASE}/scripts"
export SATURN_DIR="${HOME}/SATURN"

# ------------------------------------------------------------
# SATURN SCRIPT PATHS
# ------------------------------------------------------------
export SATURN_TRAIN="${SATURN_DIR}/train-saturn.py"
export ESM_EXTRACT="${SATURN_DIR}/protein_embeddings/extract.py"
export ESM_FILTER="${SATURN_DIR}/protein_embeddings/filter_proteins_with_embeddings.py"
export ESM_CONVERT="${SATURN_DIR}/protein_embeddings/convert_protein_embeddings_to_gene_embeddings.py"
export ESM_MODEL="esm1b_t33_650M_UR50S"

# ------------------------------------------------------------
# ACTIVE SPECIES LIST
# ------------------------------------------------------------
export SPECIES_LIST=("Mlei" "Cgig" "Spur" "Crob" "Drer")

# ------------------------------------------------------------
# PER-SPECIES CONFIGURATION
# Call: get_species_config <SPECIES>
# Sets: H5AD_IN, H5AD_LABEL_COL, PROTEOME_FASTA, GTF, EMB_STRATEGY
# ------------------------------------------------------------
get_species_config() {
    SP="$1"

    case "${SP}" in

    Mlei)
        # Mnemiopsis leidyi — ctenophore
        # var_names: g-XXXX (hyphen; FASTA renamed to match)
        # .X: raw counts
        H5AD_IN="${BASE}/00_source_data/Mlei/Mlei_raw.h5ad"
        H5AD_LABEL_COL="seurat_clusters"
        PROTEOME_FASTA="${BASE}/01_proteomes/Mlei/Mlei.protein.faa"
        GTF=""
        EMB_STRATEGY="direct"
        ;;

    Cgig)
        # Crassostrea gigas — Pacific oyster (hemocytes only)
        # var_names: LOC###### (NCBI RefSeq)
        # GTF: LOC→XP_ mapping
        H5AD_IN="${BASE}/00_source_data/Cgig/Cgig.hemocytes.h5ad"
        H5AD_LABEL_COL="seurat_clusters"
        PROTEOME_FASTA="${BASE}/01_proteomes/Cgig/Cgig.protein.faa"
        GTF="${BASE}/01_proteomes/Cgig/Cgig.genomic.gtf"
        EMB_STRATEGY="gtf_loc_to_xp"
        ;;

    Spur)
        # Strongylocentrotus purpuratus — sea urchin (full dataset)
        # var_names: LOC###### (NCBI RefSeq)
        # GTF: LOC→XP_ mapping
        H5AD_IN="${BASE}/00_source_data/Spur/Spur_full.h5ad"
        H5AD_LABEL_COL="cell_type"
        PROTEOME_FASTA="${BASE}/01_proteomes/Spur/Spur_5.0.protein.faa"
        GTF="${BASE}/01_proteomes/Spur/Spur_5.0.genomic.gtf"
        EMB_STRATEGY="gtf_loc_to_xp"
        ;;

    Crob)
        # Ciona robusta — tunicate blood (Scully et al. 2025, GSE296253)
        # var_names: KY21.Chr1.N (Ghost DB KY21 gene models)
        # FASTA headers: KY21.Chr1.N.v1.nonSL19-1 → strip suffix
        H5AD_IN="${BASE}/00_source_data/Crob/Crob_blood.h5ad"
        H5AD_LABEL_COL="cell_type"
        PROTEOME_FASTA="${BASE}/01_proteomes/Crob/Crob_KY21_ghost.faa"
        GTF=""
        EMB_STRATEGY="direct_strip_suffix"
        ;;

    Drer)
        # Danio rerio — zebrafish spleen immune cells (Jiao et al. 2024, GSE186158)
        # var_names: gene symbols (ptpn12, cd79b, etc.)
        # GTF: symbol→ENSDARG→ENSDARP mapping (Ensembl GRCz11)
        H5AD_IN="${BASE}/00_source_data/Drer/Drer.h5ad"
        H5AD_LABEL_COL="seurat_clusters"
        PROTEOME_FASTA="${BASE}/01_proteomes/Drer/Drer.protein.faa"
        GTF="${BASE}/01_proteomes/Drer/Drer.genomic.gtf"
        EMB_STRATEGY="gtf_symbol_to_ensdarg"
        ;;

    Dmel)
        # Drosophila melanogaster -- larval circulating hemocytes (Cattenoz et al. 2020, E-MTAB-8698 NI)
        # var_names: FlyBase FBgn IDs. Proteome: Ensembl r113 BDGP6 pep.all (headers carry gene:FBgn...)
        # Prepared by 07d_prepare_Dmel_Hsap.sh (no 00_source_data copy; H5AD_IN points at the prepared file)
        H5AD_IN="${BASE}/02_processed_h5ad/Dmel_saturn.h5ad"
        H5AD_LABEL_COL="cell_type"
        PROTEOME_FASTA="$(ls ${BASE}/01_proteomes/Dmel/Drosophila_melanogaster.*.pep.all.fa 2>/dev/null | head -1)"
        GTF=""
        EMB_STRATEGY="ensembl_gene_id"
        ;;

    Hsap)
        # Homo sapiens -- CellTypist cross-tissue immune atlas, 10x 3' v3 subset (Dominguez Conde et al. 2022)
        # var_names: Ensembl ENSG IDs. Proteome: Ensembl r113 GRCh38 pep.all (headers carry gene:ENSG....version)
        # Prepared by 07d_prepare_Dmel_Hsap.sh
        H5AD_IN="${BASE}/02_processed_h5ad/Hsap_saturn.h5ad"
        H5AD_LABEL_COL="cell_type"
        PROTEOME_FASTA="${BASE}/01_proteomes/Hsap/Homo_sapiens.GRCh38.pep.all.fa"
        GTF=""
        EMB_STRATEGY="ensembl_gene_id"
        ;;

    *)
        echo "ERROR: Unknown species '${SP}'. Add a config block to run_config.sh." >&2
        return 1
        ;;
    esac

    export H5AD_IN H5AD_LABEL_COL PROTEOME_FASTA GTF EMB_STRATEGY

    # Derived paths — computed from SP after case block
    export H5AD_SATURN="${BASE}/02_processed_h5ad/${SP}_saturn.h5ad"
    export EMB_DIR="${BASE}/03_embeddings/${SP}"
    export EMB_GENES_TXT="${EMB_DIR}/${SP}_genes.txt"
    export EMB_FILTERED_FASTA="${EMB_DIR}/${SP}_filtered.faa"
    export EMB_PT_DIR="${EMB_DIR}/${SP}_pt_embeddings"
    export EMB_MAP_JSON="${EMB_DIR}/${SP}_gene_to_protein.json"
    export EMB_FILTERED_JSON="${EMB_DIR}/${SP}_gene_to_protein.filtered.json"
    export EMB_H5="${EMB_DIR}/${SP}_gene_embeddings.h5"
    export EMB_PT="${EMB_DIR}/${SP}_gene_embeddings.pt"
    export EMB_LOC_MAP_TSV="${EMB_DIR}/${SP}_gene_to_protein.tsv"
}

# ------------------------------------------------------------
# SATURN HYPERPARAMETERS
# ------------------------------------------------------------
export K_MACROGENES=3000
export MACROGENE_SEED=42
export HV_GENES=8000

# Post-SATURN UMAP (settled: no StandardScaler, PCA5, cosine)
export UMAP_N_PCS=5
export UMAP_METRIC="cosine"
export UMAP_N_NEIGHBORS=80
export UMAP_MIN_DIST=0.10
export LEIDEN_RESOLUTION=1.2
export MIN_CELLS_PER_SPECIES=50   # min Mlei cells per cluster

# MIN-N macrogene filters — both run in postSATURN step
export MIN_N_FILTERS="1 3 5"        # space-separated; runs MIN3-of-5 and MIN5-of-5

# ------------------------------------------------------------
# CONDA ENVIRONMENT
# ------------------------------------------------------------
export CONDA_ENV="saturn_protein"
