#!/usr/bin/env bash
# 06a0_download_outgroup.sh
# Outgroup proteome for the 6-species OrthoFinder run (06a --species ...,Sros; 06b6):
#   Salpingoeca rosetta (choanoflagellate), assembly Proterospongia_sp_ATCC50818 (Fairclough et al. 2013).
#   Tries the NCBI RefSeq annotation first (GCF_000188695.1), then the GenBank one (GCA_000188695.1).
# Why an outgroup: the 4-species OrthoFinder 3.1 run wrote no N0.tsv (only N1/N2). With a non-metazoan
#   outgroup the root is Choanozoa and the Metazoa-level HOGs are a non-root node (expected N1).
# Writes 01_proteomes/Sros/Sros.protein.faa and Sros.genomic.gff (unzipped), plus the .gz originals.
# Runs on the 'transfers' queue (as 07e). Resumable (wget -c).
# Submit:  bsub < /scratch/dark_genes/SATURN_Mnemi/scripts/06a0_download_outgroup.sh
# If both URLs fail (NCBI path change), download the protein FASTA + GFF3 by hand from
#   https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000188695.1/  and save them under the same two names.
#BSUB -J dl_sros
#BSUB -q transfers
#BSUB -P dark_genes
#BSUB -n 1
#BSUB -W 1:00
#BSUB -o /scratch/dark_genes/SATURN_Mnemi/logs/dl_sros.%J.out
#BSUB -e /scratch/dark_genes/SATURN_Mnemi/logs/dl_sros.%J.err

set -uo pipefail
BASE="/scratch/dark_genes/SATURN_Mnemi"
DST="${BASE}/01_proteomes/Sros"
mkdir -p "${DST}" "${BASE}/logs"
ASM="Proterospongia_sp_ATCC50818"
echo "=== start $(date) on $(hostname) ==="

ok=0
for ACC in GCF_000188695.1 GCA_000188695.1; do
    D="https://ftp.ncbi.nlm.nih.gov/genomes/all/${ACC:0:3}/${ACC:4:3}/${ACC:7:3}/${ACC:10:3}/${ACC}_${ASM}"
    echo "--- trying ${D}"
    if wget -c -nv --tries=5 --timeout=60 -P "${DST}" "${D}/${ACC}_${ASM}_protein.faa.gz" && \
       wget -c -nv --tries=5 --timeout=60 -P "${DST}" "${D}/${ACC}_${ASM}_genomic.gff.gz"; then
        gzip -t "${DST}/${ACC}_${ASM}_protein.faa.gz" && gzip -t "${DST}/${ACC}_${ASM}_genomic.gff.gz" || { echo "corrupt gz"; continue; }
        gunzip -c "${DST}/${ACC}_${ASM}_protein.faa.gz" > "${DST}/Sros.protein.faa"
        gunzip -c "${DST}/${ACC}_${ASM}_genomic.gff.gz"  > "${DST}/Sros.genomic.gff"
        echo "${ACC}_${ASM}" > "${DST}/SOURCE.txt"
        ok=1; break
    fi
done
[[ ${ok} -eq 1 ]] || { echo "ERROR: both downloads failed -- see header for the manual route" >&2; exit 1; }

echo; echo "=== Sros summary (source: $(cat "${DST}/SOURCE.txt")) ==="
echo "proteins:            $(grep -c '^>' "${DST}/Sros.protein.faa")"
echo "CDS lines w/ protein_id: $(awk -F'\t' '$3=="CDS"' "${DST}/Sros.genomic.gff" | grep -c 'protein_id=')"
echo "distinct locus_tags: $(awk -F'\t' '$3=="CDS"' "${DST}/Sros.genomic.gff" | grep -o 'locus_tag=[^;]*' | sort -u | wc -l)"
ls -lh "${DST}"
echo "=== done $(date) ==="
