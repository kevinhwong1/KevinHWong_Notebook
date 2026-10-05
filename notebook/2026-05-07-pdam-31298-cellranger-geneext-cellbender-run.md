---
title: "Pdam 31298: Cell Ranger → GeneExt → CellBender run"
date: 2026-05-07
type: "Analysis"
project: "Cnidarian Stem Cells"
categories: ["Analysis", "Pocillopora damicornis", "scRNA-seq"]
tools: ["Cell Ranger", "GeneExt", "CellBender", "samtools"]
---

::: {.callout-tip}
This is the run log. The cleaned-up, reusable version of these scripts is the pipeline
[10x scRNA-seq: Cell Ranger → GeneExt → CellBender](../pipelines/scrnaseq-10x-geneext-cellbender.qmd).
:::

- **Libraries:** 31298-001 (all cells) and 31298-002 (ALDH+), *Pocillopora damicornis*.
- **Result:** GeneExt extended 13,283 / 25,422 genes (median extension 1,297 bp).
- **Hurdles:** `samtools merge` failed with the old samtools (updated and resubmitted steps 2–6); GeneExt had to be run from its install directory after clearing `tmp/` (resubmitted steps 3–6).

# Make script

```bash
cd /nethome/kxw755/scripts
nano pdam_31298_pipeline.sh
```

```bash
#!/bin/bash
# =============================================================================
# Pdam 31298 Full scRNA-seq Pipeline
# Libraries: 31298-001 (All Cells) | 31298-002 (ALDH+)
#
# Steps:
#   1a. cellranger count 001      (All Cells, original ref)
#   1b. cellranger count 002      (ALDH+,     original ref)
#   2.  samtools merge            (merge BAMs for GeneExt)
#   3.  GeneExt                   (extend gene models)
#   4.  cellranger mkref_geneext  (GeneExt-extended ref)
#   5a. cellranger count 001      (All Cells, geneext ref)
#   5b. cellranger count 002      (ALDH+,     geneext ref)
#   6a. CellBender 001            (ambient RNA removal)
#   6b. CellBender 002
#
# NOTE: cellranger mkref is SKIPPED - existing pdam ref used from:
#   /nethome/kxw755/genomes/pdam_genome/pdam
#
# Usage:
#   bash pdam_31298_pipeline.sh
#
# Each step is submitted as a dependent LSF job and will only run
# after the previous step completes successfully.
# =============================================================================

set -euo pipefail

# =============================================================================
# PATHS
# =============================================================================

GENOME_DIR="/nethome/kxw755/genomes/pdam_genome"
FASTA="${GENOME_DIR}/GCA_003704095.1_ASM370409v1_genomic.fa"
GTF="${GENOME_DIR}/pdam_exons_only.validated.gtf"

FASTQ_001="/nethome/kxw755/Ehrens-31298-001_Pdam"
FASTQ_002="/nethome/kxw755/Ehrens-31298-002_Pdam_ALDH"
SAMPLE_001="Ehrens-31298-001_GEX3"
SAMPLE_002="Ehrens-31298-002_GEX3"

SCRATCH="/scratch/projects/dark_genes/GeneExt/Pdam_31298"
GENEEXT_DIR="${SCRATCH}/GeneExt"
REF_ORIG="${GENOME_DIR}/pdam"
REF_GENEEXT="${SCRATCH}/pdam_31298_geneext"
GENEEXT_GTF="${SCRATCH}/pdam_31298_GeneExt.gtf"

GENEEXT_PY="/nethome/kxw755/GeneExt/geneext.py"

EMAIL="kxw755@earth.miami.edu"
PROJECT="dark_genes"

mkdir -p "${SCRATCH}" "${GENEEXT_DIR}"

echo "========================================================"
echo " Pdam 31298 Pipeline - submitting LSF job chain"
echo " $(date)"
echo "========================================================"

# =============================================================================
# STEP 1a: cellranger count - 31298-001 All Cells (original ref)
# mkref is SKIPPED: existing ref at ${GENOME_DIR}/pdam
# =============================================================================

cat > /tmp/step1a_count_001.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_count_001
#BSUB -q bigmem
#BSUB -P ${PROJECT}
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u ${EMAIL}
#BSUB -o ${SCRATCH}/step1a_count_001_%J.out
#BSUB -e ${SCRATCH}/step1a_count_001_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cd ${SCRATCH}

cellranger count \\
    --id=31298-001 \\
    --transcriptome=${REF_ORIG} \\
    --fastqs=${FASTQ_001} \\
    --sample=${SAMPLE_001} \\
    --create-bam=true
ENDJOB

JID_COUNT_001=$(bsub < /tmp/step1a_count_001.job | grep -oP '(?<=Job <)\d+')
echo "Step 1a submitted - count 001       job ID: ${JID_COUNT_001}"

# =============================================================================
# STEP 1b: cellranger count - 31298-002 ALDH+ (original ref)
# =============================================================================

cat > /tmp/step1b_count_002.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_count_002
#BSUB -q bigmem
#BSUB -P ${PROJECT}
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u ${EMAIL}
#BSUB -o ${SCRATCH}/step1b_count_002_%J.out
#BSUB -e ${SCRATCH}/step1b_count_002_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cd ${SCRATCH}

cellranger count \\
    --id=31298-002 \\
    --transcriptome=${REF_ORIG} \\
    --fastqs=${FASTQ_002} \\
    --sample=${SAMPLE_002} \\
    --create-bam=true
ENDJOB

JID_COUNT_002=$(bsub < /tmp/step1b_count_002.job | grep -oP '(?<=Job <)\d+')
echo "Step 1b submitted - count 002       job ID: ${JID_COUNT_002}"

# =============================================================================
# STEP 2: samtools merge (waits for both initial counts to finish)
# =============================================================================

cat > /tmp/step2_merge.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_bam_merge
#BSUB -q bigmem
#BSUB -P ${PROJECT}
#BSUB -n 16
#BSUB -W 24:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_COUNT_001}) && done(${JID_COUNT_002})"
#BSUB -u ${EMAIL}
#BSUB -o ${SCRATCH}/step2_merge_%J.out
#BSUB -e ${SCRATCH}/step2_merge_%J.err
#BSUB -B
#BSUB -N
###################################################################

source ~/anaconda3/bin/activate
conda activate samtools_env

cp ${SCRATCH}/31298-001/outs/possorted_genome_bam.bam \
   ${GENEEXT_DIR}/31298-001_possorted_genome_bam.bam

cp ${SCRATCH}/31298-002/outs/possorted_genome_bam.bam \
   ${GENEEXT_DIR}/31298-002_possorted_genome_bam.bam

cd ${GENEEXT_DIR}

samtools merge merge.bam \
    31298-001_possorted_genome_bam.bam \
    31298-002_possorted_genome_bam.bam

echo "BAM line counts:"
wc -l 31298-001_possorted_genome_bam.bam
wc -l 31298-002_possorted_genome_bam.bam
wc -l merge.bam
ENDJOB

JID_MERGE=$(bsub < /tmp/step2_merge.job | grep -oP '(?<=Job <)\d+')
echo "Step 2  submitted - bam merge       job ID: ${JID_MERGE}"

# =============================================================================
# STEP 3: GeneExt
# =============================================================================

cat > /tmp/step3_geneext.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_GeneExt
#BSUB -q bigmem
#BSUB -P ${PROJECT}
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_MERGE})"
#BSUB -u ${EMAIL}
#BSUB -o ${SCRATCH}/step3_geneext_%J.out
#BSUB -e ${SCRATCH}/step3_geneext_%J.err
#BSUB -N
###################################################################

source ~/anaconda3/bin/activate
conda activate /projectnb/pegasus/nethome/kxw755/geneext

# --clip_strand both prevents GeneExt from creating overlaps on the same strand

python ${GENEEXT_PY} \
    -g ${GTF} \
    -b ${GENEEXT_DIR}/merge.bam \
    -o ${GENEEXT_GTF} \
    -v 3 \
    -j 8 \
    --clip_strand both

echo "GeneExt log (last 20 lines):"
tail -20 ${GENEEXT_GTF}.GeneExt.log
ENDJOB

JID_GENEEXT=$(bsub < /tmp/step3_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 3  submitted - GeneExt         job ID: ${JID_GENEEXT}"

# =============================================================================
# STEP 4: cellranger mkref with GeneExt GTF
# =============================================================================

cat > /tmp/step4_mkref_geneext.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_mkref_geneext
#BSUB -q bigmem
#BSUB -P ${PROJECT}
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_GENEEXT})"
#BSUB -u ${EMAIL}
#BSUB -o ${SCRATCH}/step4_mkref_geneext_%J.out
#BSUB -e ${SCRATCH}/step4_mkref_geneext_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cd ${SCRATCH}

cellranger mkref \
    --genome=pdam_31298_geneext \
    --fasta=${FASTA} \
    --genes=${GENEEXT_GTF}
ENDJOB

JID_MKREF_GE=$(bsub < /tmp/step4_mkref_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 4  submitted - mkref geneext   job ID: ${JID_MKREF_GE}"

# =============================================================================
# STEP 5a: cellranger count - 31298-001 All Cells (geneext ref)
# =============================================================================

cat > /tmp/step5a_count_001_geneext.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_count_001_geneext
#BSUB -q bigmem
#BSUB -P ${PROJECT}
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_MKREF_GE})"
#BSUB -u ${EMAIL}
#BSUB -o ${SCRATCH}/step5a_count_001_geneext_%J.out
#BSUB -e ${SCRATCH}/step5a_count_001_geneext_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cd ${SCRATCH}

cellranger count \
    --id=31298-001_geneext \
    --transcriptome=${REF_GENEEXT} \
    --fastqs=${FASTQ_001} \
    --sample=${SAMPLE_001} \
    --create-bam=true
ENDJOB

JID_COUNT_001_GE=$(bsub < /tmp/step5a_count_001_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 5a submitted - count 001 GE    job ID: ${JID_COUNT_001_GE}"

# =============================================================================
# STEP 5b: cellranger count - 31298-002 ALDH+ (geneext ref)
# =============================================================================

cat > /tmp/step5b_count_002_geneext.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_count_002_geneext
#BSUB -q bigmem
#BSUB -P ${PROJECT}
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_MKREF_GE})"
#BSUB -u ${EMAIL}
#BSUB -o ${SCRATCH}/step5b_count_002_geneext_%J.out
#BSUB -e ${SCRATCH}/step5b_count_002_geneext_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cd ${SCRATCH}

cellranger count \
    --id=31298-002_geneext \
    --transcriptome=${REF_GENEEXT} \
    --fastqs=${FASTQ_002} \
    --sample=${SAMPLE_002} \
    --create-bam=true
ENDJOB

JID_COUNT_002_GE=$(bsub < /tmp/step5b_count_002_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 5b submitted - count 002 GE    job ID: ${JID_COUNT_002_GE}"

# =============================================================================
# STEP 6a: CellBender - 31298-001 All Cells
# =============================================================================

cat > /tmp/step6a_cellbender_001.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_cellbender_001
#BSUB -q bigmem
#BSUB -P ${PROJECT}
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_COUNT_001_GE})"
#BSUB -u ${EMAIL}
#BSUB -o ${SCRATCH}/step6a_cellbender_001_%J.out
#BSUB -e ${SCRATCH}/step6a_cellbender_001_%J.err
#BSUB -B
#BSUB -N
###################################################################

source ~/anaconda3/bin/activate
conda activate CellBender

cellbender remove-background \
    --input  ${SCRATCH}/31298-001_geneext/outs/raw_feature_bc_matrix.h5 \
    --output ${SCRATCH}/31298-001_geneext/31298-001_geneext_cellbender_output.h5 \
    --fpr 0.01 \
    --epochs 150 \
    --learning-rate 0.00005
ENDJOB

JID_CB_001=$(bsub < /tmp/step6a_cellbender_001.job | grep -oP '(?<=Job <)\d+')
echo "Step 6a submitted - CellBender 001  job ID: ${JID_CB_001}"

# =============================================================================
# STEP 6b: CellBender - 31298-002 ALDH+
# =============================================================================

cat > /tmp/step6b_cellbender_002.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_cellbender_002
#BSUB -q bigmem
#BSUB -P ${PROJECT}
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_COUNT_002_GE})"
#BSUB -u ${EMAIL}
#BSUB -o ${SCRATCH}/step6b_cellbender_002_%J.out
#BSUB -e ${SCRATCH}/step6b_cellbender_002_%J.err
#BSUB -B
#BSUB -N
###################################################################

source ~/anaconda3/bin/activate
conda activate CellBender

cellbender remove-background \
    --input  ${SCRATCH}/31298-002_geneext/outs/raw_feature_bc_matrix.h5 \
    --output ${SCRATCH}/31298-002_geneext/31298-002_geneext_cellbender_output.h5 \
    --fpr 0.01 \
    --epochs 150 \
    --learning-rate 0.00005
ENDJOB

JID_CB_002=$(bsub < /tmp/step6b_cellbender_002.job | grep -oP '(?<=Job <)\d+')
echo "Step 6b submitted - CellBender 002  job ID: ${JID_CB_002}"

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "========================================================"
echo " All jobs submitted successfully!"
echo "========================================================"
echo " Job chain:"
echo "   Step 1a count 001        : ${JID_COUNT_001}"
echo "   Step 1b count 002        : ${JID_COUNT_002}"
echo "   Step 2  bam merge        : ${JID_MERGE}  (waits on 1a && 1b)"
echo "   Step 3  GeneExt          : ${JID_GENEEXT}  (waits on 2)"
echo "   Step 4  mkref geneext    : ${JID_MKREF_GE}  (waits on 3)"
echo "   Step 5a count 001 GE     : ${JID_COUNT_001_GE}  (waits on 4)"
echo "   Step 5b count 002 GE     : ${JID_COUNT_002_GE}  (waits on 4)"
echo "   Step 6a CellBender 001   : ${JID_CB_001}  (waits on 5a)"
echo "   Step 6b CellBender 002   : ${JID_CB_002}  (waits on 5b)"
echo ""
echo " Monitor with: bjobs"
echo " Check logs in: ${SCRATCH}/"
echo ""
echo " Final outputs to scp:"
echo "   ${SCRATCH}/31298-001_geneext/31298-001_geneext_cellbender_output_filtered.h5"
echo "   ${SCRATCH}/31298-001_geneext/31298-001_geneext_cellbender_output_metrics.csv"
echo "   ${SCRATCH}/31298-001_geneext/31298-001_geneext_cellbender_output_report.html"
echo "   ${SCRATCH}/31298-002_geneext/31298-002_geneext_cellbender_output_filtered.h5"
echo "   ${SCRATCH}/31298-002_geneext/31298-002_geneext_cellbender_output_metrics.csv"
echo "   ${SCRATCH}/31298-002_geneext/31298-002_geneext_cellbender_output_report.html"
echo "========================================================"
```

`bash pdam_31298_pipeline.sh`

Samtools had an error, so I had to update samtools and re-run post count: 

`nano pdam_31298_pipeline_post_count.sh`

```bash
#!/bin/bash
# =============================================================================
# Pdam 31298 Pipeline - Steps 2-6 only (post-count resubmission)
# Run this after 31298-001 and 31298-002 cellranger counts are complete
#
# Steps:
#   2.  samtools merge            (merge BAMs for GeneExt)
#   3.  GeneExt                   (extend gene models)
#   4.  cellranger mkref_geneext  (GeneExt-extended ref)
#   5a. cellranger count 001      (All Cells, geneext ref)
#   5b. cellranger count 002      (ALDH+,     geneext ref)
#   6a. CellBender 001            (ambient RNA removal)
#   6b. CellBender 002
#
# Usage:
#   bash pdam_31298_pipeline_post_count.sh
# =============================================================================

set -euo pipefail

# =============================================================================
# PATHS
# =============================================================================

GENOME_DIR="/nethome/kxw755/genomes/pdam_genome"
FASTA="${GENOME_DIR}/GCA_003704095.1_ASM370409v1_genomic.fa"
GTF="${GENOME_DIR}/pdam_exons_only.validated.gtf"

FASTQ_001="/nethome/kxw755/Ehrens-31298-001_Pdam"
FASTQ_002="/nethome/kxw755/Ehrens-31298-002_Pdam_ALDH"
SAMPLE_001="Ehrens-31298-001_GEX3"
SAMPLE_002="Ehrens-31298-002_GEX3"

SCRATCH="/scratch/projects/dark_genes/GeneExt/Pdam_31298"
GENEEXT_DIR="${SCRATCH}/GeneExt"
REF_GENEEXT="${SCRATCH}/pdam_31298_geneext"
GENEEXT_GTF="${SCRATCH}/pdam_31298_GeneExt.gtf"

GENEEXT_PY="/nethome/kxw755/GeneExt/geneext.py"

EMAIL="kxw755@earth.miami.edu"
PROJECT="dark_genes"

mkdir -p "${SCRATCH}" "${GENEEXT_DIR}"

echo "========================================================"
echo " Pdam 31298 Pipeline - Steps 2-6 (post-count)"
echo " $(date)"
echo "========================================================"

# =============================================================================
# STEP 2: samtools merge
# =============================================================================

cat > /tmp/step2_merge.job << 'ENDJOB'
#!/bin/bash
#BSUB -J pdam31298_bam_merge
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 24:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step2_merge_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step2_merge_%J.err
#BSUB -B
#BSUB -N
###################################################################

source ~/anaconda3/bin/activate
conda activate samtools_env

cp /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001/outs/possorted_genome_bam.bam \
   /scratch/projects/dark_genes/GeneExt/Pdam_31298/GeneExt/31298-001_possorted_genome_bam.bam

cp /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002/outs/possorted_genome_bam.bam \
   /scratch/projects/dark_genes/GeneExt/Pdam_31298/GeneExt/31298-002_possorted_genome_bam.bam

cd /scratch/projects/dark_genes/GeneExt/Pdam_31298/GeneExt

samtools merge merge.bam \
    31298-001_possorted_genome_bam.bam \
    31298-002_possorted_genome_bam.bam

echo "Merge complete. BAM sizes:"
ls -lh /scratch/projects/dark_genes/GeneExt/Pdam_31298/GeneExt/*.bam
ENDJOB

JID_MERGE=$(bsub < /tmp/step2_merge.job | grep -oP '(?<=Job <)\d+')
echo "Step 2  submitted - bam merge       job ID: ${JID_MERGE}"

# =============================================================================
# STEP 3: GeneExt
# =============================================================================

cat > /tmp/step3_geneext.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_GeneExt
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_MERGE})"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step3_geneext_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step3_geneext_%J.err
#BSUB -N
###################################################################

source ~/anaconda3/bin/activate
conda activate geneext

# --clip_strand both prevents GeneExt from creating overlaps on the same strand

python /nethome/kxw755/GeneExt/geneext.py \
    -g /nethome/kxw755/genomes/pdam_genome/pdam_exons_only.validated.gtf \
    -b /scratch/projects/dark_genes/GeneExt/Pdam_31298/GeneExt/merge.bam \
    -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_GeneExt.gtf \
    -v 3 \
    -j 8 \
    --clip_strand both

echo "GeneExt log (last 20 lines):"
tail -20 /scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_GeneExt.gtf.GeneExt.log
ENDJOB

JID_GENEEXT=$(bsub < /tmp/step3_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 3  submitted - GeneExt         job ID: ${JID_GENEEXT}"

# =============================================================================
# STEP 4: cellranger mkref with GeneExt GTF
# =============================================================================

cat > /tmp/step4_mkref_geneext.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_mkref_geneext
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_GENEEXT})"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step4_mkref_geneext_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step4_mkref_geneext_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cd /scratch/projects/dark_genes/GeneExt/Pdam_31298

cellranger mkref \
    --genome=pdam_31298_geneext \
    --fasta=/nethome/kxw755/genomes/pdam_genome/GCA_003704095.1_ASM370409v1_genomic.fa \
    --genes=/scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_GeneExt.gtf
ENDJOB

JID_MKREF_GE=$(bsub < /tmp/step4_mkref_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 4  submitted - mkref geneext   job ID: ${JID_MKREF_GE}"

# =============================================================================
# STEP 5a: cellranger count - 31298-001 All Cells (geneext ref)
# =============================================================================

cat > /tmp/step5a_count_001_geneext.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_count_001_geneext
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_MKREF_GE})"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step5a_count_001_geneext_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step5a_count_001_geneext_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cd /scratch/projects/dark_genes/GeneExt/Pdam_31298

cellranger count \
    --id=31298-001_geneext \
    --transcriptome=/scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_geneext \
    --fastqs=/nethome/kxw755/Ehrens-31298-001_Pdam \
    --sample=Ehrens-31298-001_GEX3 \
    --create-bam=true
ENDJOB

JID_COUNT_001_GE=$(bsub < /tmp/step5a_count_001_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 5a submitted - count 001 GE    job ID: ${JID_COUNT_001_GE}"

# =============================================================================
# STEP 5b: cellranger count - 31298-002 ALDH+ (geneext ref)
# =============================================================================

cat > /tmp/step5b_count_002_geneext.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_count_002_geneext
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_MKREF_GE})"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step5b_count_002_geneext_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step5b_count_002_geneext_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cd /scratch/projects/dark_genes/GeneExt/Pdam_31298

cellranger count \
    --id=31298-002_geneext \
    --transcriptome=/scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_geneext \
    --fastqs=/nethome/kxw755/Ehrens-31298-002_Pdam_ALDH \
    --sample=Ehrens-31298-002_GEX3 \
    --create-bam=true
ENDJOB

JID_COUNT_002_GE=$(bsub < /tmp/step5b_count_002_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 5b submitted - count 002 GE    job ID: ${JID_COUNT_002_GE}"

# =============================================================================
# STEP 6a: CellBender - 31298-001 All Cells
# =============================================================================

cat > /tmp/step6a_cellbender_001.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_cellbender_001
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_COUNT_001_GE})"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step6a_cellbender_001_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step6a_cellbender_001_%J.err
#BSUB -B
#BSUB -N
###################################################################

source ~/anaconda3/bin/activate
conda activate CellBender

cellbender remove-background \
    --input  /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/outs/raw_feature_bc_matrix.h5 \
    --output /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/31298-001_geneext_cellbender_output.h5 \
    --fpr 0.01 \
    --epochs 150 \
    --learning-rate 0.00005
ENDJOB

JID_CB_001=$(bsub < /tmp/step6a_cellbender_001.job | grep -oP '(?<=Job <)\d+')
echo "Step 6a submitted - CellBender 001  job ID: ${JID_CB_001}"

# =============================================================================
# STEP 6b: CellBender - 31298-002 ALDH+
# =============================================================================

cat > /tmp/step6b_cellbender_002.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_cellbender_002
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_COUNT_002_GE})"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step6b_cellbender_002_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step6b_cellbender_002_%J.err
#BSUB -B
#BSUB -N
###################################################################

source ~/anaconda3/bin/activate
conda activate CellBender

cellbender remove-background \
    --input  /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/outs/raw_feature_bc_matrix.h5 \
    --output /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/31298-002_geneext_cellbender_output.h5 \
    --fpr 0.01 \
    --epochs 150 \
    --learning-rate 0.00005
ENDJOB

JID_CB_002=$(bsub < /tmp/step6b_cellbender_002.job | grep -oP '(?<=Job <)\d+')
echo "Step 6b submitted - CellBender 002  job ID: ${JID_CB_002}"

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "========================================================"
echo " All jobs submitted successfully!"
echo "========================================================"
echo " Job chain:"
echo "   Step 2  bam merge        : ${JID_MERGE}"
echo "   Step 3  GeneExt          : ${JID_GENEEXT}  (waits on 2)"
echo "   Step 4  mkref geneext    : ${JID_MKREF_GE}  (waits on 3)"
echo "   Step 5a count 001 GE     : ${JID_COUNT_001_GE}  (waits on 4)"
echo "   Step 5b count 002 GE     : ${JID_COUNT_002_GE}  (waits on 4)"
echo "   Step 6a CellBender 001   : ${JID_CB_001}  (waits on 5a)"
echo "   Step 6b CellBender 002   : ${JID_CB_002}  (waits on 5b)"
echo ""
echo " Monitor with: bjobs"
echo " Check logs in: /scratch/projects/dark_genes/GeneExt/Pdam_31298/"
echo ""
echo " Final outputs to scp:"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/31298-001_geneext_cellbender_output_filtered.h5"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/31298-001_geneext_cellbender_output_metrics.csv"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/31298-001_geneext_cellbender_output_report.html"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/31298-002_geneext_cellbender_output_filtered.h5"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/31298-002_geneext_cellbender_output_metrics.csv"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/31298-002_geneext_cellbender_output_report.html"
echo "========================================================"
```

Export prelim CellRanger outputs
```bash
scp kxw755@pegasus2.ccs.miami.edu:/scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001/outs/web_summary.html ./Pdam_31298-001_web_summary.html
scp kxw755@pegasus2.ccs.miami.edu:/scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002/outs/web_summary.html ./Pdam_31298-002_web_summary.html
```


* remove any tmp files before running geneext 

`rm -rf /nethome/kxw755/GeneExt/tmp`

`nano pdam_31298_pipeline_post_merge.sh`

```bash
#!/bin/bash
# =============================================================================
# Pdam 31298 Pipeline - Steps 3-6 only (post-merge resubmission)
# Run this after samtools merge is complete
#
# Steps:
#   3.  GeneExt                   (extend gene models)
#   4.  cellranger mkref_geneext  (GeneExt-extended ref)
#   5a. cellranger count 001      (All Cells, geneext ref)
#   5b. cellranger count 002      (ALDH+,     geneext ref)
#   6a. CellBender 001            (ambient RNA removal)
#   6b. CellBender 002
#
# Usage:
#   bash pdam_31298_pipeline_post_merge.sh
# =============================================================================

set -euo pipefail

SCRATCH="/scratch/projects/dark_genes/GeneExt/Pdam_31298"

echo "========================================================"
echo " Pdam 31298 Pipeline - Steps 3-6 (post-merge)"
echo " $(date)"
echo "========================================================"

# =============================================================================
# STEP 3: GeneExt
# =============================================================================

cat > /tmp/step3_geneext.job << 'ENDJOB'
#!/bin/bash
#BSUB -J pdam31298_GeneExt
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step3_geneext_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step3_geneext_%J.err
#BSUB -N
###################################################################

source ~/anaconda3/bin/activate
conda activate geneext

# run from GeneExt install dir so tmp/ is created in the right place
cd /nethome/kxw755/GeneExt

# --clip_strand both prevents GeneExt from creating overlaps on the same strand

python /nethome/kxw755/GeneExt/geneext.py \
    -g /nethome/kxw755/genomes/pdam_genome/pdam_exons_only.validated.gtf \
    -b /scratch/projects/dark_genes/GeneExt/Pdam_31298/GeneExt/merge.bam \
    -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_GeneExt.gtf \
    -v 3 \
    -j 8 \
    --clip_strand both

echo "GeneExt log (last 20 lines):"
tail -20 /scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_GeneExt.gtf.GeneExt.log
ENDJOB

JID_GENEEXT=$(bsub < /tmp/step3_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 3  submitted - GeneExt         job ID: ${JID_GENEEXT}"

# =============================================================================
# STEP 4: cellranger mkref with GeneExt GTF
# =============================================================================

cat > /tmp/step4_mkref_geneext.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_mkref_geneext
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_GENEEXT})"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step4_mkref_geneext_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step4_mkref_geneext_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cd /scratch/projects/dark_genes/GeneExt/Pdam_31298

cellranger mkref \
    --genome=pdam_31298_geneext \
    --fasta=/nethome/kxw755/genomes/pdam_genome/GCA_003704095.1_ASM370409v1_genomic.fa \
    --genes=/scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_GeneExt.gtf
ENDJOB

JID_MKREF_GE=$(bsub < /tmp/step4_mkref_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 4  submitted - mkref geneext   job ID: ${JID_MKREF_GE}"

# =============================================================================
# STEP 5a: cellranger count - 31298-001 All Cells (geneext ref)
# =============================================================================

cat > /tmp/step5a_count_001_geneext.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_count_001_geneext
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_MKREF_GE})"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step5a_count_001_geneext_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step5a_count_001_geneext_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cd /scratch/projects/dark_genes/GeneExt/Pdam_31298

cellranger count \
    --id=31298-001_geneext \
    --transcriptome=/scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_geneext \
    --fastqs=/nethome/kxw755/Ehrens-31298-001_Pdam \
    --sample=Ehrens-31298-001_GEX3 \
    --create-bam=true
ENDJOB

JID_COUNT_001_GE=$(bsub < /tmp/step5a_count_001_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 5a submitted - count 001 GE    job ID: ${JID_COUNT_001_GE}"

# =============================================================================
# STEP 5b: cellranger count - 31298-002 ALDH+ (geneext ref)
# =============================================================================

cat > /tmp/step5b_count_002_geneext.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_count_002_geneext
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_MKREF_GE})"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step5b_count_002_geneext_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step5b_count_002_geneext_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cd /scratch/projects/dark_genes/GeneExt/Pdam_31298

cellranger count \
    --id=31298-002_geneext \
    --transcriptome=/scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_geneext \
    --fastqs=/nethome/kxw755/Ehrens-31298-002_Pdam_ALDH \
    --sample=Ehrens-31298-002_GEX3 \
    --create-bam=true
ENDJOB

JID_COUNT_002_GE=$(bsub < /tmp/step5b_count_002_geneext.job | grep -oP '(?<=Job <)\d+')
echo "Step 5b submitted - count 002 GE    job ID: ${JID_COUNT_002_GE}"

# =============================================================================
# STEP 6a: CellBender - 31298-001 All Cells
# =============================================================================

cat > /tmp/step6a_cellbender_001.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_cellbender_001
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_COUNT_001_GE})"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step6a_cellbender_001_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step6a_cellbender_001_%J.err
#BSUB -B
#BSUB -N
###################################################################

source ~/anaconda3/bin/activate
conda activate CellBender

cellbender remove-background \
    --input  /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/outs/raw_feature_bc_matrix.h5 \
    --output /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/31298-001_geneext_cellbender_output.h5 \
    --fpr 0.01 \
    --epochs 150 \
    --learning-rate 0.00005
ENDJOB

JID_CB_001=$(bsub < /tmp/step6a_cellbender_001.job | grep -oP '(?<=Job <)\d+')
echo "Step 6a submitted - CellBender 001  job ID: ${JID_CB_001}"

# =============================================================================
# STEP 6b: CellBender - 31298-002 ALDH+
# =============================================================================

cat > /tmp/step6b_cellbender_002.job << ENDJOB
#!/bin/bash
#BSUB -J pdam31298_cellbender_002
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -w "done(${JID_COUNT_002_GE})"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o /scratch/projects/dark_genes/GeneExt/Pdam_31298/step6b_cellbender_002_%J.out
#BSUB -e /scratch/projects/dark_genes/GeneExt/Pdam_31298/step6b_cellbender_002_%J.err
#BSUB -B
#BSUB -N
###################################################################

source ~/anaconda3/bin/activate
conda activate CellBender

cellbender remove-background \
    --input  /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/outs/raw_feature_bc_matrix.h5 \
    --output /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/31298-002_geneext_cellbender_output.h5 \
    --fpr 0.01 \
    --epochs 150 \
    --learning-rate 0.00005
ENDJOB

JID_CB_002=$(bsub < /tmp/step6b_cellbender_002.job | grep -oP '(?<=Job <)\d+')
echo "Step 6b submitted - CellBender 002  job ID: ${JID_CB_002}"

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "========================================================"
echo " All jobs submitted successfully!"
echo "========================================================"
echo " Job chain:"
echo "   Step 3  GeneExt          : ${JID_GENEEXT}"
echo "   Step 4  mkref geneext    : ${JID_MKREF_GE}  (waits on 3)"
echo "   Step 5a count 001 GE     : ${JID_COUNT_001_GE}  (waits on 4)"
echo "   Step 5b count 002 GE     : ${JID_COUNT_002_GE}  (waits on 4)"
echo "   Step 6a CellBender 001   : ${JID_CB_001}  (waits on 5a)"
echo "   Step 6b CellBender 002   : ${JID_CB_002}  (waits on 5b)"
echo ""
echo " Monitor with: bjobs"
echo " Check logs in: /scratch/projects/dark_genes/GeneExt/Pdam_31298/"
echo ""
echo " Final outputs to scp:"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/31298-001_geneext_cellbender_output_filtered.h5"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/31298-001_geneext_cellbender_output_metrics.csv"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/31298-001_geneext_cellbender_output_report.html"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/31298-002_geneext_cellbender_output_filtered.h5"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/31298-002_geneext_cellbender_output_metrics.csv"
echo "   /scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/31298-002_geneext_cellbender_output_report.html"
echo "========================================================"
```

GeneExt results
```bash
╭───────────╮
│ All done! │
╰───────────╯
Extended 13283/25422 genes
Median extension length: 1297.0 bp
Running:
        Rscript geneext/plot_extensions.R tmp/extensions.tsv /scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_GeneExt.gtf.extension_length.pdf
Running:
        Rscript geneext/peak_density.R tmp/genic_peaks.bed tmp/allpeaks_noov.bed /scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_GeneExt.gtf.peak_coverage.pdf 25 
Removing tmp/_genes_peaks_closest
Removing tmp/minus.bam
Removing tmp/plus.bam
Removing tmp/_genes_tmp
Removing tmp/_peaks_tmp
Removing tmp/_peaks_tmp_sorted
Removing tmp/_genes_tmp_sorted
GeneExt log (last 20 lines):
mRNA with the most downstream exon: unassigned_transcript_25422
adding unassigned_transcript_25422: [13283/13284]
        Extended genes written: /scratch/projects/dark_genes/GeneExt/Pdam_31298/pdam_31298_GeneExt.gtf
done
```

Export prelim CellRanger outputs
```bash
scp kxw755@pegasus2.ccs.miami.edu:/scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/outs/web_summary.html ./Pdam_31298-001_web_summary_geneext.html
scp kxw755@pegasus2.ccs.miami.edu:/scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/outs/web_summary.html ./Pdam_31298-002_web_summary_geneext.html
```

Export CellBender Outputs

```bash
scp kxw755@pegasus2.ccs.miami.edu:/scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/31298-001_geneext_cellbender_output.pdf ./
scp kxw755@pegasus2.ccs.miami.edu:/scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/31298-001_geneext_cellbender_output_filtered.h5 ./
scp kxw755@pegasus2.ccs.miami.edu:/scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-001_geneext/31298-001_geneext_cellbender_output_metrics.csv ./
```

```bash
scp kxw755@pegasus2.ccs.miami.edu:/scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/31298-002_geneext_cellbender_output.pdf ./
scp kxw755@pegasus2.ccs.miami.edu:/scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/31298-002_geneext_cellbender_output_filtered.h5 ./
scp kxw755@pegasus2.ccs.miami.edu:/scratch/projects/dark_genes/GeneExt/Pdam_31298/31298-002_geneext/31298-002_geneext_cellbender_output_metrics.csv ./
```