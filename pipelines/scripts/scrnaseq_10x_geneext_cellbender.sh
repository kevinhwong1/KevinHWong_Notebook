#!/bin/bash
# =============================================================================
# 10x scRNA-seq: Cell Ranger -> GeneExt -> Cell Ranger -> CellBender  (LSF clusters, e.g. Pegasus)
# Pipeline page: https://kevinhwong1.github.io/KevinHWong_Notebook/pipelines/scrnaseq-10x-geneext-cellbender.html
#
# Steps (each is an LSF job that waits on the one before it):
#   0.  cellranger mkref   original GTF            (optional: RUN_MKREF_ORIG=true)
#   1.  cellranger count   every library, original ref   (BAMs for GeneExt)
#   2.  samtools merge     all library BAMs -> merge.bam
#   3.  GeneExt            extend 3' ends of gene models using merge.bam
#   4.  cellranger mkref   GeneExt GTF
#   5.  cellranger count   every library, GeneExt ref
#   6.  CellBender         ambient RNA removal on each library
#
# Usage:
#   bash scrnaseq_10x_geneext_cellbender.sh --init my_run.config    # 1. write a config template
#   nano my_run.config                                              # 2. replace every placeholder
#   bash scrnaseq_10x_geneext_cellbender.sh my_run.config           # 3. submit the whole chain
#   START_AT=3 bash scrnaseq_10x_geneext_cellbender.sh my_run.config  # resume from a step (0-6)
#
# Nothing in this file needs editing; all run-specific settings live in the config file.
# =============================================================================
set -euo pipefail

usage() { sed -n '15,19p' "$0" | sed 's/^# \{0,1\}//'; }

write_template() {
cat > "$1" <<'TEMPLATE'
# Config for scrnaseq_10x_geneext_cellbender.sh
# Replace every "<fill in: ...>" value. The script refuses to run while any are left.
# This file is sourced by bash: no spaces around "=", and quote values.

# --- Run and cluster ----------------------------------------------------------
RUN="<fill in: short run name, used in job and reference names, e.g. species_runID>"
LSF_PROJECT="<fill in: LSF project/allocation passed to bsub -P>"
EMAIL="<fill in: email for LSF job notifications>"
QUEUE="bigmem"          # LSF queue
CORES=16                # cores per job (bsub -n)
MEM_MB=15000            # bsub -R "rusage[mem=...]"

# --- Genome -------------------------------------------------------------------
FASTA="<fill in: /path/to/genome.fa>"
GTF="<fill in: /path/to/genes.gtf (exons only, consistent gene_id/transcript_id)>"
REF_ORIG="<fill in: /path/to/existing/cellranger/reference (or where step 0 should build it)>"
RUN_MKREF_ORIG=false    # true = build REF_ORIG from FASTA + GTF first (step 0)

# --- Libraries: one entry per library, same order in all three lists -----------
LIB_IDS=(     "<fill in: library ID 1, e.g. lib-001>"     "<fill in: library ID 2>" )
LIB_SAMPLES=( "<fill in: FASTQ sample prefix 1 (cellranger --sample)>" "<fill in: FASTQ sample prefix 2>" )
LIB_FASTQS=(  "<fill in: /path/to/fastq/folder/1>"         "<fill in: /path/to/fastq/folder/2>" )

# --- Output and software --------------------------------------------------------
OUTDIR="<fill in: /path/to/scratch/output/folder>"
GENEEXT_INSTALL="<fill in: /path/to/GeneExt (git clone containing geneext.py)>"

CELLRANGER_MODULE="cellranger/8.0.0"
CONDA_SH="${HOME}/anaconda3/bin/activate"   # conda activation script
SAMTOOLS_ENV="samtools_env"                 # conda env names
GENEEXT_ENV="geneext"
CELLBENDER_ENV="CellBender"

# --- Tool settings ----------------------------------------------------------------
GENEEXT_ARGS="-v 3 -j 8 --clip_strand both"        # --clip_strand both: no same-strand overlaps
CELLBENDER_ARGS="--fpr 0.01 --epochs 150 --learning-rate 0.00005"
TEMPLATE
}

# ----------------------------------------------------------------------------- read config
if [[ "${1:-}" == "--init" ]]; then
  cfg="${2:-run.config}"
  if [[ -e "${cfg}" ]]; then echo "${cfg} already exists; not overwriting." >&2; exit 1; fi
  write_template "${cfg}"
  echo "Wrote ${cfg}. Fill in every <fill in: ...> value, then run:"
  echo "  bash $0 ${cfg}"
  exit 0
fi

CONFIG="${1:-}"
if [[ -z "${CONFIG}" || ! -f "${CONFIG}" ]]; then usage; exit 1; fi

if grep -nE '^[^#]*<fill in' "${CONFIG}"; then
  echo "" >&2
  echo "Replace the placeholder values above in ${CONFIG} first." >&2
  exit 1
fi

# shellcheck source=/dev/null
source "${CONFIG}"
START_AT="${START_AT:-0}"

# ----------------------------------------------------------------------------- sanity checks
problems=0
need_file() { [[ -e "$1" ]] || { echo "Not found: $1  ($2)" >&2; problems=$((problems+1)); }; }
if [[ ${#LIB_IDS[@]} -ne ${#LIB_SAMPLES[@]} || ${#LIB_IDS[@]} -ne ${#LIB_FASTQS[@]} ]]; then
  echo "LIB_IDS, LIB_SAMPLES and LIB_FASTQS must have the same number of entries." >&2; problems=$((problems+1))
fi
need_file "${FASTA}" FASTA
need_file "${GTF}" GTF
need_file "${GENEEXT_INSTALL}/geneext.py" "GeneExt install"
if (( START_AT <= 1 )) && [[ "${RUN_MKREF_ORIG}" != true ]]; then need_file "${REF_ORIG}" "REF_ORIG (or set RUN_MKREF_ORIG=true)"; fi
if (( START_AT <= 5 )); then
  for d in "${LIB_FASTQS[@]}"; do need_file "${d}" "FASTQ folder"; done
fi
if (( problems > 0 )); then echo "Fix the ${problems} problem(s) above in ${CONFIG}." >&2; exit 1; fi

GENEEXT_DIR="${OUTDIR}/GeneExt"
GENEEXT_GTF="${OUTDIR}/${RUN}_GeneExt.gtf"
REF_GENEEXT="${OUTDIR}/${RUN}_geneext"
LOGDIR="${OUTDIR}/logs"
mkdir -p "${OUTDIR}" "${GENEEXT_DIR}" "${LOGDIR}"

# deps JID [JID ...] -> 'done(1) && done(2)', skipping empty IDs (steps not run this time)
deps() {
  local out="" j
  for j in "$@"; do
    [[ -z "${j}" ]] && continue
    out+="${out:+ && }done(${j})"
  done
  echo "${out}"
}

# submit NAME WALLTIME "DEPS" <<EOF ... EOF   -> prints the LSF job ID
# The job file is kept in ${LOGDIR} (not /tmp) so every run has a record of what was submitted.
submit() {
  local name=$1 walltime=$2 wait=$3 jobfile="${LOGDIR}/${1}.lsf"
  {
    echo "#!/bin/bash"
    echo "#BSUB -J ${RUN}_${name}"
    echo "#BSUB -q ${QUEUE}"
    echo "#BSUB -P ${LSF_PROJECT}"
    echo "#BSUB -n ${CORES}"
    echo "#BSUB -W ${walltime}"
    echo "#BSUB -R \"rusage[mem=${MEM_MB}]\""
    [[ -n "${wait}" ]] && echo "#BSUB -w \"${wait}\""
    echo "#BSUB -u ${EMAIL}"
    echo "#BSUB -o ${LOGDIR}/${name}_%J.out"
    echo "#BSUB -e ${LOGDIR}/${name}_%J.err"
    echo "#BSUB -N"
    echo "set -eo pipefail   # no -u: conda activate scripts reference unset variables"
    cat                                   # job body from the heredoc
  } > "${jobfile}"
  bsub < "${jobfile}" | grep -oP '(?<=Job <)\d+'
}

echo "== ${RUN}: submitting LSF job chain from step ${START_AT} ($(date)) =="

# ----------------------------------------------------------------------------- STEP 0
JID_MKREF_ORIG=""
if (( START_AT <= 0 )) && [[ "${RUN_MKREF_ORIG}" == true ]]; then
  JID_MKREF_ORIG=$(submit step0_mkref_orig 24:00 "" <<EOF
module load ${CELLRANGER_MODULE}
cd $(dirname "${REF_ORIG}")
cellranger mkref --genome=$(basename "${REF_ORIG}") --fasta=${FASTA} --genes=${GTF}
EOF
)
  echo "step 0  mkref (original)      ${JID_MKREF_ORIG}"
fi

# ----------------------------------------------------------------------------- STEP 1
JIDS_COUNT_ORIG=()
if (( START_AT <= 1 )); then
  for i in "${!LIB_IDS[@]}"; do
    jid=$(submit "step1_count_${LIB_IDS[$i]}" 120:00 "$(deps "${JID_MKREF_ORIG}")" <<EOF
module load ${CELLRANGER_MODULE}
cd ${OUTDIR}
cellranger count \\
  --id=${LIB_IDS[$i]} \\
  --transcriptome=${REF_ORIG} \\
  --fastqs=${LIB_FASTQS[$i]} \\
  --sample=${LIB_SAMPLES[$i]} \\
  --create-bam=true
EOF
)
    JIDS_COUNT_ORIG+=("${jid}")
    echo "step 1  count ${LIB_IDS[$i]} (original) ${jid}"
  done
fi

# ----------------------------------------------------------------------------- STEP 2
JID_MERGE=""
if (( START_AT <= 2 )); then
  BAMS=""
  for id in "${LIB_IDS[@]}"; do BAMS+=" ${OUTDIR}/${id}/outs/possorted_genome_bam.bam"; done
  JID_MERGE=$(submit step2_bam_merge 24:00 "$(deps "${JIDS_COUNT_ORIG[@]:-}")" <<EOF
source ${CONDA_SH}
conda activate ${SAMTOOLS_ENV}
samtools --version | head -1
cd ${GENEEXT_DIR}
samtools merge -f -@ ${CORES} merge.bam ${BAMS}
samtools quickcheck merge.bam && echo "merge.bam OK"
ls -lh ${BAMS} merge.bam
EOF
)
  echo "step 2  samtools merge        ${JID_MERGE}"
fi

# ----------------------------------------------------------------------------- STEP 3
JID_GENEEXT=""
if (( START_AT <= 3 )); then
  JID_GENEEXT=$(submit step3_geneext 120:00 "$(deps "${JID_MERGE}")" <<EOF
source ${CONDA_SH}
conda activate ${GENEEXT_ENV}
# GeneExt writes tmp/ and calls its R scripts with paths relative to the install dir,
# so it has to run from there, and a stale tmp/ from an earlier run must be removed first.
cd ${GENEEXT_INSTALL}
rm -rf tmp
python geneext.py -g ${GTF} -b ${GENEEXT_DIR}/merge.bam -o ${GENEEXT_GTF} ${GENEEXT_ARGS}
tail -20 ${GENEEXT_GTF}.GeneExt.log
EOF
)
  echo "step 3  GeneExt               ${JID_GENEEXT}"
fi

# ----------------------------------------------------------------------------- STEP 4
JID_MKREF_GE=""
if (( START_AT <= 4 )); then
  JID_MKREF_GE=$(submit step4_mkref_geneext 120:00 "$(deps "${JID_GENEEXT}")" <<EOF
module load ${CELLRANGER_MODULE}
cd ${OUTDIR}
cellranger mkref --genome=$(basename "${REF_GENEEXT}") --fasta=${FASTA} --genes=${GENEEXT_GTF}
EOF
)
  echo "step 4  mkref (GeneExt)       ${JID_MKREF_GE}"
fi

# ----------------------------------------------------------------------------- STEPS 5 + 6
for i in "${!LIB_IDS[@]}"; do
  id="${LIB_IDS[$i]}"
  JID_COUNT_GE=""
  if (( START_AT <= 5 )); then
    JID_COUNT_GE=$(submit "step5_count_${id}_geneext" 120:00 "$(deps "${JID_MKREF_GE}")" <<EOF
module load ${CELLRANGER_MODULE}
cd ${OUTDIR}
cellranger count \\
  --id=${id}_geneext \\
  --transcriptome=${REF_GENEEXT} \\
  --fastqs=${LIB_FASTQS[$i]} \\
  --sample=${LIB_SAMPLES[$i]} \\
  --create-bam=true
EOF
)
    echo "step 5  count ${id} (GeneExt)  ${JID_COUNT_GE}"
  fi

  JID_CB=$(submit "step6_cellbender_${id}" 120:00 "$(deps "${JID_COUNT_GE}")" <<EOF
source ${CONDA_SH}
conda activate ${CELLBENDER_ENV}
cd ${OUTDIR}/${id}_geneext
cellbender remove-background \\
  --input  ${OUTDIR}/${id}_geneext/outs/raw_feature_bc_matrix.h5 \\
  --output ${OUTDIR}/${id}_geneext/${id}_geneext_cellbender_output.h5 \\
  ${CELLBENDER_ARGS}
EOF
)
  echo "step 6  CellBender ${id}      ${JID_CB}"
done

cat <<EOF

All jobs submitted. Monitor with:  bjobs -w | grep ${RUN}
Logs and the exact job scripts:     ${LOGDIR}/
Outputs to copy back, per library:  ${OUTDIR}/<lib>_geneext/<lib>_geneext_cellbender_output{_filtered.h5,_metrics.csv,_report.html,.pdf}
                                    ${OUTDIR}/<lib>_geneext/outs/web_summary.html
EOF
