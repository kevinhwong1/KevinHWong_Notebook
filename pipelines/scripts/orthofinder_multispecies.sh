#!/bin/bash
# =============================================================================
# Multi-species orthology with OrthoFinder (LSF clusters, e.g. Pegasus)
# Pipeline page: https://kevinhwong1.github.io/KevinHWong_Notebook/pipelines/orthofinder-multispecies.html
#
# What it does:
#   1. Prepares one clean protein FASTA per species (download if a URL is given; drop any text before
#      the first header and blank lines; shorten headers to the first word; check for duplicate IDs).
#      Each file is named <SPECIES>.fa, so species codes become the column names in Orthogroups.tsv.
#   2. Submits ONE OrthoFinder job with all species together (MSA gene trees by default).
#   3. When it finishes, copies the key tables to <RUN_DIR>/summary/ and counts ortholog pairs.
#
# Usage:
#   bash orthofinder_multispecies.sh --init my_orthofinder.config   # 1. write a config template
#   nano my_orthofinder.config                                      # 2. replace every placeholder
#   bash orthofinder_multispecies.sh my_orthofinder.config          # 3. prepare FASTAs + submit
#   PREP_ONLY=1 bash orthofinder_multispecies.sh my_orthofinder.config  # only prepare/check FASTAs
#
# Nothing in this file needs editing; all run-specific settings live in the config file.
# =============================================================================
set -euo pipefail

usage() { sed -n '13,17p' "$0" | sed 's/^# \{0,1\}//'; }

write_template() {
cat > "$1" <<'TEMPLATE'
# Config for orthofinder_multispecies.sh
# Replace every "<fill in: ...>" value. The script refuses to run while any are left.
# This file is sourced by bash: no spaces around "=", and quote values.

# --- Run and cluster ----------------------------------------------------------
RUN="<fill in: short run name, e.g. Pdam_Amil_Spis_Nvec>"
RUN_DIR="<fill in: /path/to/scratch/orthofinder/RUN_NAME>"
LSF_PROJECT="<fill in: LSF project/allocation passed to bsub -P>"
EMAIL="<fill in: email for LSF job notifications>"
QUEUE="bigmem"
CORES=16
MEM_MB=15000
WALLTIME="120:00"

# --- Species: one entry per species, same order in both lists ----------------------
# SPECIES codes become the column names in Orthogroups.tsv (short, no spaces, e.g. Pdam).
# PROTEOMES can be a local path or an https:// URL (downloaded on the login node).
# Include an outgroup (e.g. Nematostella vectensis) to help OrthoFinder root gene trees.
SPECIES=(   "<fill in: code 1, e.g. Pdam>"   "<fill in: code 2>"   "<fill in: code 3>" )
PROTEOMES=( "<fill in: /path/or/URL/to/proteins_1.fasta>" "<fill in: proteins_2>" "<fill in: proteins_3>" )

# --- Software ---------------------------------------------------------------------
CONDA_SH="${HOME}/anaconda3/etc/profile.d/conda.sh"
ORTHOFINDER_ENV="orthofinder_env"     # OrthoFinder + MAFFT + FastTree; keep separate from other envs

# MSA gene trees (more accurate orthologs than the default DendroBLAST trees, but slower)
ORTHOFINDER_ARGS="-M msa -A mafft -T fasttree"
TEMPLATE
}

# ----------------------------------------------------------------------------- read config
if [[ "${1:-}" == "--init" ]]; then
  cfg="${2:-orthofinder.config}"
  if [[ -e "${cfg}" ]]; then echo "${cfg} already exists; not overwriting." >&2; exit 1; fi
  write_template "${cfg}"
  echo "Wrote ${cfg}. Fill in every <fill in: ...> value, then run:"
  echo "  bash $0 ${cfg}"
  exit 0
fi

CONFIG="${1:-}"
if [[ -z "${CONFIG}" || ! -f "${CONFIG}" ]]; then usage; exit 1; fi
if grep -nE '^[^#]*<fill in' "${CONFIG}"; then
  echo "" >&2; echo "Replace the placeholder values above in ${CONFIG} first." >&2; exit 1
fi
# shellcheck source=/dev/null
source "${CONFIG}"

# ----------------------------------------------------------------------------- checks
if [[ ${#SPECIES[@]} -ne ${#PROTEOMES[@]} ]]; then
  echo "SPECIES and PROTEOMES must have the same number of entries." >&2; exit 1
fi
if [[ ${#SPECIES[@]} -lt 2 ]]; then echo "Give at least 2 species." >&2; exit 1; fi
if [[ $(printf '%s\n' "${SPECIES[@]}" | sort | uniq -d | wc -l) -gt 0 ]]; then
  echo "SPECIES codes must be unique." >&2; exit 1
fi
for s in "${SPECIES[@]}"; do
  [[ "${s}" =~ ^[A-Za-z0-9_]+$ ]] || { echo "Species code '${s}': use only letters, numbers and _." >&2; exit 1; }
done

FASTA_DIR="${RUN_DIR}/proteomes"
RESULTS_PARENT="${RUN_DIR}/results"
RESULTS_DIR="${RESULTS_PARENT}/Results_${RUN}"
LOGDIR="${RUN_DIR}/logs"
SUMMARY_DIR="${RUN_DIR}/summary"
mkdir -p "${FASTA_DIR}" "${LOGDIR}"

if [[ -e "${RESULTS_DIR}" ]]; then
  echo "${RESULTS_DIR} already exists. Use a new RUN name (or move the old results)." >&2; exit 1
fi

# ----------------------------------------------------------------------------- 1. prepare FASTAs
echo "== Preparing proteomes in ${FASTA_DIR} =="
printf "species\tsequences\tsource\n" > "${RUN_DIR}/proteome_summary.tsv"
problems=0
for i in "${!SPECIES[@]}"; do
  sp="${SPECIES[$i]}"; src="${PROTEOMES[$i]}"; out="${FASTA_DIR}/${sp}.fa"
  tmp="${FASTA_DIR}/.${sp}.download"
  if [[ "${src}" =~ ^https?:// ]]; then
    echo "  ${sp}: downloading ${src}"
    curl -fsSL -o "${tmp}" "${src}" || { echo "  ${sp}: download failed" >&2; problems=$((problems+1)); continue; }
    in="${tmp}"
  else
    [[ -f "${src}" ]] || { echo "  ${sp}: not found: ${src}" >&2; problems=$((problems+1)); continue; }
    in="${src}"
  fi
  # unzip if needed; drop text before the first '>' and blank lines; keep the first word of each header
  if [[ "${in}" == *.gz ]] || gzip -t "${in}" 2>/dev/null; then cat_cmd="gzip -dc"; else cat_cmd="cat"; fi
  ${cat_cmd} "${in}" | tr -d '\r' | awk '
    /^>/ { seen=1; split(substr($0,2), a, /[ \t]/); print ">" a[1]; next }
    seen && NF { print }' > "${out}"
  rm -f "${tmp}"

  n=$(grep -c '^>' "${out}" || true)
  dups=$(grep '^>' "${out}" | sort | uniq -d | wc -l)
  if (( n == 0 )); then echo "  ${sp}: no sequences found!" >&2; problems=$((problems+1)); fi
  if (( dups > 0 )); then echo "  ${sp}: ${dups} duplicate protein IDs after shortening headers" >&2; problems=$((problems+1)); fi
  printf "%s\t%s\t%s\n" "${sp}" "${n}" "${src}" >> "${RUN_DIR}/proteome_summary.tsv"
  echo "  ${sp}: ${n} proteins"
done
if (( problems > 0 )); then echo "Fix the ${problems} problem(s) above before running OrthoFinder." >&2; exit 1; fi

echo "Check that each species has one protein per gene (longest isoform). If a proteome has many"
echo "isoforms per gene, filter it first; OrthoFinder treats every isoform as a separate gene."

if [[ "${PREP_ONLY:-0}" == 1 ]]; then echo "PREP_ONLY=1: stopping before submission."; exit 0; fi

# ----------------------------------------------------------------------------- 2. submit OrthoFinder
JOBFILE="${LOGDIR}/orthofinder_${RUN}.lsf"
cat > "${JOBFILE}" <<EOF
#!/bin/bash
#BSUB -J OF_${RUN}
#BSUB -q ${QUEUE}
#BSUB -P ${LSF_PROJECT}
#BSUB -n ${CORES}
#BSUB -W ${WALLTIME}
#BSUB -R "rusage[mem=${MEM_MB}]"
#BSUB -u ${EMAIL}
#BSUB -o ${LOGDIR}/orthofinder_${RUN}_%J.out
#BSUB -e ${LOGDIR}/orthofinder_${RUN}_%J.err
#BSUB -B
#BSUB -N
set -eo pipefail   # no -u: conda activate scripts reference unset variables

source ${CONDA_SH}
conda activate ${ORTHOFINDER_ENV}
orthofinder --help | head -2 || true

orthofinder -f ${FASTA_DIR} -o ${RESULTS_PARENT} -n ${RUN} -t ${CORES} -a ${CORES} ${ORTHOFINDER_ARGS}

# ---- summary tables
mkdir -p ${SUMMARY_DIR}
cp ${RESULTS_DIR}/Orthogroups/Orthogroups.tsv \\
   ${RESULTS_DIR}/Orthogroups/Orthogroups.GeneCount.tsv \\
   ${RESULTS_DIR}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt \\
   ${RESULTS_DIR}/Comparative_Genomics_Statistics/Statistics_Overall.tsv \\
   ${SUMMARY_DIR}/

# ortholog pairs per species pair: total rows and strict 1:1 rows (one gene on each side)
printf "species_a\tspecies_b\tortholog_rows\tone_to_one\n" > ${SUMMARY_DIR}/ortholog_pair_summary.tsv
for f in ${RESULTS_DIR}/Orthologues/Orthologues_*/*__v__*.tsv; do
  pair=\$(basename "\$f" .tsv)
  awk -F'\t' -v a="\${pair%%__v__*}" -v b="\${pair##*__v__}" '
    NR>1 { rows++; if (\$2 !~ /,/ && \$3 !~ /,/) one++ }
    END  { printf "%s\t%s\t%d\t%d\n", a, b, rows, one }' "\$f" >> ${SUMMARY_DIR}/ortholog_pair_summary.tsv
done
echo "Done. Key tables in ${SUMMARY_DIR}/; full results in ${RESULTS_DIR}/"
EOF

JID=$(bsub < "${JOBFILE}" | grep -oP '(?<=Job <)\d+')
cat <<EOF

Submitted OrthoFinder job ${JID} for ${#SPECIES[@]} species (${SPECIES[*]}).
Monitor with:     bjobs ${JID}      (or bpeek ${JID})
Job script:       ${JOBFILE}
Results will be:  ${RESULTS_DIR}/
Key tables:       ${SUMMARY_DIR}/  (Orthogroups.tsv, ortholog_pair_summary.tsv, ...)
EOF
