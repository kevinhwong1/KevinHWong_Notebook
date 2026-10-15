#!/usr/bin/env bash
# ============================================================
# 02_embed.sh  —  ESM protein embedding for one species
#
# Usage (submit directly with bsub):
#   SPECIES=Mlei  bsub [flags] < scripts/02_embed.sh
#   SPECIES=Cgig  bsub [flags] < scripts/02_embed.sh
#   SPECIES=Spur  bsub [flags] < scripts/02_embed.sh
#   SPECIES=Crob  bsub [flags] < scripts/02_embed.sh
#   SPECIES=Drer  bsub [flags] < scripts/02_embed.sh
#   SPECIES=Dmel  bsub [flags] < scripts/02_embed.sh
#   SPECIES=Hsap  bsub [flags] < scripts/02_embed.sh
# ============================================================

#BSUB -J embed_SPECIES
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=20000]"
#BSUB -o logs/embed_SPECIES_%J.out
#BSUB -e logs/embed_SPECIES_%J.err
#BSUB -u kxw755@earth.miami.edu
#BSUB -N

set -euo pipefail

SCRIPT_DIR="/scratch/dark_genes/SATURN_Mnemi/scripts"
source "${SCRIPT_DIR}/run_config.sh"

: "${SPECIES:?ERROR: SPECIES not set. Run as: SPECIES=Mlei bsub < 02_embed.sh}"

get_species_config "${SPECIES}"

# conda activation scripts reference unset vars in batch jobs: relax -u here (as in 05a)
set +u
source ~/anaconda3/bin/activate
conda activate "${CONDA_ENV}"
set -u

echo "=== 02_embed.sh  SPECIES=${SPECIES}  STRATEGY=${EMB_STRATEGY} ==="
date
hostname

mkdir -p "${EMB_DIR}" "${EMB_PT_DIR}"

# Write all Python helper scripts to the embeddings dir
# This avoids heredoc expansion issues entirely

PY_DIR="${EMB_DIR}/py_scripts"
mkdir -p "${PY_DIR}"

# ============================================================
# Write step0.py — export gene list
# ============================================================
cat > "${PY_DIR}/step0.py" << 'PYEOF'
import sys, scanpy as sc
from pathlib import Path
h5ad  = sys.argv[1]
genes_txt = sys.argv[2]
ad = sc.read_h5ad(h5ad, backed="r")
genes = sorted(set(ad.var_names.astype(str)))
Path(genes_txt).write_text("\n".join(genes) + "\n")
print("Genes in h5ad   : " + str(ad.n_vars))
print("Unique var_names: " + str(len(genes)))
print("First 5         : " + str(genes[:5]))
PYEOF

# ============================================================
# Write step1_direct.py
# ============================================================
cat > "${PY_DIR}/step1_direct.py" << 'PYEOF'
import sys, json
from pathlib import Path
genes_txt, fasta_in, out_json, out_fasta = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
genes = set(x.strip() for x in open(genes_txt) if x.strip())
fasta_ids = set()
with open(fasta_in) as f:
    for line in f:
        if line.startswith(">"): fasta_ids.add(line[1:].split()[0].strip())
gene_to_prot = {g: [g] for g in genes if g in fasta_ids}
print("Genes in h5ad    : " + str(len(genes)))
print("Proteins in FASTA: " + str(len(fasta_ids)))
print("Mapped (direct)  : " + str(len(gene_to_prot)))
print("Unmapped         : " + str(len(genes) - len(gene_to_prot)))
with open(out_json, "w") as f: json.dump(gene_to_prot, f, indent=2)
prot_keep = set(p for ps in gene_to_prot.values() for p in ps)
seen = set(); kept = 0
with open(fasta_in) as fin, open(out_fasta, "w") as fout:
    write = False
    for line in fin:
        if line.startswith(">"):
            pid = line[1:].split()[0].strip()
            if pid in prot_keep and pid not in seen:
                write = True; seen.add(pid)
                fout.write(">" + pid + "\n"); kept += 1
            else: write = False
        elif write: fout.write(line)
print("Filtered FASTA   : " + str(kept) + " proteins")
PYEOF

# ============================================================
# Write step1_gtf_loc_to_xp.py
# ============================================================
cat > "${PY_DIR}/step1_gtf_loc_to_xp.py" << 'PYEOF'
import sys, re, json, csv
from collections import defaultdict
from pathlib import Path
genes_txt, gtf_in, fasta_in, out_tsv, out_json, out_fasta = \
    sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]
genes = set(x.strip() for x in open(genes_txt) if x.strip())
print("Genes in h5ad : " + str(len(genes)))
loc_to_xp = defaultdict(set)
n_cds = 0
for line in open(gtf_in):
    if line.startswith("#"): continue
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 9 or parts[2] != "CDS": continue
    attr = parts[8]
    m_loc = re.search(r'gene_id "([^"]+)"', attr)
    m_xp  = re.search(r'protein_id "([^"]+)"', attr)
    if m_loc and m_xp:
        loc, xp = m_loc.group(1), m_xp.group(1)
        if loc in genes: loc_to_xp[loc].add(xp); n_cds += 1
print("GTF CDS matched : " + str(n_cds))
print("LOC->XP mapped  : " + str(len(loc_to_xp)) + " of " + str(len(genes)))
print("Unmapped LOCs   : " + str(len(genes) - len(loc_to_xp)))
with open(out_tsv, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["loc_id", "protein_id"])
    for loc, xps in sorted(loc_to_xp.items()):
        for xp in sorted(xps): w.writerow([loc, xp])
gene_to_prot = {loc: sorted(xps) for loc, xps in loc_to_xp.items()}
with open(out_json, "w") as f: json.dump(gene_to_prot, f, indent=2)
xp_needed = set(xp for xps in loc_to_xp.values() for xp in xps)
print("XP IDs needed   : " + str(len(xp_needed)))
seq_buf = {}; cur_pid = None; cur_lines = []
with open(fasta_in) as f:
    for line in f:
        if line.startswith(">"):
            if cur_pid and cur_pid in xp_needed: seq_buf[cur_pid] = "".join(cur_lines)
            cur_pid = line[1:].split()[0].strip(); cur_lines = [line]
        else: cur_lines.append(line)
    if cur_pid and cur_pid in xp_needed: seq_buf[cur_pid] = "".join(cur_lines)
with open(out_fasta, "w") as f:
    for pid, block in seq_buf.items():
        lines = block.rstrip("\n").split("\n")
        f.write(">" + pid + "\n")
        for l in lines[1:]: f.write(l + "\n")
print("Filtered FASTA  : " + str(len(seq_buf)) + " proteins")
PYEOF

# ============================================================
# Write step1_direct_strip_suffix.py
# ============================================================
cat > "${PY_DIR}/step1_direct_strip_suffix.py" << 'PYEOF'
import sys, re, json
genes_txt, fasta_in, out_json, out_fasta = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
genes = set(x.strip() for x in open(genes_txt) if x.strip())
print("Genes in h5ad    : " + str(len(genes)))
strip_re = re.compile(r'^(KY21\.[^.]+\.\d+)')
stripped_to_full = {}
with open(fasta_in) as f:
    for line in f:
        if line.startswith(">"):
            full_id = line[1:].split()[0].strip()
            m = strip_re.match(full_id)
            if m: stripped_to_full.setdefault(m.group(1), []).append(full_id)
gene_to_prot = {g: stripped_to_full[g] for g in genes if g in stripped_to_full}
print("Unique base IDs  : " + str(len(stripped_to_full)))
print("Mapped           : " + str(len(gene_to_prot)))
print("Unmapped         : " + str(len(genes) - len(gene_to_prot)))
with open(out_json, "w") as f: json.dump(gene_to_prot, f, indent=2)
prot_keep = set(p for ps in gene_to_prot.values() for p in ps)
seen = set(); kept = 0
with open(fasta_in) as fin, open(out_fasta, "w") as fout:
    write = False
    for line in fin:
        if line.startswith(">"):
            pid = line[1:].split()[0].strip()
            if pid in prot_keep and pid not in seen:
                write = True; seen.add(pid)
                fout.write(">" + pid + "\n"); kept += 1
            else: write = False
        elif write: fout.write(line)
print("Filtered FASTA   : " + str(kept) + " proteins")
PYEOF

# ============================================================
# Write step1_gtf_symbol_to_ensdarg.py
# ============================================================
cat > "${PY_DIR}/step1_gtf_symbol_to_ensdarg.py" << 'PYEOF'
import sys, re, json, csv
from collections import defaultdict
genes_txt, gtf_in, fasta_in, out_tsv, out_json, out_fasta = \
    sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]
genes = set(x.strip() for x in open(genes_txt) if x.strip())
print("Genes (symbols)  : " + str(len(genes)))
sym_to_ensdarg = defaultdict(set)
n_gene = 0
for line in open(gtf_in):
    if line.startswith("#"): continue
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 9 or parts[2] != "gene": continue
    attr = parts[8]
    m_id  = re.search(r'gene_id "([^"]+)"', attr)
    m_sym = re.search(r'gene_name "([^"]+)"', attr)
    if m_id and m_sym:
        sym, gid = m_sym.group(1), m_id.group(1)
        n_gene += 1
        if sym in genes: sym_to_ensdarg[sym].add(gid)
ensdarg_needed = set(g for gs in sym_to_ensdarg.values() for g in gs)
print("GTF gene lines   : " + str(n_gene))
print("Symbols->ENSDARG : " + str(len(sym_to_ensdarg)))
gene_re = re.compile(r'gene:(ENSDARG\d+)')
ensdarg_to_ensdarp = defaultdict(list)
with open(fasta_in) as f:
    for line in f:
        if not line.startswith(">"): continue
        header = line[1:].strip()
        ensdarp = header.split()[0]
        m = gene_re.search(header)
        if m:
            ensdarg = m.group(1)
            if ensdarg in ensdarg_needed: ensdarg_to_ensdarp[ensdarg].append(ensdarp)
print("ENSDARG->ENSDARP : " + str(len(ensdarg_to_ensdarp)))
gene_to_prot = {}
for sym, ensdargs in sym_to_ensdarg.items():
    endarps = []
    for ensdarg in ensdargs: endarps.extend(ensdarg_to_ensdarp.get(ensdarg, []))
    if endarps: gene_to_prot[sym] = sorted(set(endarps))
print("Symbols->ENSDARP : " + str(len(gene_to_prot)))
print("Unmapped         : " + str(len(genes) - len(gene_to_prot)))
with open(out_tsv, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene_symbol", "ensdarp_id"])
    for sym, endarps in sorted(gene_to_prot.items()):
        for ep in endarps: w.writerow([sym, ep])
with open(out_json, "w") as jf: json.dump(gene_to_prot, jf, indent=2)
prot_keep = set(p for ps in gene_to_prot.values() for p in ps)
seen = set(); kept = 0
with open(fasta_in) as fin, open(out_fasta, "w") as fout:
    write = False
    for line in fin:
        if line.startswith(">"):
            pid = line[1:].split()[0].strip()
            if pid in prot_keep and pid not in seen:
                write = True; seen.add(pid)
                fout.write(">" + pid + "\n"); kept += 1
            else: write = False
        elif write: fout.write(line)
print("Filtered FASTA   : " + str(kept) + " proteins")
PYEOF

# ============================================================
# Write step1_ensembl_gene_id.py  (Dmel, Hsap; added 2026-10-15)
# var_names are Ensembl gene IDs (FBgn / ENSG); map each gene to ALL its
# proteins via the "gene:" field of Ensembl pep.all headers (version suffix
# stripped), i.e. the same all-isoforms mapping used for Drer.
# ============================================================
cat > "${PY_DIR}/step1_ensembl_gene_id.py" << 'PYEOF'
import sys, re, json, csv
from collections import defaultdict
genes_txt, fasta_in, out_tsv, out_json, out_fasta = \
    sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
genes = set(x.strip() for x in open(genes_txt) if x.strip())
print("Genes in h5ad    : " + str(len(genes)))
gene_re = re.compile(r'gene:([^ .]+)')
gene_to_prot = defaultdict(list)
n_prot = 0
with open(fasta_in) as f:
    for line in f:
        if not line.startswith(">"): continue
        n_prot += 1
        header = line[1:].strip()
        pid = header.split()[0]
        m = gene_re.search(header)
        if m and m.group(1) in genes: gene_to_prot[m.group(1)].append(pid)
gene_to_prot = {g: sorted(set(ps)) for g, ps in gene_to_prot.items()}
n_map = sum(len(v) for v in gene_to_prot.values())
print("Proteins in FASTA: " + str(n_prot))
print("Genes mapped     : " + str(len(gene_to_prot)))
print("Proteins mapped  : " + str(n_map))
print("Unmapped genes   : " + str(len(genes) - len(gene_to_prot)) + "  (mostly non-coding)")
with open(out_tsv, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene_id", "protein_id"])
    for g, ps in sorted(gene_to_prot.items()):
        for pid in ps: w.writerow([g, pid])
with open(out_json, "w") as jf: json.dump(gene_to_prot, jf, indent=2)
prot_keep = set(p for ps in gene_to_prot.values() for p in ps)
seen = set(); kept = 0
with open(fasta_in) as fin, open(out_fasta, "w") as fout:
    write = False
    for line in fin:
        if line.startswith(">"):
            pid = line[1:].split()[0].strip()
            if pid in prot_keep and pid not in seen:
                write = True; seen.add(pid)
                fout.write(">" + pid + "\n"); kept += 1
            else: write = False
        elif write: fout.write(line)
print("Filtered FASTA   : " + str(kept) + " proteins")
PYEOF

# ============================================================
# Write step5_validate.py
# ============================================================
cat > "${PY_DIR}/step5_validate.py" << 'PYEOF'
import sys, torch
from pathlib import Path
pt_path, genes_txt = sys.argv[1], sys.argv[2]
genes_in = set(x.strip() for x in open(genes_txt) if x.strip())
pt = Path(pt_path)
if not pt.exists():
    print("WARNING: embedding file not found: " + pt_path); sys.exit(1)
data = torch.load(str(pt))
gene_ids = list(data.keys())
emb_dim  = next(iter(data.values())).shape[-1]
overlap  = len(set(gene_ids) & genes_in)
pct = 100.0 * overlap / len(genes_in) if genes_in else 0
print("Embedding file   : " + pt.name)
print("Genes embedded   : " + str(len(gene_ids)))
print("Genes in h5ad    : " + str(len(genes_in)))
print("Overlap          : " + str(overlap) + "  (" + str(round(pct,1)) + "%)")
print("Embedding dim    : " + str(emb_dim))
print("First 5 IDs      : " + str(gene_ids[:5]))
PYEOF

# ============================================================
# STEP 0: Export gene list
# ============================================================
echo ""
echo "--- Step 0: Export gene list ---"
python "${PY_DIR}/step0.py" "${H5AD_SATURN}" "${EMB_GENES_TXT}"

# ============================================================
# STEP 1: Build gene→protein map + filtered FASTA
# ============================================================
echo ""
echo "--- Step 1: Build gene->protein map (strategy: ${EMB_STRATEGY}) ---"

if [[ "${EMB_STRATEGY}" == "direct" ]]; then
    python "${PY_DIR}/step1_direct.py" \
        "${EMB_GENES_TXT}" "${PROTEOME_FASTA}" "${EMB_MAP_JSON}" "${EMB_FILTERED_FASTA}"

elif [[ "${EMB_STRATEGY}" == "gtf_loc_to_xp" ]]; then
    python "${PY_DIR}/step1_gtf_loc_to_xp.py" \
        "${EMB_GENES_TXT}" "${GTF}" "${PROTEOME_FASTA}" \
        "${EMB_LOC_MAP_TSV}" "${EMB_MAP_JSON}" "${EMB_FILTERED_FASTA}"

elif [[ "${EMB_STRATEGY}" == "direct_strip_suffix" ]]; then
    python "${PY_DIR}/step1_direct_strip_suffix.py" \
        "${EMB_GENES_TXT}" "${PROTEOME_FASTA}" "${EMB_MAP_JSON}" "${EMB_FILTERED_FASTA}"

elif [[ "${EMB_STRATEGY}" == "gtf_symbol_to_ensdarg" ]]; then
    python "${PY_DIR}/step1_gtf_symbol_to_ensdarg.py" \
        "${EMB_GENES_TXT}" "${GTF}" "${PROTEOME_FASTA}" \
        "${EMB_LOC_MAP_TSV}" "${EMB_MAP_JSON}" "${EMB_FILTERED_FASTA}"

elif [[ "${EMB_STRATEGY}" == "ensembl_gene_id" ]]; then
    python "${PY_DIR}/step1_ensembl_gene_id.py" \
        "${EMB_GENES_TXT}" "${PROTEOME_FASTA}" \
        "${EMB_LOC_MAP_TSV}" "${EMB_MAP_JSON}" "${EMB_FILTERED_FASTA}"

else
    echo "ERROR: Unknown EMB_STRATEGY '${EMB_STRATEGY}'" >&2
    exit 1
fi

# ============================================================
# STEP 2: Run ESM extract.py
# ============================================================
echo ""
echo "--- Step 2: ESM extract.py ---"
N_PROT=$(grep -c '^>' "${EMB_FILTERED_FASTA}" || true)
echo "Proteins to embed: ${N_PROT}"

python "${ESM_EXTRACT}" \
    "${ESM_MODEL}" \
    "${EMB_FILTERED_FASTA}" \
    "${EMB_PT_DIR}" \
    --repr_layers 33 \
    --include mean \
    --truncation_seq_length 1022

echo "ESM complete. PT files: $(find ${EMB_PT_DIR} -name '*.pt' | wc -l)"

# ============================================================
# STEP 3: Filter map to proteins with .pt files
# ============================================================
echo ""
echo "--- Step 3: Filter map to embedded proteins ---"
python "${ESM_FILTER}" \
    --input_json    "${EMB_MAP_JSON}" \
    --embedding_dir "${EMB_PT_DIR}" \
    --output_json   "${EMB_FILTERED_JSON}"

# ============================================================
# STEP 4: Convert .pt files → gene-level embedding
# ============================================================
echo ""
echo "--- Step 4: Convert to gene embeddings ---"
python "${ESM_CONVERT}" \
    --embedding_dir "${EMB_PT_DIR}" \
    --gene_symbol_to_protein_ids_path "${EMB_FILTERED_JSON}" \
    --embedding_model "ESM1b" \
    --save_path "${EMB_PT}"

# ============================================================
# STEP 5: Validate
# ============================================================
echo ""
echo "--- Step 5: Validate ---"
python "${PY_DIR}/step5_validate.py" "${EMB_PT}" "${EMB_GENES_TXT}"

echo ""
echo "=== 02_embed.sh DONE  SPECIES=${SPECIES} ==="
date
