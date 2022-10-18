---
layout: post
title: Troubleshooting WGBS Scripts
date: '2022-10-18'
categories: Analysis
tags: WGBS, Porites astreoides
---

This this my previous WGBS scripts to archive. This analysis includes all samples, which I am not currently using for my analysis anymore.


visual# *Porites astreoides* Thermal Transplant WGBS Workflow
Script modified from [D Becker-Polinski](https://github.com/hputnam/Becker_E5/blob/master/Bioinformatics/Workflows/Becker_WGBS_Workflow.md)

## 1) Methylation Quantification with Methylseq for *P. astreoides*

A few testing parameters for the settings can be found in this [notebook post](https://kevinhwong1.github.io/KevinHWong_Notebook/Methylseq-trimming-test-to-remove-m-bias/). I decided to use the Trim_3 parameters as it keeps the most amount of data while still reducing the m-bias.

I also used the *Porites astreoides* genome from [this repository](https://github.com/hputnam/Past_Genome).

### methylseq_trim3.sh

```bash
#!/bin/bash
#SBATCH --job-name="methylseq"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3
#SBATCH --exclusive

# load modules needed

module load Nextflow/21.03.0

# run nextflow methylseq

nextflow run nf-core/methylseq \
-profile singularity \
--aligner bismark \
--igenomes_ignore \
--fasta /data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta \
--save_reference \
--input '/data/putnamlab/KITT/hputnam/20211008_Past_ThermalTransplant_WGBS/*_R{1,2}_001.fastq.gz' \
--clip_r1 15 \
--clip_r2 30 \
--three_prime_clip_r1 30 \
--three_prime_clip_r2 15 \
--non_directional \
--cytosine_report \
--relax_mismatches \
--unmapped \
--outdir WGBS_methylseq
```

```
scp kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/methylseq_trim3/WGBS_methylseq/MultiQC/multiqc_report.html MyProjects/Thermal_Transplant_Molecular/output/multiqc_report_trim3_full.html
```

## 2) Merge strands

The Bismark coverage2cytosine command re-reads the genome-wide report and merges methylation evidence of both top and bottom strand.

`mkdir cov_to_cyto`

### cov_to_cyto.sh

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto
#SBATCH --exclusive

# load modules needed

module load Bismark/0.20.1-foss-2018b

# run coverage2cytosine merge of strands

 find /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/methylseq_trim3/WGBS_methylseq/bismark_methylation_calls/methylation_coverage/*deduplicated.bismark.cov.gz \
 | xargs basename -s _L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz \
 | xargs -I{} coverage2cytosine \
 --genome_folder /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/methylseq_trim3/WGBS_methylseq/reference_genome/BismarkIndex \
 -o {} \
 --merge_CpG \
 --zero_based \
/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/methylseq_trim3/WGBS_methylseq/bismark_methylation_calls/methylation_coverage/{}_L004_R1_001_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
```

## 3) Sorting the merged files so scaffolds are all in the same order and multiIntersectBed will run correctly. Run for loop using bedtools to sort all .tab files

`nano bedtools.sort.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto
#SBATCH --cpus-per-task=3

module load BEDTools/2.27.1-foss-2018b

for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  bedtools sort -i "${f}" \
  > "${STEM}"_sorted.cov
done
```

### Create files for statistical analysis

#### Run loop to filter CpGs for 5x coverage, creating tab files with raw count for glms

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_5x_sorted.tab
done
```

`wc -l *5x_sorted.tab`

```
4527111 18-106_S163_5x_sorted.tab
5597804 18-118_S162_5x_sorted.tab
4104763 18-130_S172_5x_sorted.tab
4884343 18-142_S189_5x_sorted.tab
2463417 18-167_S166_5x_sorted.tab
5188393 18-178_S191_5x_sorted.tab
5235506 18-190_S186_5x_sorted.tab
3310274 18-202_S188_5x_sorted.tab
4670049 18-20_S202_5x_sorted.tab
4704858 18-227_S170_5x_sorted.tab
3954330 18-239_S185_5x_sorted.tab
6100768 18-250_S195_5x_sorted.tab
6854723 18-262_S179_5x_sorted.tab
5632004 18-311_S187_5x_sorted.tab
5462756 18-322_S180_5x_sorted.tab
5424011 18-32_S178_5x_sorted.tab
 377146 18-334_S164_5x_sorted.tab
7094939 18-346_S193_5x_sorted.tab
4928861 18-358_S201_5x_sorted.tab
5907039 18-370_S171_5x_sorted.tab
4421506 18-394_S192_5x_sorted.tab
5986501 18-406_S177_5x_sorted.tab
5640553 18-418_S196_5x_sorted.tab
7830166 18-442_S165_5x_sorted.tab
6501465 18-44_S198_5x_sorted.tab
7627601 18-454_S197_5x_sorted.tab
7787309 18-466_S199_5x_sorted.tab
5835415 18-55_S190_5x_sorted.tab
4083785 18-67_S176_5x_sorted.tab
7046855 18-79_S181_5x_sorted.tab
1035017 18-91_S160_5x_sorted.tab
 725233 18-9_S159_5x_sorted.tab
4867549 L-1029_S183_5x_sorted.tab
4555292 L-1038_S184_5x_sorted.tab
6352836 L-1053_S167_5x_sorted.tab
7551302 L-1059_S175_5x_sorted.tab
7108688 L-1093_S168_5x_sorted.tab
8781711 L-1257_S205_5x_sorted.tab
8085738 L-1263_S173_5x_sorted.tab
4471821 L-562_S174_5x_sorted.tab
7147304 L-571_S194_5x_sorted.tab
6156923 L-661_S182_5x_sorted.tab
6763146 L-704_S169_5x_sorted.tab
2165197 L-728_S161_5x_sorted.tab
4272213 L-862_S200_5x_sorted.tab
5074266 L-924_S204_5x_sorted.tab
6627315 L-933_S203_5x_sorted.tab
250925802 total
```

### Run loop to filter CpGs for 10x coverage, creating tab files with raw count for glms

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4, $5, $6}}' \
  > "${STEM}"_10x_sorted.tab
done
```

`wc -l *10x_sorted.tab`

```
676580 18-106_S163_10x_sorted.tab
1118168 18-118_S162_10x_sorted.tab
847700 18-130_S172_10x_sorted.tab
1142935 18-142_S189_10x_sorted.tab
210115 18-167_S166_10x_sorted.tab
1315924 18-178_S191_10x_sorted.tab
1055468 18-190_S186_10x_sorted.tab
411100 18-202_S188_10x_sorted.tab
825928 18-20_S202_10x_sorted.tab
1134939 18-227_S170_10x_sorted.tab
696483 18-239_S185_10x_sorted.tab
1787233 18-250_S195_10x_sorted.tab
1943403 18-262_S179_10x_sorted.tab
1543061 18-311_S187_10x_sorted.tab
1427802 18-322_S180_10x_sorted.tab
1068870 18-32_S178_10x_sorted.tab
  8494 18-334_S164_10x_sorted.tab
2079079 18-346_S193_10x_sorted.tab
1154658 18-358_S201_10x_sorted.tab
1390085 18-370_S171_10x_sorted.tab
865273 18-394_S192_10x_sorted.tab
1688362 18-406_S177_10x_sorted.tab
1464387 18-418_S196_10x_sorted.tab
2583003 18-442_S165_10x_sorted.tab
1917367 18-44_S198_10x_sorted.tab
2678858 18-454_S197_10x_sorted.tab
2749194 18-466_S199_10x_sorted.tab
1629133 18-55_S190_10x_sorted.tab
597189 18-67_S176_10x_sorted.tab
2694546 18-79_S181_10x_sorted.tab
 53547 18-91_S160_10x_sorted.tab
 27615 18-9_S159_10x_sorted.tab
976316 L-1029_S183_10x_sorted.tab
900254 L-1038_S184_10x_sorted.tab
1551870 L-1053_S167_10x_sorted.tab
3312304 L-1059_S175_10x_sorted.tab
2415616 L-1093_S168_10x_sorted.tab
3295825 L-1257_S205_10x_sorted.tab
3449982 L-1263_S173_10x_sorted.tab
1051976 L-562_S174_10x_sorted.tab
2590771 L-571_S194_10x_sorted.tab
1920044 L-661_S182_10x_sorted.tab
1790820 L-704_S169_10x_sorted.tab
205111 L-728_S161_10x_sorted.tab
642700 L-862_S200_10x_sorted.tab
1495979 L-924_S204_10x_sorted.tab
2332830 L-933_S203_10x_sorted.tab
68718897 total
```

## 4) Create a file with positions found in all samples at specified coverage

`nano 5x_intersect.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

multiIntersectBed -i *_5x_sorted.tab > CpG.all.samps.5x_sorted.bed

cat CpG.all.samps.5x_sorted.bed | awk '$4 ==47' > CpG.filt.all.samps.5x_sorted.bed

```

`nano 10x_intersect.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

multiIntersectBed -i *_10x_sorted.tab > CpG.all.samps.10x_sorted.bed

cat CpG.all.samps.10x_sorted.bed | awk '$4 ==47' > CpG.filt.all.samps.10x_sorted.bed

```

## 5) Create bedgraphs post merge

5x bedgraph

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4}}' \
  > "${STEM}"_5x_sorted.bedgraph
done
```

10 x bedgraph

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4}}' \
  > "${STEM}"_10x_sorted.bedgraph
done
```

## 6) Use intersectBed to find where loci and genes intersect, allowing loci to be mapped to annotated genes

#### Filtered gff to only include gene positions modified.gff > gene.gff

```bash
cd /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1
awk '{if ($3 == "gene") {print}}' Pastreoides_all_v1.gff  > Pastreoides_all_v1.gene.gff
```

`nano 5x_intersectBed.sh`

```
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

for i in *5x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gene.gff \
  > ${i}_gene
done
```

`nano 10x_intersectBed.sh`

```
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

for i in *10x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gene.gff \
  > ${i}_gene
done
```

## 7) Intersect with file to subset only those positions found in all samples

`nano 5x_intersect_final.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

for i in *_5x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.all.samps.5x_sorted.bed \
  > ${i}_CpG_5x_enrichment.bed
done
```

`wc -l *5x_enrichment.bed`

```
5231 18-106_S163_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-118_S162_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-130_S172_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-142_S189_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-167_S166_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-178_S191_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-190_S186_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-202_S188_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-20_S202_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-227_S170_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-239_S185_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-250_S195_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-262_S179_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-311_S187_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-322_S180_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-32_S178_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-334_S164_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-346_S193_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-358_S201_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-370_S171_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-394_S192_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-406_S177_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-418_S196_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-442_S165_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-44_S198_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-454_S197_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-466_S199_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-55_S190_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-67_S176_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-79_S181_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-91_S160_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 18-9_S159_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-1029_S183_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-1038_S184_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-1053_S167_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-1059_S175_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-1093_S168_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-1257_S205_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-1263_S173_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-562_S174_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-571_S194_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-661_S182_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-704_S169_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-728_S161_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-862_S200_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-924_S204_5x_sorted.tab_gene_CpG_5x_enrichment.bed
5231 L-933_S203_5x_sorted.tab_gene_CpG_5x_enrichment.bed
245857 total
```

`nano 10x_intersect_final.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

for i in *_10x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.all.samps.10x_sorted.bed \
  > ${i}_CpG_10x_enrichment.bed
done
```

`wc -l *10x_enrichment.bed`

```
467 18-106_S163_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-118_S162_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-130_S172_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-142_S189_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-167_S166_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-178_S191_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-190_S186_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-202_S188_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-20_S202_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-227_S170_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-239_S185_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-250_S195_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-262_S179_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-311_S187_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-322_S180_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-32_S178_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-334_S164_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-346_S193_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-358_S201_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-370_S171_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-394_S192_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-406_S177_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-418_S196_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-442_S165_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-44_S198_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-454_S197_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-466_S199_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-55_S190_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-67_S176_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-79_S181_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-91_S160_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 18-9_S159_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-1029_S183_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-1038_S184_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-1053_S167_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-1059_S175_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-1093_S168_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-1257_S205_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-1263_S173_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-562_S174_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-571_S194_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-661_S182_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-704_S169_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-728_S161_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-862_S200_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-924_S204_10x_sorted.tab_gene_CpG_10x_enrichment.bed
467 L-933_S203_10x_sorted.tab_gene_CpG_10x_enrichment.bed
```

## Notes

From here, we are losing alot of data at step 4. I am going to run another analysis where I remove samples 18-334, 18-91 and 18-9 as they have low number of reads compared to the other samples.

# Starting at step 4 with the removal of samples with low reads

```bash
mkdir cov_to_cyto_reduced

cd cov_to_cyto

cp *_sorted.cov ../cov_to_cyto_reduced
cp *_5x_sorted.tab ../cov_to_cyto_reduced
cp *_10x_sorted.tab ../cov_to_cyto_reduced

cd ../cov_to_cyto_reduced

rm -r 18-334_S164_*x_sorted.tab
rm -r 18-91_S160_*x_sorted.tab
rm -r 18-9_S159_*x_sorted.tab

rm -r 18-334_S164_sorted.cov
rm -r 18-91_S160_sorted.cov
rm -r 18-9_S159_sorted.cov
```

## 4) Create a file with positions found in all samples at specified coverage

`nano 5x_intersect.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

multiIntersectBed -i *_5x_sorted.tab > CpG.all.samps.5x_sorted.bed

cat CpG.all.samps.5x_sorted.bed | awk '$4 ==44' > CpG.filt.all.samps.5x_sorted.bed
```

`nano 10x_intersect.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

multiIntersectBed -i *_10x_sorted.tab > CpG.all.samps.10x_sorted.bed

cat CpG.all.samps.10x_sorted.bed | awk '$4 ==44' > CpG.filt.all.samps.10x_sorted.bed
```

## 5) Create bedgraphs post merge

5x bedgraph

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4}}' \
  > "${STEM}"_5x_sorted.bedgraph
done
```

10 x bedgraph

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4}}' \
  > "${STEM}"_10x_sorted.bedgraph
done
```

## 6) Use intersectBed to find where loci and genes intersect, allowing loci to be mapped to annotated genes

#### Filtered gff to only include gene positions modified.gff > gene.gff

```bash
cd /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1
awk '{if ($3 == "gene") {print}}' Pastreoides_all_v1.gff  > Pastreoides_all_v1.gene.gff
```

`nano 5x_intersectBed.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

for i in *5x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gene.gff \
  > ${i}_gene
done
```

`nano 10x_intersectBed.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

for i in *10x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gene.gff \
  > ${i}_gene
done
```

## 7) Intersect with file to subset only those positions found in all samples

`nano 5x_intersect_final.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

for i in *_5x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.all.samps.5x_sorted.bed \
  > ${i}_CpG_5x_enrichment.bed
done
```

`wc -l *5x_enrichment.bed`

```
78052 18-106_S163_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-118_S162_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-130_S172_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-142_S189_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-167_S166_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-178_S191_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-190_S186_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-202_S188_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-20_S202_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-227_S170_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-239_S185_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-250_S195_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-262_S179_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-311_S187_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-322_S180_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-32_S178_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-346_S193_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-358_S201_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-370_S171_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-394_S192_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-406_S177_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-418_S196_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-442_S165_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-44_S198_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-454_S197_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-466_S199_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-55_S190_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-67_S176_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 18-79_S181_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-1029_S183_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-1038_S184_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-1053_S167_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-1059_S175_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-1093_S168_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-1257_S205_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-1263_S173_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-562_S174_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-571_S194_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-661_S182_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-704_S169_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-728_S161_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-862_S200_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-924_S204_5x_sorted.tab_gene_CpG_5x_enrichment.bed
78052 L-933_S203_5x_sorted.tab_gene_CpG_5x_enrichment.bed
3434288 total
```

`nano 10x_intersect_final.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_reduced
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

for i in *_10x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.all.samps.10x_sorted.bed \
  > ${i}_CpG_10x_enrichment.bed
done
```

`wc -l *10x_enrichment.bed`

```
2143 18-106_S163_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-118_S162_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-130_S172_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-142_S189_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-167_S166_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-178_S191_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-190_S186_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-202_S188_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-20_S202_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-227_S170_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-239_S185_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-250_S195_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-262_S179_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-311_S187_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-322_S180_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-32_S178_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-346_S193_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-358_S201_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-370_S171_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-394_S192_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-406_S177_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-418_S196_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-442_S165_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-44_S198_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-454_S197_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-466_S199_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-55_S190_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-67_S176_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 18-79_S181_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-1029_S183_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-1038_S184_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-1053_S167_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-1059_S175_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-1093_S168_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-1257_S205_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-1263_S173_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-562_S174_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-571_S194_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-661_S182_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-704_S169_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-728_S161_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-862_S200_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-924_S204_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2143 L-933_S203_10x_sorted.tab_gene_CpG_10x_enrichment.bed
94292 total
```

# Pipeline for each group comparison

## Lifestage comparison

Samples: Paired adult and larval samples

`mkdir cov_to_cyto_lifestage`


## Running all steps above in one script

#### 5x coverage

```bash
cp ../cov_to_cyto_reduced/18-118_S162_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-202_S188_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-322_S180_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-358_S201_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-106_S163_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-190_S186_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-370_S171_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-454_S197_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-130_S172_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-142_S189_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-418_S196_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-55_S190_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-167_S166_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-227_S170_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-406_S177_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1029_S183_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-728_S161_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-924_S204_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-933_S203_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1053_S167_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1257_S205_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-704_S169_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-862_S200_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1038_S184_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1263_S173_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-562_S174_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-661_S182_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1059_S175_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1093_S168_5x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-571_S194_5x_sorted.tab ./
```

`nano runall_meth_5x.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_lifestage
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

#Intersect positions found in all samples at a specific coverage

multiIntersectBed -i *_5x_sorted.tab > CpG.lifestage.samps.5x_sorted.bed

cat CpG.lifestage.samps.5x_sorted.bed | awk '$4 ==30' > CpG.filt.lifestage.samps.5x_sorted.bed

# Use intersectBed to find where loci and genes intersect, allowing loci to be mapped to annotated genes
#wb: Print all lines in the second file
#a: file
#b: annotated gene list
#Save output in a new file that has the same base name

for i in *5x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gene.gff \
  > ${i}_gene
done

#Intersect with file to subset only those positions found in all samples

for i in *_5x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.lifestage.samps.5x_sorted.bed \
  > ${i}_CpG_5x_enrichment.bed
done
```

## 5) Create bedgraphs post merge

5x bedgraph

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4}}' \
  > "${STEM}"_5x_sorted.bedgraph
done
```

10 x bedgraph

```
for f in *_sorted.cov
do
  STEM=$(basename "${f}" _sorted.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 10) {print $1, $2, $3, $4}}' \
  > "${STEM}"_10x_sorted.bedgraph
done




`wc -l *5x_enrichment.bed`

```
103636 18-106_S163_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-118_S162_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-130_S172_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-142_S189_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-167_S166_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-190_S186_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-202_S188_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-227_S170_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-322_S180_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-358_S201_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-370_S171_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-406_S177_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-418_S196_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-454_S197_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 18-55_S190_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1029_S183_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1038_S184_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1053_S167_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1059_S175_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1093_S168_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1257_S205_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-1263_S173_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-562_S174_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-571_S194_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-661_S182_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-704_S169_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-728_S161_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-862_S200_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-924_S204_5x_sorted.tab_gene_CpG_5x_enrichment.bed
103636 L-933_S203_5x_sorted.tab_gene_CpG_5x_enrichment.bed
 3109080 total
```

### 10x coverage

```bash
cp ../cov_to_cyto_reduced/18-118_S162_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-202_S188_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-322_S180_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-358_S201_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-106_S163_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-190_S186_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-370_S171_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-454_S197_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-130_S172_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-142_S189_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-418_S196_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-55_S190_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-167_S166_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-227_S170_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/18-406_S177_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1029_S183_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-728_S161_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-924_S204_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-933_S203_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1053_S167_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1257_S205_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-704_S169_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-862_S200_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1038_S184_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1263_S173_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-562_S174_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-661_S182_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1059_S175_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-1093_S168_10x_sorted.tab ./
cp ../cov_to_cyto_reduced/L-571_S194_10x_sorted.tab ./
```

`nano runall_meth_10x.sh`

```bash
#!/bin/bash
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=500GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_lifestage
#SBATCH --cpus-per-task=3

# load modules needed

module load BEDTools/2.27.1-foss-2018b

#Intersect positions found in all samples at a specific coverage

multiIntersectBed -i *_10x_sorted.tab > CpG.lifestage.samps.10x_sorted.bed

cat CpG.lifestage.samps.10x_sorted.bed | awk '$4 ==30' > CpG.filt.lifestage.samps.10x_sorted.bed

# Use intersectBed to find where loci and genes intersect, allowing loci to be mapped to annotated genes
#wb: Print all lines in the second file
#a: file
#b: annotated gene list
#Save output in a new file that has the same base name

for i in *10x_sorted.tab
do
  intersectBed \
  -wb \
  -a ${i} \
  -b /data/putnamlab/kevin_wong1/Past_Genome/past_struc_annotations_v1/Pastreoides_all_v1.gene.gff \
  > ${i}_gene
done

#Intersect with file to subset only those positions found in all samples

for i in *_10x_sorted.tab_gene
do
  intersectBed \
  -a ${i} \
  -b CpG.filt.lifestage.samps.10x_sorted.bed \
  > ${i}_CpG_10x_enrichment.bed
done
```


`wc -l *10x_enrichment.bed`

```
2687 18-106_S163_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-118_S162_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-130_S172_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-142_S189_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-167_S166_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-190_S186_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-202_S188_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-227_S170_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-322_S180_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-358_S201_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-370_S171_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-406_S177_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-418_S196_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-454_S197_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 18-55_S190_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1029_S183_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1038_S184_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1053_S167_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1059_S175_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1093_S168_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1257_S205_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-1263_S173_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-562_S174_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-571_S194_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-661_S182_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-704_S169_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-728_S161_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-862_S200_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-924_S204_10x_sorted.tab_gene_CpG_10x_enrichment.bed
2687 L-933_S203_10x_sorted.tab_gene_CpG_10x_enrichment.bed
80610 total
```

Exporting bedfiles to local computer:

```
scp 'kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_lifestage/*_5x_enrichment.bed' ~/MyProjects/Thermal_Transplant_Molecular/output/WGBS/meth_counts_5x

scp 'kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/Thermal_Transplant_WGBS/Past_WGBS/cov_to_cyto_lifestage/*_10x_enrichment.bed' ~/MyProjects/Thermal_Transplant_Molecular/output/WGBS/meth_counts_10x

```
