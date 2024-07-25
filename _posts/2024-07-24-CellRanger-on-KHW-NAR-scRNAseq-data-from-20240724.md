---
layout: post
title: CellRanger on KHW NAR scRNAseq data from 20240724
date: '2024-07-24'
categories: Analysis
tags: scRNAseq
---

# Sample information

| Sample # | Customer Sample Name |       BaseSpace Sample ID       | Flow Cell Type |    Read Type   | Total number of clusters (reads) per sample |
|:--------:|:--------------------:|:-------------------------------:|:--------------:|:--------------:|:-------------------------------------------:|
|     1    |     gfas_control     | AndradeRodriguez-19876-001_GEX3 |     10B-300    | paired-end,151 |                 581,421,834                 |
|     2    |      gfas_bleach     | AndradeRodriguez-19876-002_GEX3 |     10B-300    | paired-end,151 |                 636,439,053                 |


# Upload files

```bash
scp -r /Users/cnidarianimmunity/Desktop/KHW/Andrade_Rodriguez-01-07-19876-425120288/2024* kxw755@pegasus.ccs.miami.edu:/nethome/kxw755
```

# Count 001

`nano count_001.job`

```bash
#BSUB -J count_19876-001
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_19876-001.out
#BSUB -e count_19876-001.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=19876-001 \
 --transcriptome=/nethome/kxw755/Gfas_v1/gfas \
 --fastqs=/nethome/kxw755/20240724_AndradeRodriguez-19876-001 \
 --sample=AndradeRodriguez-19876-001_GEX3 \
 --create-bam=true
```

`bsub < count_001.job`

Export

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-001/19876-001/outs/web_summary.html /Users/kxw755/Desktop/DarkGenes/20240724_CellRanger/20240724_AndradeRodriguez-19876-001_web_summary.html

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-001/19876-001/outs/filtered_feature_bc_matrix.h5 /Users/kxw755/Desktop/DarkGenes/20240724_CellRanger/20240724_AndradeRodriguez-19876-001_filtered_feature_bc_matrix.h5

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-001/19876-001/outs/metrics_summary.csv /Users/kxw755/Desktop/DarkGenes/20240724_CellRanger/20240724_AndradeRodriguez-19876-001_metrics_summary.csv
```

# Count 002

`nano count_002.job`

```bash
#BSUB -J count_19876-002
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_19876-002.out
#BSUB -e count_19876-002.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=19876-002 \
 --transcriptome=/nethome/kxw755/Gfas_v1/gfas \
 --fastqs=/nethome/kxw755/20240724_AndradeRodriguez-19876-002 \
 --sample=AndradeRodriguez-19876-002_GEX3 \
 --create-bam=true
```

`bsub < count_002.job`

Export

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-002/19876-002/outs/web_summary.html /Users/kxw755/Desktop/DarkGenes/20240724_CellRanger/20240724_AndradeRodriguez-19876-002_web_summary.html

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-002/19876-002/outs/filtered_feature_bc_matrix.h5 /Users/kxw755/Desktop/DarkGenes/20240724_CellRanger/20240724_AndradeRodriguez-19876-002_filtered_feature_bc_matrix.h5

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-002/19876-002/outs/metrics_summary.csv /Users/kxw755/Desktop/DarkGenes/20240724_CellRanger/20240724_AndradeRodriguez-19876-002_metrics_summary.csv
```
