---
layout: post
title: CellRanger on Ehrens scRNAseq data
date: '2024-04-24'
categories: Analysis
tags: scRNAseq
---

# Background

* This post describes the initial QC analysis for scRNAseq on Acropora cervicornis and Pocillopora damicornis using the new 10X GEMX cel capture technology. FYI to analyze this data, you must use `cellranger/8.0.0`. 

* For each species, there are two separate libraries: 001 = All cells, 002 = ALDH+ cells

# Acropora cervicornis analysis

## Downloading Acropora cervicornis genome

https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_032359415.1/

`scp -r GCA_032359415.1.zip kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/`

`unzip GCA_032359415.1.zip`

I need to convert the .fna file to .fa for STAR to work properly 

```bash
cd GCA_032359415.1
cp GCA_032359415.1_NEU_Acer_K2_genomic.fna GCA_032359415.1_NEU_Acer_K2_genomic.fa
```

## Make reference 

`nano mkref2.job`

```bash
#!/bin/bash

#BSUB -J scSeq_mkref_acer
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 12:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o acer_mkref_2.out
#BSUB -e acer_mkref_2.err

###################################################################

module load cellranger/8.0.0

cellranger mkref \
--genome=acer_K2_2 --fasta=GCA_032359415.1_NEU_Acer_K2_genomic.fa --genes=genomic.gtf

```

`bsub < mkref2.job `

/nethome/kxw755/GCA_032359415.1/acer_K2


## Count

### Ehrens-18003-001 library (All Cells)

`cd /nethome/kxw755/20240418_Ehrens-18003-001_GEX3_GEMX`

`nano count_001.job`

```bash
#BSUB -J count_18003-001
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_18003-001.out
#BSUB -e count_18003-001.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=18003-001 \
 --transcriptome=/nethome/kxw755/GCA_032359415.1/acer_K2_2   \
 --fastqs=/nethome/kxw755/20240418_Ehrens-18003-001_GEX3_GEMX \
 --sample=Ehrens-18003-001_GEX3_GEMX \
 --create-bam=true
 ```

 `bsub < count_001.job`

Export 

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240418_Ehrens-18003-001_GEX3_GEMX/18003-001/outs/web_summary.html /Users/kxw755/Desktop/Ehrens_StemCells/20240418_Ehrens-18003-001_web_summary.html

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240418_Ehrens-18003-001_GEX3_GEMX/18003-001/outs/filtered_feature_bc_matrix.h5 /Users/kxw755/Desktop/Ehrens_StemCells/20240418_Ehrens-18003-001_filtered_feature_bc_matrix.h5
```

![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/20240503_18003-001_Stats.png)
![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/20240503_18003-001_tsne.png)

### Ehrens-18003-002 library (ALDH+)

`cd /nethome/kxw755/20240418_Ehrens-18003-002_GEX3_GEMX`

`nano count_002.job`

```bash
#BSUB -J count_18003-002
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_18003-002.out
#BSUB -e count_18003-002.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=18003-002 \
 --transcriptome=/nethome/kxw755/GCA_032359415.1/acer_K2_2 \
 --fastqs=/nethome/kxw755/20240418_Ehrens-18003-002_GEX3_GEMX \
 --sample=Ehrens-18003-002_GEX3_GEMX \
 --create-bam=true
 ```

`bsub < count_002.job`

Export

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240418_Ehrens-18003-002_GEX3_GEMX/18003-002/outs/web_summary.html /Users/kxw755/Desktop/Ehrens_StemCells/20240418_Ehrens-18003-002_web_summary.html

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240418_Ehrens-18003-002_GEX3_GEMX/18003-002/outs/filtered_feature_bc_matrix.h5 /Users/kxw755/Desktop/Ehrens_StemCells/20240418_Ehrens-18003-002_filtered_feature_bc_matrix.h5
```

![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/20240503_18003-002_Stats.png)
![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/20240503_18003-002_tsne.png)


# Pocillopora damicornis analysis

## Download Pocillopora damicornis genome

https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003704095.1/

`scp -r pdam_genome kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/`

I need to convert the .fna file to .fa for STAR to work properly 

```bash
cd GCA_032359415.1
cp GCA_003704095.1_ASM370409v1_genomic.fna GCA_003704095.1_ASM370409v1_genomic.fa
```

## Make reference 

`nano mkref.job`

```bash
#!/bin/bash

#BSUB -J scSeq_mkref_pdam
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 12:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o pdam_mkref_2.out
#BSUB -e pdam_mkref_2.err
#BSUB -B
#BSUB -N

###################################################################

module load cellranger/8.0.0

cellranger mkref \
--genome=pdam --fasta=GCA_003704095.1_ASM370409v1_genomic.fa --genes=genomic.gtf

```

`bsub < mkref.job`


## Count

### Ehrens-18114-001 library (All Cells)

`cd /nethome/kxw755/20240501_Ehrens-18114-001_GEX3_GEMX`

`nano count_001.job`

```bash
#BSUB -J count_18114-001
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_18114-001.out
#BSUB -e count_18114-001.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=18114-001 \
 --transcriptome=/nethome/kxw755/pdam_genome/pdam   \
 --fastqs=/nethome/kxw755/20240501_Ehrens-18114-001_GEX3_GEMX \
 --sample=Ehrens-18114-001_GEX3_GEMX \
 --create-bam=false
 ```

 `bsub < count_001.job`

Export

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240501_Ehrens-18114-001_GEX3_GEMX/18114-001/outs/web_summary.html /Users/kxw755/Desktop/Ehrens_StemCells/20240501_Ehrens-18114-001_web_summary.html

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240501_Ehrens-18114-001_GEX3_GEMX/18114-001/outs/filtered_feature_bc_matrix.h5 /Users/kxw755/Desktop/Ehrens_StemCells/20240501_Ehrens-18114-001_filtered_feature_bc_matrix.h5
```

![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/20240503_18114-001_Stats.png)
![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/20240503_18114-001_tsne.png)

### Ehrens-18114-002 library (ALDH+)

`cd /nethome/kxw755/20240501_Ehrens-18114-002_GEX3_GEMX`

`nano count_002.job`

```bash
#BSUB -J count_18114-002
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_18114-002.out
#BSUB -e count_18114-002.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=18114-002 \
 --transcriptome=/nethome/kxw755/pdam_genome/pdam   \
 --fastqs=/nethome/kxw755/20240501_Ehrens-18114-002_GEX3_GEMX \
 --sample=Ehrens-18114-002_GEX3_GEMX \
 --create-bam=false
 ```

 `bsub < count_002.job`

Export

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240501_Ehrens-18114-002_GEX3_GEMX/18114-002/outs/web_summary.html /Users/kxw755/Desktop/Ehrens_StemCells/20240501_Ehrens-18114-002_web_summary.html

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240501_Ehrens-18114-002_GEX3_GEMX/18114-002/outs/filtered_feature_bc_matrix.h5 /Users/kxw755/Desktop/Ehrens_StemCells/20240501_Ehrens-18114-002_filtered_feature_bc_matrix.h5
```

![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/20240503_18114-002_Stats.png)
![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/20240503_18114-002_tsne.png)
