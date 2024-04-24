---
layout: post
title: CellRanger on Ehrens scRNAseq data
date: '2024-04-24'
categories: Analysis
tags: scRNAseq
---


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

`nano mkref.job`

```bash
#!/bin/bash

#BSUB -J scSeq_mkref_acer
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 12:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o acer_mkref_cat.out
#BSUB -e acer_mkref_cat.err

###################################################################

cellranger mkref \
--genome=acer_K2 --fasta=GCA_032359415.1_NEU_Acer_K2_genomic.fa --genes=genomic.gtf

```


## Count on 18003-001 library

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

cellranger count \
 --id=18003-001 \
 --transcriptome=/nethome/kxw755/GCA_032359415.1/acer_K2 \
 --fastqs=/nethome/kxw755/20240418_Ehrens-18003-001_GEX3_GEMX \
 --sample=Ehrens-18003-001_GEX3_GEMX
 ```

`bsub < count_001.job`

got the following error
 ```
 [error] Pipestance failed. Error log at:
18003-001/SC_RNA_COUNTER_CS/SC_MULTI_CORE/MULTI_CHEMISTRY_DETECTOR/_GEM_WELL_CHEMISTRY_DETECTOR/DETECT_COUNT_CHEMISTRY/fork0/chnk0-u55412970cb/_errors

Log message:
An extremely low rate of correct barcodes was observed for all the candidate chemistry choices for the input: Sample Ehrens-18003-001_GEX3_GEMX in "/projectnb/pegasus/nethome/kxw755/20240418_Ehrens-18003-001_GEX3_GEMX". Please check your input data.
- 1.1% for chemistry SC3Pv3
- 1.1% for chemistry SC3Pv3HT
- 0.0% for chemistry SC5P-R2
- 0.0% for chemistry SC3Pv2
- 0.0% for chemistry SC3Pv3LT
```

## Count on 18003-002 library

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

cellranger count \
 --id=18003-002 \
 --transcriptome=/nethome/kxw755/GCA_032359415.1/acer_K2 \
 --fastqs=/nethome/kxw755/20240418_Ehrens-18003-002_GEX3_GEMX \
 --sample=Ehrens-18003-002_GEX3_GEMX
 ```

`bsub < count_002.job`

Again, I got this error. 
```
Log message:
An extremely low rate of correct barcodes was observed for all the candidate chemistry choices for the input: Sample Ehrens-18003-002_GEX3_GEMX in "/projectnb/pegasu
s/nethome/kxw755/20240418_Ehrens-18003-002_GEX3_GEMX". Please check your input data.
- 0.7% for chemistry SC3Pv3
- 0.7% for chemistry SC3Pv3HT
- 0.1% for chemistry SC5P-R2
- 0.1% for chemistry SC3Pv2
- 0.0% for chemistry SC3Pv3LT
```
