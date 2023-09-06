---
layout: post
title: Gfas Dark Genes Single cell preliminary analysis
date: '2023-09-05'
categories: Analysis
tags: 
---

**Goal:**
- To assess the quality of our single cell data on a control (ambient) sample

**Data:**
- Single cell suspension prepartion protocol (add)
- This is coral W-045
- From the 10X Chromium sequencing (10X 3' v3.1), we performed shallow sequencing (500 reads per cell)
- This produced 10 files

**Reference Genomes:**
- Galaxea fasicularis v1 
    - From Reef Genomics
    - [Download site](http://gfas.reefgenomics.org/)
- Durisdinum trenchii (SCF082)
    - [Dougan et al. 2022](https://www.biorxiv.org/content/10.1101/2022.04.10.487810v1.full.pdf)
    - [Download site](https://espace.library.uq.edu.au/view/UQ:27da3e7)

## Transfer data from Box to Pegasus 

1. I manually downloaded the files from Box to my hard drive. 

2. md5 on hard drive files

#### All Cells:

`$ cd 20230905_SingleCell_DG/`

`$ md5 *.gz`

```bash 
MD5 (AndradeRodriguez-15275-001_GEX3_S9_L001_R1_001.fastq.gz) = a41d311fed6542869cbfca5d933a3d5d

MD5 (AndradeRodriguez-15275-001_GEX3_S9_L001_R2_001.fastq.gz) = 8a70aaa04c6bffbf2bcb4576f4aeb9b3

MD5 (AndradeRodriguez-15275-001_GEX3_S9_L002_R1_001.fastq.gz) = 06f5e9b6b9b6713fe65376dc53539488

MD5 (AndradeRodriguez-15275-001_GEX3_S9_L002_R2_001.fastq.gz) = 7dcffd1668d4b335061bdaded9d9e898

MD5 (AndradeRodriguez-15275-001_GEX3_S9_L003_R1_001.fastq.gz) = 8e801a9bb82cedac7e6b3470e27e4095

MD5 (AndradeRodriguez-15275-001_GEX3_S9_L003_R2_001.fastq.gz) = 4eff72dca6601f22134982a6c0a3577f

MD5 (AndradeRodriguez-15275-001_GEX3_S9_L004_R1_001.fastq.gz) = f9b85b014e123dcadf854ea037cf7f75

MD5 (AndradeRodriguez-15275-001_GEX3_S9_L004_R2_001.fastq.gz) = b953a5d66ba7b47e19d3485cdc3a774c

MD5 (AndradeRodriguez-15275-001_GEX3_S9_L005_R1_001.fastq.gz) = 3df3e6b965af17f451df16dab1604110

MD5 (AndradeRodriguez-15275-001_GEX3_S9_L005_R2_001.fastq.gz) = 4a33ee1e8dda001a9b1d7463afbe726b
```

3. Transfer data from hard drive to Pegasus

Raw fastq files:

`$ scp -r /Volumes/Passport/Kevin/Andrade_Rodriguez-21-07-15275-396479269/20230905_SingleCell_DG kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/`

Galaxea fasicularis genome files:

`$ scp -r /Volumes/Passport/Kevin/Gfas_v1 kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/`

Durisdinum trenchii genome files:

`$ scp -r /Volumes/Passport/Kevin/Dtrenchii_SCF082 kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/`

4. md5 file for Pegasus files

`$ cd 20230905_SingleCell_DG/`

`$ md5sum *.gz `

```bash
a41d311fed6542869cbfca5d933a3d5d  AndradeRodriguez-15275-001_GEX3_S9_L001_R1_001.fastq.gz
8a70aaa04c6bffbf2bcb4576f4aeb9b3  AndradeRodriguez-15275-001_GEX3_S9_L001_R2_001.fastq.gz
06f5e9b6b9b6713fe65376dc53539488  AndradeRodriguez-15275-001_GEX3_S9_L002_R1_001.fastq.gz
7dcffd1668d4b335061bdaded9d9e898  AndradeRodriguez-15275-001_GEX3_S9_L002_R2_001.fastq.gz
8e801a9bb82cedac7e6b3470e27e4095  AndradeRodriguez-15275-001_GEX3_S9_L003_R1_001.fastq.gz
4eff72dca6601f22134982a6c0a3577f  AndradeRodriguez-15275-001_GEX3_S9_L003_R2_001.fastq.gz
f9b85b014e123dcadf854ea037cf7f75  AndradeRodriguez-15275-001_GEX3_S9_L004_R1_001.fastq.gz
b953a5d66ba7b47e19d3485cdc3a774c  AndradeRodriguez-15275-001_GEX3_S9_L004_R2_001.fastq.gz
3df3e6b965af17f451df16dab1604110  AndradeRodriguez-15275-001_GEX3_S9_L005_R1_001.fastq.gz
4a33ee1e8dda001a9b1d7463afbe726b  AndradeRodriguez-15275-001_GEX3_S9_L005_R2_001.fastq.gz
```


## Build a custom reference

- [Tutorial](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_mr)

1. Decompress files

```bash
$ gunzip gfas_1.0.genes.gff3.gz

$ gunzip gfas_final_1.0.fasta.gz
```

I have to convert this GFF3 file to GTF:

- [Resource link](https://www.biostars.org/p/45791/)

`$ nano convert.job`

```bash
#!/bin/bash

#BSUB -J scSeq_convertGTF
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 12:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o gfas_convertGTF.out
#BSUB -e gfas_convertGTF.err

###################################################################

module load cufflinks/2.2.1 

gffread gfas_1.0.genes.gff3 -T -o gfas_1.0.gtf

```

`$ bsub < convert.job`


2. Filter GTF

`$ nano mkgtf.job`

```bash
#!/bin/bash

#BSUB -J scSeq_mkgtf_dg1
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 12:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o gfas_mkgtf_v1.out
#BSUB -e gfas_mkgtf_v1.err

###################################################################

cellranger mkgtf \
  gfas_1.0.gtf \
  gfas_1.0.filtered.gtf \
  --attribute=gene_biotype:protein_coding
```

`$ bsub < mkgtf.job`

3. Make Reference

`$ nano mkref.job`

```bash
#!/bin/bash

#BSUB -J scSeq_mkref_dg1
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 12:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o gfas_mkref_v1.out
#BSUB -e gfas_mkref_v1.err

###################################################################

cellranger mkref \
--genome=gfas \
--fasta=gfas_final_1.0.fasta \
--genes=gfas_1.0.filtered.gtf
```

`$ bsub < mkref.job`

## Count 

`$ cd /nethome/kxw755/20230905_SingleCell_DG/raw_data`

`$ nano count_W045_1.job`

```bash
#BSUB -J count_W045
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 120:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count.out
#BSUB -e count.err
#BSUB -B
#BSUB -N
###################################################################

cellranger count \
 --id=W-045_1 \
 --transcriptome=/nethome/kxw755/Gfas_v1/gfas \
 --fastqs=/nethome/kxw755/20230905_SingleCell_DG/raw_data \
 --sample=AndradeRodriguez-15275-001_GEX3
```

`$ bsub < count_W045_1.job`

## Exporting files

`scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20230905_SingleCell_DG/raw_data/W-045_1/outs/filtered_feature_bc_matrix.h5 /Users/kevinwong/MyProjects/DarkGenes_Bleaching_Comparison/output/CellRanger/W-045_1_raw_feature_bc_matrix.h5`

`scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20230905_SingleCell_DG/raw_data/W-045_1/outs/web_summary.html /Users/kevinwong/MyProjects/DarkGenes_Bleaching_Comparison/output/CellRanger/W-045_1_web_summary.html`

## Summary

![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/20230905_W045_scRNAseq_Summary.png?raw=true)

![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/20230905_W045_scRNAseq_Summary2.png?raw=true)