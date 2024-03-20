---
layout: post
title: 20240122 Gfas single cell fresh vs frozen analysis
date: '2024-01-22'
categories: Analysis
tags: single cell
---

# Analyzing fresh vs frozen single cell preprations 

**Goal:**
- To assess the quality of our single cell data on a control (ambient) samples that were frozen vs fresh

**Data:**
- Single cell suspension prepartion protocol (add)
- This is coral W-045
- From the 10X Chromium sequencing (10X 3' v3.1)
- We also did an experimental approach to overload the beads during GEM capture (15000 when targeting 10000 cells)
- Sample numbers: 
    - 1. Frozen control
    - 2. Fresh sample 

**Reference Genomes:**
- Galaxea fasicularis v1 
    - From Reef Genomics
    - [Download site](http://gfas.reefgenomics.org/)
- Durisdinum trenchii (SCF082)
    - [Dougan et al. 2022](https://www.biorxiv.org/content/10.1101/2022.04.10.487810v1.full.pdf)
    - [Download site](https://espace.library.uq.edu.au/view/UQ:27da3e7)

## Transfer data to Pegasus 

1. I downloaded the files to the lab computer though the BaseSpace downloader and uploaded to pegasus (dark genes folder). 

`scp -r /Users/cnidarianimmunity/Desktop/KHW/Andrade_Rodriguez-25-09-16001-406992761/BCL_Convert_01_19_2024_9_38_28-714586393 kxw755@pegasus.ccs.miami.edu:/scratch/projects/dark_genes`

## Moving sample files into the same folder

```bash
mkdir AndradeRodriguez-16001-001

mv AndradeRodriguez-16001-001_GEX3_L1-ds.74ad8ffaed9f491181d8662a2ab2f80c/*gz AndradeRodriguez-16001-001/

mv AndradeRodriguez-16001-001_GEX3_L2-ds.fab93f6d2b9f46238661b92fcadf5fd6/*gz AndradeRodriguez-16001-001/

mv AndradeRodriguez-16001-001_GEX3_L3-ds.508128be6fcc421eb4afe9eabf53e973/*gz AndradeRodriguez-16001-001/
```

```bash
mkdir AndradeRodriguez-16001-002

mv AndradeRodriguez-16001-002_GEX3_L1-ds.24b3e322d8b14a8b8969dce5df5aa989/*gz AndradeRodriguez-16001-002/

mv AndradeRodriguez-16001-002_GEX3_L2-ds.38c544ab4d6a4f29a7504698bbe8f264/*gz AndradeRodriguez-16001-002/

mv AndradeRodriguez-16001-002_GEX3_L3-ds.065684cd381348cdaac690ec6ce589c2/*gz AndradeRodriguez-16001-002/
```


# Count 1 (Frozen)

`nano count_1.job`

```bash
#BSUB -J count_1
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
 --id=1_frozen \
 --transcriptome=/nethome/kxw755/Gfas_v1/gfas \
 --fastqs=/scratch/projects/dark_genes/BCL_Convert_01_19_2024_9_38_28-714586393/AndradeRodriguez-16001-001 \
 --sample=AndradeRodriguez-16001-001_GEX3
```

`bsub < count_1.job`

export

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/scratch/projects/dark_genes/BCL_Convert_01_19_2024_9_38_28-714586393/AndradeRodriguez-16001-001/1_frozen/outs/filtered_feature_bc_matrix.h5 /Users/kevinwong/MyProjects/DarkGenes_Bleaching_Comparison/output/CellRanger/20240123_1_frozen_filt_feature_bc_matrix.h5

scp -r kxw755@pegasus.ccs.miami.edu:/scratch/projects/dark_genes/BCL_Convert_01_19_2024_9_38_28-714586393/AndradeRodriguez-16001-001/1_frozen/outs/web_summary.html /Users/kevinwong/MyProjects/DarkGenes_Bleaching_Comparison/output/CellRanger/20240123_1_frozen_web_summary.html

scp -r kxw755@pegasus.ccs.miami.edu:/scratch/projects/dark_genes/BCL_Convert_01_19_2024_9_38_28-714586393/AndradeRodriguez-16001-001/1_frozen/outs/metrics_summary.csv /Users/kevinwong/MyProjects/DarkGenes_Bleaching_Comparison/output/CellRanger/20240123_1_frozen_metrics_summary.csv
```


# Count 2 (Fresh)

`nano count_2.job`

```bash
#BSUB -J count_2
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
 --id=2_fresh \
 --transcriptome=/nethome/kxw755/Gfas_v1/gfas \
 --fastqs=/scratch/projects/dark_genes/BCL_Convert_01_19_2024_9_38_28-714586393/AndradeRodriguez-16001-002 \
 --sample=AndradeRodriguez-16001-002_GEX3
```

`bsub < count_2.job`

export

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/scratch/projects/dark_genes/BCL_Convert_01_19_2024_9_38_28-714586393/AndradeRodriguez-16001-002/2_fresh/outs/filtered_feature_bc_matrix.h5 /Users/kevinwong/MyProjects/DarkGenes_Bleaching_Comparison/output/CellRanger/20240123_2_fresh_filt_feature_bc_matrix.h5

scp -r kxw755@pegasus.ccs.miami.edu:/scratch/projects/dark_genes/BCL_Convert_01_19_2024_9_38_28-714586393/AndradeRodriguez-16001-002/2_fresh/outs/web_summary.html /Users/kevinwong/MyProjects/DarkGenes_Bleaching_Comparison/output/CellRanger/20240123_2_fresh_web_summary.html

scp -r kxw755@pegasus.ccs.miami.edu:/scratch/projects/dark_genes/BCL_Convert_01_19_2024_9_38_28-714586393/AndradeRodriguez-16001-002/2_fresh/outs/metrics_summary.csv /Users/kevinwong/MyProjects/DarkGenes_Bleaching_Comparison/output/CellRanger/20240123_2_fresh_metrics_summary.csv
```

# D.trenchii reference 

`gunzip Dtrenchii_SCF082_ANNOT_gff.gz`

`mv Dtrenchii_SCF082_ANNOT_gff Dtrenchii_SCF082_ANNOT.gff3`


I have to convert this GFF3 file to GTF:

Resource link

`nano convert.job`

```bash
#!/bin/bash

#BSUB -J scSeq_convertGTF
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 12:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o dtre_convertGTF.out
#BSUB -e dtre_convertGTF.err

###################################################################

module load cufflinks/2.2.1 

gffread Dtrenchii_SCF082_ANNOT.gff3 -T -o Dtrenchii_SCF082_ANNOT.gtf
```

`bsub < convert.job`

# 

`gunzip Dtrenchii_SCF082_ASSEMBLY_fasta.gz`

`mv Dtrenchii_SCF082_ASSEMBLY_fasta Dtrenchii_SCF082_ASSEMBLY.fasta`

```bash
cd ..
mkdir mkref_gfas_dtre
```

`nano mkgtf.job`

```bash
#!/bin/bash

#BSUB -J scSeq_mkgtf_dg1
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 12:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o gfas_dtre_mkgtf_v1.out
#BSUB -e gfas_dtre_mkgtf_v1.err

###################################################################

cellranger mkref 
--genome=gfas_1.0 \
--fasta=/nethome/kxw755/Gfas_v1/gfas_final_1.0.fasta \
--genes=/nethome/kxw755/Gfas_v1/gfas_1.0.filtered.gtf \
--genome=Dtrenchii_SCF082 \
--fasta=/nethome/kxw755/Dtrenchii_SCF082/Dtrenchii_SCF082_ASSEMBLY.fasta \ 
--genes=/nethome/kxw755/Dtrenchii_SCF082/Dtrenchii_SCF082_ANNOT.gtf
```

`bsub < mkgtf.job`

/nethome/kxw755/Dtrenchii_SCF082/Dtrenchii_SCF082_ANNOT.gtf
/nethome/kxw755/Dtrenchii_SCF082/Dtrenchii_SCF082_ASSEMBLY.fasta

/nethome/kxw755/Gfas_v1/gfas_1.0.filtered.gtf
/nethome/kxw755/Gfas_v1/gfas_final_1.0.fasta