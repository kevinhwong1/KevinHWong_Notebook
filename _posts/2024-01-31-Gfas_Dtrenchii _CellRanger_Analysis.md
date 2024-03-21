---
layout: post
title: Gfas_Dtrenchii _CellRanger_Analysis
date: '2024-01-31'
categories: Analysis
tags: Cnidarin immunity lab
---

# Goal: To re-analyze the 10X single cell data with the combined genomic references of Galaxea fasicularis and Durisdinium trenchii

**Reference Genomes:**
- Galaxea fasicularis v1 
    - From Reef Genomics
    - [Download site](http://gfas.reefgenomics.org/)
- Durisdinum trenchii (SCF082)
    - [Dougan et al. 2022](https://www.biorxiv.org/content/10.1101/2022.04.10.487810v1.full.pdf)
    - [Download site](https://espace.library.uq.edu.au/view/UQ:27da3e7)
 
# D.trenchii reference 

`gunzip Dtrenchii_SCF082_ANNOT_gff.gz`

`mv Dtrenchii_SCF082_ANNOT_gff Dtrenchii_SCF082_ANNOT.gff3`


# Converting this GFF3 file to GTF:

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

# Making combined genome reference

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

cellranger mkref \
--genome=gfas_1.0 --fasta=/nethome/kxw755/Gfas_v1/gfas_final_1.0.fasta --genes=/nethome/kxw755/Gfas_v1/gfas_1.0.filtered.gtf \
--genome=Dtrenchii_SCF082 --fasta=/nethome/kxw755/Dtrenchii_SCF082/Dtrenchii_SCF082_ASSEMBLY.fasta --genes=/nethome/kxw755/Dtrenchii_SCF082/Dtrenchii_SCF082_ANNOT.gtf
```

`bsub < mkgtf.job`

Combined fasta genome path:

`/nethome/kxw755/mkref_gfas_dtre/gfas_1.0_and_Dtrenchii_SCF082/fasta/genome.fa`

# Rerunning Count with the combined genome

## Run 1

cd /nethome/kxw755/20231004_SingleCell_DG/

`nano count_W045_deep_combgenome.job`

```bash
#BSUB -J count_W045_deep_combgenome
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
 --id=W-045_1_deep_combgenome \
 --transcriptome=/nethome/kxw755/mkref_gfas_dtre/gfas_1.0_and_Dtrenchii_SCF082/ \
 --fastqs=/nethome/kxw755/20231004_SingleCell_DG \
 --sample=AndradeRodriguez-15275-001_GEX3
```

`bsub < count_W045_deep_combgenome.job`
