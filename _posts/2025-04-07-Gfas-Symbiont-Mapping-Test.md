---
layout: post
title: Gfas Symbiont Mapping Test
date: '2025-04-07'
categories: Analysis
tags: scRNAseq, gfas
---

# Test for which symbiont type is in the Gfas dataset

So I cannot combine all the genomes together (it is too big for cell ranger to process) so I will run the following 3 tests: 

1. cgor (C)
2. dtre (D)
3. Smic (A)


# Download Genomes

* Cladocopium goreaui 
    * by [Lu et al. 2018](https://pubmed.ncbi.nlm.nih.gov/30271976/)
    * [SCF055-01](https://genome.jgi.doe.gov/portal/Clago1/Clago1.download.html)
* Durisdinium trenchii 
    * https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1381693
* Symbiodinium microadriaticum
    * https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Symmic1

for each genome I will need a fasta and a GTF file

```bash
scp -r Clago1 kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/genomes
scp -r gfas_mitogenome kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/genomes
scp -r * kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/genomes/Symmic1
```

Activate conda environment

```bash
source anaconda3/bin/activate 
conda activate agat
```

unzip clago genome files and convert gff3 to gtf
```bash
cd/nethome/kxw755/genomes/Clago1
gunzip Clago1_AssemblyScaffolds.fasta.gz 
gunzip Clago1_GeneCatalog_genes_20200812.gff3.gz 

agat_convert_sp_gff2gtf.pl -i Clago1_GeneCatalog_genes_20200812.gff3 -o Clago1_GeneCatalog_genes_20200812.gtf
```

unzip dtre genome files and convert gff3 to gtf
```bash
/nethome/kxw755/genomes/Dtrenchii_SCF082
gunzip Dtrenchii_SCF082_ASSEMBLY_fasta.gz
gunzip Dtrenchii_SCF082_ANNOT_gff.gz 

mv Dtrenchii_SCF082_ASSEMBLY_fasta Dtrenchii_SCF082_ASSEMBLY.fasta
mv Dtrenchii_SCF082_ANNOT_gff Dtrenchii_SCF082_ANNOT.gff3

agat_convert_sp_gff2gtf.pl -i Dtrenchii_SCF082_ANNOT.gff3 -o Dtrenchii_SCF082_ANNOT.gtf
```

Convert Smic gff2 to gtf

```bash
agat_convert_sp_gff2gtf.pl -i Symmic1_GeneCatalog_20180603.gff3 -o Symmic1_GeneCatalog_20180603.gtf
```

# Index reference genomes with Cell Ranger *mkref*

`nano mkref_Smic.job`

```bash
#!/bin/bash

#BSUB -J smic_mkref
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o smic_mkref_%J.out
#BSUB -e smic_mkref_%J.err

###################################################################

module load cellranger/8.0.0

cellranger mkref \
--genome=smic_mkref \
--fasta=/nethome/kxw755/genomes/Smic_JGI/Symmic1_AssemblyScaffolds.fasta \
--genes=/nethome/kxw755/genomes/Smic_JGI/Symmic1_GeneCatalog_20180603.gtf
```

`nano mkref_dtre.job`

```bash
#!/bin/bash

#BSUB -J dtre_mkref
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o dtre_mkref_%J.out
#BSUB -e dtre_mkref_%J.err

###################################################################

module load cellranger/8.0.0

cellranger mkref \
--genome=dtre_mkref \
--fasta=/nethome/kxw755/genomes/Dtrenchii_SCF082/Dtrenchii_SCF082_ASSEMBLY.fasta \
--genes=/nethome/kxw755/genomes/carnegie_gfas/gfas_mito_Cgor_Dtre/processed_gtfs/dtre.prefixed.gtf
```

`nano mkref_cgor.job`

```bash
#!/bin/bash

#BSUB -J mkref_cgor
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o mkref_cgor_%J.out
#BSUB -e mkref_cgor_%J.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger mkref \
--genome=cgor_mkref \
--fasta=/nethome/kxw755/genomes/Clago1/Clago1_AssemblyScaffolds.fasta \
--genes=/nethome/kxw755/genomes/carnegie_gfas/gfas_mito_Cgor_Dtre/processed_gtfs/cgor.validated.gtf
```

# Run count on control sample (19876-001)

`nano count_001_smic.job`

```bash
#BSUB -J count_19876-001_smic
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_19876-001_smic.out
#BSUB -e count_19876-001_smic.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=19876-001_Smic \
 --transcriptome=/nethome/kxw755/genomes/Smic_JGI/smic_mkref \
 --fastqs=/nethome/kxw755/20240724_AndradeRodriguez-19876-001 \
 --sample=AndradeRodriguez-19876-001_GEX3 \
 --create-bam=true
```

`bsub < count_001_cgor.job`


`nano count_001_cgor.job`

```bash
#BSUB -J count_19876-001_cgor
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_19876-001_cgor.out
#BSUB -e count_19876-001_cgor.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=19876-001_Cgor \
 --transcriptome=/nethome/kxw755/genomes/Clago1/cgor_mkref \
 --fastqs=/nethome/kxw755/20240724_AndradeRodriguez-19876-001 \
 --sample=AndradeRodriguez-19876-001_GEX3 \
 --create-bam=true
```

`bsub < count_001_cgor.job`


`nano count_001_dtre.job`

```bash
#BSUB -J count_19876-001_dtre
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_19876-001_dtre.out
#BSUB -e count_19876-001_dtre.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=19876-001_Dtre \
 --transcriptome=/nethome/kxw755/genomes/Dtrenchii_SCF082/dtre_mkref \
 --fastqs=/nethome/kxw755/20240724_AndradeRodriguez-19876-001 \
 --sample=AndradeRodriguez-19876-001_GEX3 \
 --create-bam=true
```

`bsub < count_001_dtre.job`


# Eport files

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-001/19876-001_Cgor/outs/web_summary.html ./19876-001_Cgor_web_summary.html
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-001/19876-001_Cgor/outs/metrics_summary.csv ./19876-001_Cgor_metrics_summary.csv

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-001/19876-001_Dtre/outs/web_summary.html ./19876-001_Dtre_web_summary.html
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-001/19876-001_Dtre/outs/metrics_summary.csv ./19876-001_Dtre_metrics_summary.csv

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-001/19876-001_Smic/outs/web_summary.html ./19876-001_Smic_web_summary.html
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20240724_AndradeRodriguez-19876-001/19876-001_Smic/outs/metrics_summary.csv ./19876-001_Smic_metrics_summary.csv
```




# Run count on a control sample (22901-001)

`nano count_22901-001_smic.job`

```bash
#BSUB -J count_22901-001_smic
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_22901-001_smic.out
#BSUB -e count_22901-001_smic.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=22901-001_Smic \
 --transcriptome=/nethome/kxw755/genomes/Smic_JGI/smic_mkref \
 --fastqs=/nethome/kxw755/20250306_AndradeRodriguez-22901-001 \
 --sample=AndradeRodriguez-22901-001_GEX3 \
 --create-bam=true
```

`bsub < count_22901-001_smic.job`

`nano count_22901-001_cgor.job`

```bash
#BSUB -J count_22901-001_cgor
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_22901-001_cgor.out
#BSUB -e count_22901-001_cgor.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=22901-001_Cgor \
 --transcriptome=/nethome/kxw755/genomes/Clago1/cgor_mkref \
 --fastqs=/nethome/kxw755/20250306_AndradeRodriguez-22901-001 \
 --sample=AndradeRodriguez-22901-001_GEX3 \
 --create-bam=true
```

`bsub < count_22901-001_cgor.job`


`nano count_22901-001_dtre.job`

```bash
#BSUB -J count_22901-001_dtre
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_22901-001_dtre.out
#BSUB -e count_22901-001_dtre.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=22901-001_Dtre \
 --transcriptome=/nethome/kxw755/genomes/Dtrenchii_SCF082/dtre_mkref \
 --fastqs=/nethome/kxw755/20250306_AndradeRodriguez-22901-001 \
 --sample=AndradeRodriguez-22901-001_GEX3 \
 --create-bam=true
```

`bsub < count_22901-001_dtre.job`

# Export files

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-22901-001/22901-001_Cgor/outs/web_summary.html ./22901-001_Cgor_web_summary.html
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-22901-001/22901-001_Cgor/outs/metrics_summary.csv ./22901-001_Cgor_metrics_summary.csv

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-22901-001/22901-001_Dtre/outs/web_summary.html ./22901-001_Dtre_web_summary.html
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-22901-001/22901-001_Dtre/outs/metrics_summary.csv ./22901-001_Dtre_metrics_summary.csv

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-22901-001/22901-001_Smic/outs/web_summary.html ./22901-001_Smic_web_summary.html
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-22901-001/22901-001_Smic/outs/metrics_summary.csv ./22901-001_Smic_metrics_summary.csv
```

# Run count on a control sample (23362_001)

`nano count_23362-001_smic.job`

```bash
#BSUB -J count_23362-001_smic
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_23362-001_smic.out
#BSUB -e count_23362-001_smic.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=23362-001_Smic \
 --transcriptome=/nethome/kxw755/genomes/Smic_JGI/smic_mkref \
 --fastqs=/nethome/kxw755/20250306_AndradeRodriguez-23362_001 \
 --sample=AndradeRodriguez-23362-001_GEX3 \
 --create-bam=true
```

`bsub < count_23362-001_smic.job`


`nano count_23362-001_cgor.job`

```bash
#BSUB -J count_23362-001_cgor
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_23362-001_cgor.out
#BSUB -e count_23362-001_cgor.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=23362-001_Cgor \
 --transcriptome=/nethome/kxw755/genomes/Clago1/cgor_mkref \
 --fastqs=/nethome/kxw755/20250306_AndradeRodriguez-23362_001 \
 --sample=AndradeRodriguez-23362-001_GEX3 \
 --create-bam=true
```

`bsub < count_23362-001_cgor.job`


`nano count_23362-001_dtre.job` 

```bash
#BSUB -J count_23362-001_dtre
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_23362-001_dtre.out
#BSUB -e count_23362-001_dtre.err
#BSUB -B
#BSUB -N
###################################################################

module load cellranger/8.0.0

cellranger count \
 --id=23362-001_Dtre \
 --transcriptome=/nethome/kxw755/genomes/Dtrenchii_SCF082/dtre_mkref \
 --fastqs=/nethome/kxw755/20250306_AndradeRodriguez-23362_001 \
 --sample=AndradeRodriguez-23362-001_GEX3 \
 --create-bam=true
```

`bsub < count_23362-001_dtre.job`


```bash
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-23362_001/23362-001_Cgor/outs/web_summary.html ./23362-001_Cgor_web_summary.html
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-23362_001/23362-001_Cgor/outs/metrics_summary.csv ./23362-001_Cgor_metrics_summary.csv

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-23362_001/23362-001_Dtre/outs/web_summary.html ./23362-001_Dtre_web_summary.html
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-23362_001/23362-001_Dtre/outs/metrics_summary.csv ./23362-001_Dtre_metrics_summary.csv

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-23362_001/23362-001_Smic/outs/web_summary.html ./23362-001_Smic_web_summary.html
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/20250306_AndradeRodriguez-23362_001/23362-001_Smic/outs/metrics_summary.csv ./23362-001_Smic_metrics_summary.csv
```

# Summary

![](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/20250407_Gfas_Symbiont_Mapping_test.png)
