---
layout: post
title: Genome Annotation on Galaxea fasicularis
date: '2023-09-11'
categories: Processing, Analysis
tags: KoFamScan, Galaxea
---

# Map KEGG terms to genome

This post is inspired by [E. Chille's](https://github.com/echille/E.-Chille-Open-Lab-Notebook/blob/master/_posts/2020-10-08-M-capitata-functional-annotation-pipeline.md), [D. Becker's](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-06-KofamScan-Workflow.md) and my [post](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/_posts/2023-03-23-KofamScan-Workflow.md) on genomic functional annotations. 

### 1. Download and inflate the Kofam database

Already completed in [this post](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/_posts/2023-03-23-KofamScan-Workflow.md)


### 2. Run KofamScan

Download Gfas protein file:

```bash
wget http://gfas.reefgenomics.org/download/gfas_1.0.proteins.fasta.gz
gunzip gfas_1.0.proteins.fasta.gz 
```

`nano kofamscan_gfas.sh`

```bash
#!/bin/bash
#SBATCH --job-name="KofamScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mem=100GB
#SBATCH -D /data/putnamlab/kevin_wong1/kofamscan

echo "Loading modules" $(date)
module load kofam_scan/1.3.0-foss-2019b
module load libyaml/0.1.5
module unload HMMER/3.3.1-foss-2019b
module load HMMER/3.3.2-gompi-2019b
module list

#echo "Starting analysis... downloading KO database" $(date)
#wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz #download KO database
#wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz
#gunzip ko_list.gz
#tar xf profiles.tar.gz

echo "Beginning mapping" $(date)
/opt/software/kofam_scan/1.3.0-foss-2019b/exec_annotation \
-o Gfas_KO_annot.txt \
-k ./ko_list \
-p ./profiles/eukaryote.hal \
-E 0.00001 \
-f detail-tsv \
--report-unannotated ./gfas_1.0.proteins.fasta

echo "Analysis complete!" $(date)
```

```bash
scp -r kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/kofamscan/Gfas_KO_annot.txt /Users/kevinwong/MyProjects/DarkGenes_Bleaching_Comparison/output/Annotations/Gfas_KO_annot.txt
```


# InterProScan

Resources:
* [User Manual](https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html)
* [My annotation pipeline on Porites astreoides](https://github.com/hputnam/Past_Genome/blob/master/genome_annotation_pipeline.md#12-interproscan)


`nano interproscan_gfas.sh`

```bash
#!/bin/bash
#SBATCH --job-name="InterProScan"
#SBATCH -t 30-00:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH --mem=100GB
#SBATCH --error="interproscan_out_error"
#SBATCH --output="interproscan_out"
#SBATCH -D /data/putnamlab/kevin_wong1/interproscan

echo "START $(date)"

# Load module
module load InterProScan/5.60-92.0-foss-2021b
module load Java/11.0.2
java -version

interproscan.sh --cpu $SLURM_CPUS_ON_NODE ...
interproscan.sh -version
interproscan.sh -f TSV -i gfas_1.0.proteins.fasta -b Gfas.interpro.20231004 -iprlookup -goterms -pa
#interproscan.sh -mode convert -f GFF3 -i Past.interpro.20220113.xml -b Past.interpro.20220113

# -i is the input data
# -b is the output file base
# -f is formats
# -iprlookup enables mapping
# -goterms is GO Term
# -pa is pathway mapping
# -version displays version number

echo "DONE $(date)"

```


# eggNOG-Mapper

http://eggnog-mapper.embl.de/

used the metazoan and all database and default parameters. 


