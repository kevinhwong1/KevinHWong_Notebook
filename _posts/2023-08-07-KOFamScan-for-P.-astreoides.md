---
layout: post
title: KOFamScan for P. astreoides
date: '2023-08-07'
categories: Analysis
tags: KOFamScan, Porites astreoides
---

# KOFamScan for *P. astreoides*

`cp ../Past_Genome/past_struc_annotations_v1/Pastreoides_proteins_v1.fasta ./`  

`nano kofamscan_Past.sh`

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
-o Past_KO_annot.txt \
-k ./ko_list \
-p ./profiles/eukaryote.hal \
-E 0.00001 \
-f detail-tsv \
--report-unannotated ./Pastreoides_proteins_v1.fasta

echo "Analysis complete!" $(date)
```

`sbatch kofamscan_Past.sh`

Took 1 hour to complete

`scp -r kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/kofamscan/Past_KO_annot.txt /Users/kevinwong/MyProjects/Porites_Rim_Bleaching_2019/data/Molecular/Past_KO_annot.txt`