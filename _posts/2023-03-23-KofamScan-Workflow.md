---
layout: post
title: KofamScan Workflow
date: '2023-03-23'
categories: Analysis
tags: KofamScan, bioinformatics
---

This post is inspired by [E. Chille's](https://github.com/echille/E.-Chille-Open-Lab-Notebook/blob/master/_posts/2020-10-08-M-capitata-functional-annotation-pipeline.md) and [D. Becker's](https://github.com/daniellembecker/DanielleBecker_Lab_Notebook/blob/master/_posts/2021-09-06-KofamScan-Workflow.md) posts on genomic functional annotations. 

# Map KEGG terms to genome

KofamScan is the command-line version of the popular KofamKOALA web-based tool, used to map Kegg terms (containing pathway information) to a genes. KofamScan and KofamKoala work by using HMMER/HMMSEARCH to search against KOfam (a customized HMM database of KEGG Orthologs (KOs). Mappings are considered robust because each Kegg term has an individual pre-defined threshold that a score has to exceed in order to map to a gene. While all mappings are outputted, high scoring (significant) assignments are highlighted with an asterisk.

The commands that I used are below. In order to run KofamScan, you will need a fasta file of predicted protein sequences (preferably the same one used to run InterProScan).

General Protocol:

### 1. Download and inflate the Kofam database

To get the most up-to-date Kofam database, download it just before running KofamScan. You will also need to download the profiles associated with the Kofam database containing threshold information.

```bash
cd /nethome/kxw755/opt
curl -O ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz  #download and unzip KO database
curl -O ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz  #download and inflate profiles

gunzip ko_list.gz
tar xf profiles.tar.gz
```
*Downloaded 20230323*


Install Anaconda on pegasus (only if you have not done this before): 

```bash
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh
source /nethome/kxw755/anaconda3/bin/activate
```
*Installed 20230323*

Install KofamScan:

```bash
cd /nethome/kxw755/opt/
conda install -c bioconda kofamscan
```
*something did not install correctly here. Need to revisit*

*I am going to excute this on andromeda instead*


```bash
mkdir kofamscan
cd ../../data/putnamlab/kevin_wong1/kofamscan

curl -O ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz  #download and unzip KO database
curl -O ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz  #download and inflate profiles

gunzip ko_list.gz
tar xf profiles.tar.gz
```
*Downloaded 20230323*



**Installed 20230323**

### 2. Run KofamScan

Download Mnemi protein file:

```bash
wget https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz
gunzip ML2.2.aa.gz 
```

`nano kofamscan_ML.sh`

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
-o Mnemi_KO_annot.txt \
-k ./ko_list \
-p ./profiles/eukaryote.hal \
-E 0.00001 \
-f detail-tsv \
--report-unannotated ./ML2.2.aa


echo "Analysis complete!" $(date)
```
