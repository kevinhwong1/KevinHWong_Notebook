---
layout: post
title: Training AUGUSTUS with BUSCO
date: '2021-10-21'
categories: Analysis
tags: porites astreoides, genome, augustus, busco, bioinfiormatics
---

First attempt at training AUGUSTUS with BUSCO

References:
- [1](https://scilifelab.github.io/courses/annotation/2017/practical_session/TrainingAbInitionpredictor.html)
- [2](https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html)
- [3](https://www.biostars.org/p/207801/)

## Config parameters

```
# This is the BUSCOv5 configuration file template.
# It is not necessary to use this, as BUSCO will use the dependencies available on your PATH by default.
# The busco run parameters can all be set on the command line. See the help prompt (busco -h) for details.
#
# To use this file for an alternative configuration, or to specify particular versions of dependencies:
# 1) edit the path and command values to match your desired dependency versions.
#    WARNING: passing a parameter through the command line overrides the value specified in this file.
#
# 2) Enable a parameter by removing ";"
#
# 3) Make this config file available to BUSCO either by setting an environment variable
#
#                   export BUSCO_CONFIG_FILE="/path/to/myconfig.ini"
#
#    or by passing it as a command line argument
#
#                   busco <args> --config /path/to/config.ini
#
[busco_run]
# Input file
;in = /path/to/input_file.fna
# Run name, used in output files and folder
;out = BUSCO_run
# Where to store the output directory
;out_path = /path/to/output_folder
# Path to the BUSCO dataset
;lineage_dataset = /data/putnamlab/shared/busco/downloads
# Which mode to run (genome / proteins / transcriptome)
;mode = genome
# Run lineage auto selector
;auto-lineage = True
# Run auto selector only for non-eukaryote datasets
;auto-lineage-prok = True
# Run auto selector only for eukaryote datasets
;auto-lineage-euk = True
# How many threads to use for multithreaded steps
;cpu = 16
# Force rewrite if files already exist (True/False)
;force = False
# Restart a previous BUSCO run (True/False)
;restart = False
# Blast e-value
;evalue = 1e-3
# How many candidate regions (contigs, scaffolds) to consider for each BUSCO
;limit = 3
# Metaeuk parameters for initial run
;metaeuk_parameters='--param1=value1,--param2=value2'
# Metaeuk parameters for rerun
;metaeuk_rerun_parameters=""
# Augustus parameters
;augustus_parameters='--param1=value1,--param2=value2'
# Quiet mode (True/False)
;quiet = False
# Local destination path for downloaded lineage datasets
;download_path = /data/putnamlab/shared/busco/downloads/
# Run offline
;offline=True
# Ortho DB Datasets version
;datasets_version = odb10
# URL to BUSCO datasets
;download_base_url = https://busco-data.ezlab.org/v4/data/
# Download most recent BUSCO data and files
;update-data = True
# Use Augustus gene predictor instead of metaeuk
;use_augustus = True

[tblastn]
path = /opt/software/BLAST+/2.11.0-gompi-2020b/bin/
command = tblastn

[makeblastdb]
path = /opt/software/BLAST+/2.11.0-gompi-2020b/bin/
command = makeblastdb

[metaeuk]
path = /opt/software/MetaEuk/4-GCC-10.2.0/bin/
command = metaeuk

[augustus]
path = /opt/software/AUGUSTUS/3.4.0-foss-2020b/bin/
command = augustus

[etraining]
path = /opt/software/AUGUSTUS/3.4.0-foss-2020b/bin/
command = etraining

[gff2gbSmallDNA.pl]
path = /opt/software/AUGUSTUS/3.4.0-foss-2020b/scripts/
command = gff2gbSmallDNA.pl

[new_species.pl]
path = /opt/software/AUGUSTUS/3.4.0-foss-2020b/scripts/
command = new_species.pl

[optimize_augustus.pl]
path = /opt/software/AUGUSTUS/3.4.0-foss-2020b/scripts/
command = optimize_augustus.pl

[hmmsearch]
path = /opt/software/HMMER/3.3.2-gompi-2020b/bin/

[sepp]
path = /opt/software/SEPP/4.4.0-foss-2020b/bin
command = run_sepp.py

[prodigal]
path = /opt/software/prodigal/2.6.3-GCCcore-10.2.0/bin
command = prodigal

```


## script

`nano AUGUSTUS_BUSCO.sh`

```
#!/bin/bash
#SBATCH --job-name="BUSCO"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome/AUGUSTUS_BUSCO
#SBATCH --mem=500GB

echo "Starting BUSCO" $(date)

#load modules
module load BUSCO/5.2.2-foss-2020b

#run BUSCO
busco \
--config config.ini \
--in ../past_filtered_assembly.fasta \
--out past_geome_AUG_BUSCO \
-l busco_downloads/metazoa_odb10 \
-m genome \
-f \
--long \
--augustus_parameters='--progress=true' \
--offline

echo "BUSCO Mission complete!" $(date)
```

`sbatch /data/putnamlab/kevin_wong1/Past_Genome/BUSCO_AUGUSTUS_past.sh`
