---
layout: post
title: Past genome MAKER rnd 1 troubleshooting
date: '2021-08-19'
categories: Analysis, Processing
tags: genome, annotation, porites astreoides, MAKER
---

# 20210817

In the desired folder for your MAKER run type:

`maker -CTL`

This will create three files:

- **maker_bopts.ctl** containing settings for BLAST and Exonerate.
- **masker_exe.ctl** with all the paths to different executables used by MAKER on your system.
- **maker_opts.ctl** is the file controlling MAKERs running behavior.

These files contain paths to files that are used by `maker` and our choices for analysis. Edit the maker_opts.ctl file to specify the genome assembly sequence, experimental alignment evidence and which gene finding method to use.

#### Modifications to `maker_opts.ctl`

For MAKER to run, modify the following with the appropriate paths:
- genome=**PATH_TO_GENOME**
- est=**PATH_TO_TRANSCRIPTOME**
- protein=**PATH_TO_PROTEIN**
- maxdnalength= (optional, default is 100000)

```
#-----Genome (these are always required)
genome=/data/putnamlab/kevin_wong1/Past_Genome/past_filtered_assembly.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=/data/putnamlab/kevin_wong1/Past_Genome/refs/Kenkel2013_past_transcriptome.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=/data/putnamlab/kevin_wong1/Past_Genome/refs/plut2v1.1.proteins.fasta  #protein sequence file in fasta format (i.e. from mutiple or$
protein_gff=  #aligned protein homology evidence from an external GFF3 file

maxdnalength=300000 #previously 1000000
```

#### Shell script: maker_rnd1.sh

```
#!/bin/bash
#SBATCH --job-name="MAKER_RND1"
#SBATCH -t 100:00:00
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome
#SBATCH --mem=100GB

module load maker/3.01.03

maker -base Rnd1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

echo "Mission complete." $(date)

```

JobID: 75609

Started: August 17, 2021 at 12pm

Ended: August 19, 2021 - CANCELLED TO MODIFY TIME/MEMORY/CPU PARAMETERS


# 20210819


From Kevin Bryan:

- I don’t know that I can gauge how long it should take, however I did notice that it is not running as efficiently as it could.
- First, remove the --ntasks-per-node=20, as it is not used here. This parameter is intended for MPI jobs.
- Then, it looks like maker takes a -cpus parameter that you did not include in your script. You probably should cancel it and resubmit with:

`maker -cpus $SLURM_CPUS_ON_NODE -base Rnd1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl`

- As most nodes have 36 cores, this should speed up the calculations.
- The help also has an -again parameter, which you should not specify, because it looks like it will not recalculate outputs that are already done by default, which is great.


#### Shell script: maker_rnd1.sh

```
#!/bin/bash
#SBATCH --job-name="MAKER_RND1"
#SBATCH -t 120:00:00
#SBATCH --export=NONE
#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/Past_Genome
#SBATCH --mem=250GB

module load maker/3.01.03

maker -cpus $SLURM_CPUS_ON_NODE -base Rnd1 maker_opts.ctl maker_bopts.ctl maker_exe.ctl

echo "Mission complete." $(date)

```

JobID: 77090

Started: August 19, 2021 at 11:24am

Ended:
