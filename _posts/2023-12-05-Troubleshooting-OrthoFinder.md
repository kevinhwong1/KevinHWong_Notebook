---
layout: post
title: Troubleshooting OrthoFinder
date: '2023-12-05'
categories: Analysis
tags: Bioinformatics
---

# Troubleshooing OrthoFinder to compare species orthologs

In this post I will be executing a few different combinations to compare orthologs between species: 

1. Mnemiopsis ledyi vs Drosophila melanogaster
* This will allow us to determine cell cycle phases based off Drosophila melanogaster marker genes

2. Galaxea fascicularis vs Stylophora pistillata vs Drosophila melanogaster
* This will allow us to compare dark genes across coral-related species

3. Xenia spp. vs Hvul vs Nematostella vs (maybe)

## References: 

OrthoFinder needs protein inputs in an fa format

* Mnemiopsis ledyi 
    * https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz
* Drosophila melanogaster
    *  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000803/UP000000803_7227.fasta.gz
* Galaxea fascicularis
* Stylophora pistillata
* Hvul
* Xenia spp.
* Nematostella

# Mnemi vs Dros

Make directories

```bash
mkdir orthofinder_comps
cd orthofinder_comps/
mkdir Mnemi_Dros
cd Mnemi_Dros/
```

Upload protein files

```bash
wget https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000000803/UP000000803_7227.fasta.gz


mv ML2.2.aa ML2.2.fasta
mv UP000000803_7227.fasta.gz dmel_prot.fasta.gz
```

https://ftp.flybase.net/releases/FB2023_05/dmel_r6.54/fasta/dmel-all-predicted-r6.54.fasta.gz


Unzip files

```bash
gunzip gunzip ML2.2.aa.gz 

gunzip dmel_prot.fasta.gz 
```

# Install OrthoFinder

```bash
module load anaconda3
source /share/apps/anaconda/anaconda3_build/bin/activate
conda create -n orthoENV -c bioconda orthofinder
```


`nano orthofinder_MD.job`

```bash
#BSUB -J orthofinder_MD
#BSUB -q general
#BSUB -P dark_genes
#BSUB -n 6
#BSUB -W 120:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o orthofinder_MD.out
#BSUB -e orthofinder_MD.err
#BSUB -B
#BSUB -N
###################################################################

orthofinder -f /scratch/projects/dark_genes/orthofinder_comps/Mnemi_Dros

```

```bash
module load anaconda3
source /share/apps/anaconda/anaconda3_build/bin/activate
conda activate orthoENV
```

`bsub < orthofinder_MD.job`







/nethome/kxw755/opt/OrthoFinder
/scratch/projects/dark_genes/modules/OrthoFinder
/nethome/kxw755

/projects/lsf_spool/1702051340.28323103.shell: line 13: OrthoFinder/orthofinder: No such file or directory

# Change bash profile

`nano .bash_profile`

```bash
# .bash_profile
# User specific environment and startup programs

PATH=$PATH:$HOME/bin
export PATH="$PATH:/nethome/kxw755/opt/cellranger-7.1.0"
export PATH="$PATH:/nethome/kxw755/opt/OrthoFinder"

export PATH
```



#### Running on andromeda


`nano of_MLDM.sh`

```bash
#!/bin/bash
#SBATCH --job-name="of_MLDM"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/of_comps

# load modules needed

module load OrthoFinder/2.5.2-intel-2019b-Python-3.7.4

orthofinder -f /data/putnamlab/kevin_wong1/of_comps/Mnemi_Dros

```


```bash
Writing orthogroups to file
---------------------------
OrthoFinder assigned 20340 genes (67.0% of total) to 5610 orthogroups. Fifty percent of all genes were in orthogroups with 2 or more genes (G50 was 2) and were cont
ained in the largest 3033 orthogroups (O50 was 3033). There were 3738 orthogroups with all species present and 2761 of these consisted entirely of single-copy genes
.
```


`$ scp -r kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/of_comps/Mnemi_Dros/OrthoFinder/Results_Dec10_1/Phylogenetic_Hierarchical_Orthogroups/N0.tsv /Users/kevinwong/MyProjects/Mnemi_Phagocyte/output/orthofiner_MLDM.tsv`


`$ scp -r kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/of_comps/Mnemi_Dros/ML2.2.fasta /Users/kevinwong/MyProjects/Mnemi_Phagocyte/output/ML2.2.fasta`

/data/putnamlab/kevin_wong1/of_comps/Mnemi_Dros/OrthoFinder/Results_Dec10_1/Phylogenetic_Hierarchical_Orthogroups




#### testing blastx 


`nano ml_dm_blastx.sh`

```bash
#!/bin/bash
#SBATCH --job-name="blastx"
#SBATCH -t 500:00:00
#SBATCH --nodes=1 --ntasks-per-node=10
#SBATCH --mem=120GB
#SBATCH --account=putnamlab
#SBATCH --export=NONE
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=kevin_wong1@uri.edu
#SBATCH -D /data/putnamlab/kevin_wong1/of_comps/Mnemi_Dros

#load module
module load BLAST+/2.11.0-gompi-2020b #load blast module

# Make Mnemi database
makeblastdb -in ML2.2.fasta -dbtype prot -out ML_db

#run blastx
blastx -query dmel-all-gene-r6.54.fasta -db ML_db -out ML_DM.out -outfmt 6 -evalue 1e-6

```

`$ scp -r kevin_wong1@ssh3.hac.uri.edu:/data/putnamlab/kevin_wong1/of_comps/Mnemi_Dros/ML_DM.out /Users/kevinwong/MyProjects/Mnemi_Phagocyte/output/MLDM_blastx.tsv`