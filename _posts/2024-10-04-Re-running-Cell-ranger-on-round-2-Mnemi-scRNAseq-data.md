---
layout: post
title: Re-running Cell ranger on round 2 Mnemi scRNAseq data
date: '2024-10-04'
categories: Analysis
tags: mnemi, scRNAseq
---

# Copy files to pegasus

scp -r 10X_round2/ kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/Mnemi_phagocyte_10x/

# Count 

### Run_2_AllCells_1

`nano count_R2_AllCells1_v3_cat.job`

```bash
#BSUB -J count_R2_AllCells1
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_R2_AllCells1.out
#BSUB -e count_R2_AllCells1.err
#BSUB -B
#BSUB -N
###################################################################

cellranger count \
 --id=R2_AllCells1_v3_cat \
 --transcriptome=/nethome/kxw755/Mnemi_phagocyte_10x/Mle_F31_T2T_genome/mnemi_v3_cat \
 --fastqs=/nethome/kxw755/Mnemi_phagocyte_10x/10X_round2/Round2_AllCells_1 \
 --sample=lib99405_GEX
```

`bsub < count_R2_AllCells1_v3_cat.job`


### Run_2_AllCells_2

`nano count_R2_AllCells2_v3_cat.job`

```bash
#BSUB -J count_R2_AllCells2
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_R2_AllCells2.out
#BSUB -e count_R2_AllCells2.err
#BSUB -B
#BSUB -N
###################################################################

cellranger count \
 --id=R2_AllCells2_v3_cat \
 --transcriptome=/nethome/kxw755/Mnemi_phagocyte_10x/Mle_F31_T2T_genome/mnemi_v3_cat \
 --fastqs=/nethome/kxw755/Mnemi_phagocyte_10x/10X_round2/Round2_AllCells_2 \
 --sample=lib99406_GEX
```

`bsub < count_R2_AllCells2_v3_cat.job`

### Run_2_AllCells_3

`nano count_R2_AllCells3_v3_cat.job`

```bash
#BSUB -J count_R2_AllCells3
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -u kxw755@earth.miami.edu
#BSUB -o count_R2_AllCells3.out
#BSUB -e count_R2_AllCells3.err
#BSUB -B
#BSUB -N
###################################################################

cellranger count \
 --id=R2_AllCells3_v3_cat \
 --transcriptome=/nethome/kxw755/Mnemi_phagocyte_10x/Mle_F31_T2T_genome/mnemi_v3_cat \
 --fastqs=/nethome/kxw755/Mnemi_phagocyte_10x/10X_round2/Round2_AllCells_3 \
 --sample=lib99407_GEX
```

`bsub < count_R2_AllCells3_v3_cat.job`


# Aggrgate

`mkdir R2_AllCells_Aggr`

```bash
cd /nethome/kxw755/Mnemi_phagocyte_10x/10X_round2/Round2_AllCells_1/R2_AllCells1_v3_cat/outs
cp molecule_info.h5 ../../../R2_AllCells_Aggr/R2_AllCells1_molecule_info.h5

cd /nethome/kxw755/Mnemi_phagocyte_10x/10X_round2/Round2_AllCells_2/R2_AllCells2_v3_cat/outs
cp molecule_info.h5 ../../../R2_AllCells_Aggr/R2_AllCells2_molecule_info.h5

cd /nethome/kxw755/Mnemi_phagocyte_10x/10X_round2/Round2_AllCells_3/R2_AllCells3_v3_cat/outs
cp molecule_info.h5 ../../../R2_AllCells_Aggr/R2_AllCells3_molecule_info.h5
```

`nano AllCells_libraries.csv`

```bash
sample_id,molecule_h5
R2-AllCells1,R2_AllCells1_molecule_info.h5
R2-AllCells2,R2_AllCells2_molecule_info.h5
R2-AllCells3,R2_AllCells3_molecule_info.h5
```

## Run aggr in CellRanger

`nano aggr_allcells.job`

```bash
#BSUB -J R2_AllCells_aggr
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 16
#BSUB -W 120:00
#BSUB -u kxw755@earth.miami.edu
#BSUB -o AllCells_aggr.out
#BSUB -e AllCells_aggr.err
#BSUB -B
#BSUB -N
###################################################################

cellranger aggr --id=R2_AllCells \
 --csv=/nethome/kxw755/Mnemi_phagocyte_10x/10X_round2/R2_AllCells_Aggr/AllCells_libraries.csv \
 --normalize=mapped
```

`bsub < aggr_allcells.job`

# Export

```bash
scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/Mnemi_phagocyte_10x/10X_round2/R2_AllCells_Aggr/R2_AllCells/outs/count/filtered_feature_bc_matrix.h5 R2_AllCells_ONLY_Aggr_filtered_feature_bc_matrix.h5

scp -r kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/Mnemi_phagocyte_10x/10X_round2/R2_AllCells_Aggr/R2_AllCells/outs/web_summary.html R2_AllCells_ONLY_Aggr_web_summary.html
```