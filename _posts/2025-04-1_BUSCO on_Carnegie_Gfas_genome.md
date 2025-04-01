---
layout: post
title: BUSCO on Carnegie Gfas genome
date: '2025-04-01'
categories: Analysis
tags: genome
---

# Transfer files to pegasus

```bash
scp -r carnegie_gfas_v1_pkg-version-0.0.tar.gz kxw755@pegasus.ccs.miami.edu:/nethome/kxw755/genomes/

gunzip /nethome/kxw755/genomes/carnegie_gfas_v1_pkg/02-associated/genome_assembly/canu_associated.fasta.gz
gunzip /nethome/kxw755/genomes/carnegie_gfas_v1_pkg/02-associated/gene_annotations/carnegie_gfas_associated_OGS_v1.0.gtf.gz
```

# Create and activate a new conda environment for BUSCO

```bash
conda create -n busco_env -c conda-forge -c bioconda busco=5.4.7 -y
conda activate busco_env
```
# Download BUSCO databases

```bash
mkdir BUSCO
cd BUSCO/
busco --download eukaryota_odb10
busco --download metazoa_odb10
```

# Run BUSCO

`nano busco.job `

```bash
#!/bin/bash
#BSUB -J BUSCO_Genome
#BSUB -q bigmem
#BSUB -P dark_genes
#BSUB -n 8
#BSUB -W 120:00
#BSUB -R "rusage[mem=15000]"
#BSUB -o busco_%J.out
#BSUB -e busco_%J.err
#BSUB -N

# Set paths
GENOME="/nethome/kxw755/genomes/carnegie_gfas_v1_pkg/02-associated/genome_assembly/canu_associated.fasta"
LINEAGE="/nethome/kxw755/BUSCO/busco_downloads/lineages/metazoa_odb10"  # e.g., eukaryota_odb10
OUT_DIR="busco_output"
RUN_NAME="02-associated_gfas_busco"

# Run BUSCO
busco -i "$GENOME" \
      -l "$LINEAGE" \
      -o "$RUN_NAME" \
      -m genome \
      --cpu 8 \
      --out_path "$OUT_DIR"
```

# Analyze output
```bash
        --------------------------------------------------
        |Results from dataset metazoa_odb10               |
        --------------------------------------------------
        |C:90.9%[S:87.8%,D:3.1%],F:1.8%,M:7.3%,n:954      |
        |868    Complete BUSCOs (C)                       |
        |838    Complete and single-copy BUSCOs (S)       |
        |30     Complete and duplicated BUSCOs (D)        |
        |17     Fragmented BUSCOs (F)                     |
        |69     Missing BUSCOs (M)                        |
        |954    Total BUSCO groups searched               |
        --------------------------------------------------
```

# Investigate 02-associated GTF for potential gene overlaps present in 01-primary genome assembly

```bash
associated_contig_1     AUGUSTUS        gene    531126  536859  .       +       .       carnegie_gfas-associated-1.0_g50
associated_contig_1     AUGUSTUS        gene    553887  587426  .       +       .       carnegie_gfas-associated-1.0_g51
associated_contig_1     AUGUSTUS        gene    591159  598460  .       +       .       carnegie_gfas-associated-1.0_g52
associated_contig_1     AUGUSTUS        gene    605382  611641  .       +       .       carnegie_gfas-associated-1.0_g53
associated_contig_1     AUGUSTUS        gene    617204  622463  .       +       .       carnegie_gfas-associated-1.0_g54

associated_contig_250   AUGUSTUS        gene    509285  518766  .       +       .       carnegie_gfas-associated-1.0_g19409
associated_contig_250   AUGUSTUS        gene    525504  526502  .       +       .       carnegie_gfas-associated-1.0_g19410
associated_contig_250   AUGUSTUS        gene    531386  534241  .       -       .       carnegie_gfas-associated-1.0_g19411
```
