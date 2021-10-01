---
layout: post
title: BUSCO on Pastreoides reference transcriptomes
date: '2021-09-30'
categories: Analysis
tags: Porites astreoides, transcriptome,BUSCO
---


# Creating new ref folders
mkdir Past_mansour
cd Past_mansour

# Uploading reference transcriptome and annotations from [Mansour](http://gigadb.org/dataset/100207)
wget https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100207/p_ast2016.fasta

wget https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100207/Por_ast_ann_allTrans.fasta

# BUSCO shell script for Mansour

```
#!/bin/bash

#SBATCH --job-name="busco"
#SBATCH --time="100:00:00"
#SBATCH --nodes 1 --ntasks-per-node=20
#SBATCH --mem=250G
##SBATCH --output="busco-%u-%x-%j"
##SBATCH --account=putnamlab
##SBATCH --export=NONE

echo "START" $(date)

labbase=/data/putnamlab
busco_shared="${labbase}/shared/busco"
[ -z "$query" ] && query="/data/putnamlab/kevin_wong1/REFS/Past_Mansour/p_ast2016.fasta" # set this to the query (genome/transcriptome) you are running
[ -z "$db_to_compare" ] && db_to_compare="${busco_shared}/downloads/lineages/metazoa_odb10"

source "${busco_shared}/scripts/busco_init.sh"  # sets up the modules required for this in the right order

# we require the agustus_config/ directory copied to a "writetable" location for
# busco to run and AUGUSTUS_CONFIG_PATH set to that

if [ ! -d "${labbase}/${USER}/agustus_config" ] ; then
    echo -e "Copying agustus_config/ to ${labbase}/${USER} .. "
    tar -C "${labbase}/${USER}" -xzf "${busco_shared}/agustus_config.tgz"
    echo done
fi

export AUGUSTUS_CONFIG_PATH="${labbase}/${USER}/agustus_config"
# This will generate output under your $HOME/busco_output
cd "${labbase}/${USER}"
busco --config "${busco_shared}/scripts/busco-config.ini"  -f -c 20 --long -i "${query}" -l "${db_to_compare}" -o busco_output_Past_Mansour -m transcriptome

echo "STOP" $(date)

```

# Running BUSCO for Mansour Transcriptome
```
sbatch -o ~/%u-%x.%j.out -e ~/%u-%x.%j.err --mail-type=BEGIN,END,FAIL --mail-user=kevin_wong1@uri.edu \
--export query=/data/putnamlab/kevin_wong1/REFS/Past_Mansour/p_ast2016.fasta  \
/data/putnamlab/kevin_wong1/scripts/run-busco-transcriptome_Past_mansour.sh
```

# Mansour Output
```
# BUSCO version is: 4.0.6
# The lineage dataset is: metazoa_odb10 (Creation date: 2019-11-20, number of species: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/kevin_wong1/REFS/Past_Mansour/p_ast2016.fasta
# BUSCO was run in mode: transcriptome

        ***** Results: *****

        C:30.5%[S:19.7%,D:10.8%],F:5.1%,M:64.4%,n:954
        291     Complete BUSCOs (C)
        188     Complete and single-copy BUSCOs (S)
        103     Complete and duplicated BUSCOs (D)
        49	Fragmented BUSCOs (F)
        614     Missing BUSCOs (M)
        954     Total BUSCO groups searched
```


# Running BUSCO for Kenkel

```
#!/bin/bash

#SBATCH --job-name="busco"
#SBATCH --time="100:00:00"
#SBATCH --nodes 1 --ntasks-per-node=20
#SBATCH --mem=250G
##SBATCH --output="busco-%u-%x-%j"
##SBATCH --account=putnamlab
##SBATCH --export=NONE

echo "START" $(date)

labbase=/data/putnamlab
busco_shared="${labbase}/shared/busco"
[ -z "$query" ] && query="/data/putnamlab/kevin_wong1/REFS/Past/Kenkel2013_past_transcriptome.fasta" # set this to the query (genome/transcriptome) you are running
[ -z "$db_to_compare" ] && db_to_compare="${busco_shared}/downloads/lineages/metazoa_odb10"

source "${busco_shared}/scripts/busco_init.sh"  # sets up the modules required for this in the right order

# we require the agustus_config/ directory copied to a "writetable" location for
# busco to run and AUGUSTUS_CONFIG_PATH set to that

if [ ! -d "${labbase}/${USER}/agustus_config" ] ; then
    echo -e "Copying agustus_config/ to ${labbase}/${USER} .. "
    tar -C "${labbase}/${USER}" -xzf "${busco_shared}/agustus_config.tgz"
    echo done
fi

export AUGUSTUS_CONFIG_PATH="${labbase}/${USER}/agustus_config"
# This will generate output under your $HOME/busco_output
cd "${labbase}/${USER}"
busco --config "${busco_shared}/scripts/busco-config.ini"  -f -c 20 --long -i "${query}" -l "${db_to_compare}" -o busco_output_Past_kenkel -m transcriptome

echo "STOP" $(date)
```

# Running BUSCO for Kenkel Transcriptome
```
sbatch -o ~/%u-%x.%j.out -e ~/%u-%x.%j.err --mail-type=BEGIN,END,FAIL --mail-user=kevin_wong1@uri.edu \
--export query=/data/putnamlab/kevin_wong1/REFS/Past/Kenkel2013_past_transcriptome.fasta  \
/data/putnamlab/kevin_wong1/scripts/run-busco-transcriptome_Past_Kenkel.sh
```

# Kenkel output

```
# BUSCO version is: 4.0.6
# The lineage dataset is: metazoa_odb10 (Creation date: 2019-11-20, number of species: 65, number of BUSCOs: 954)
# Summarized benchmarking in BUSCO notation for file /data/putnamlab/kevin_wong1/REFS/Past/Kenkel2013_past_transcriptome.fasta
# BUSCO was run in mode: transcriptome

        ***** Results: *****

        C:23.6%[S:22.7%,D:0.9%],F:30.9%,M:45.5%,n:954
        226     Complete BUSCOs (C)
        217     Complete and single-copy BUSCOs (S)
        9	Complete and duplicated BUSCOs (D)
        295     Fragmented BUSCOs (F)
        433     Missing BUSCOs (M)
        954     Total BUSCO groups searched
```
