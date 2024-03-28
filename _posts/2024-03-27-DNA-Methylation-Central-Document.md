---
layout: post
title: DNA Methylation Central Document
date: '2024-03-27'
categories: Protocols
tags: Central Document
---

# DNA Methylation Sequencing Analysis Central Document

This will be my working document to connect all the datasets and different analysis pipelines together. 

## While Genome Bisulfite library prepration

Whole genome bisulfite sequencing is an approach to obtain single base-pair resolution of DNA methylation within your sample. I have developed a library prepration protocol with the Zymo Pico Methyl Seq kit for the coral *Porites astreoides* linked [here](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/).

## Bioinformatic pipelines

* [DNA methylation HPC workflow](https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/Past_WGBS_Workflow.md)
    * This workflow describes how to use NF-core Methylseq to quanitfy methylation from raw reads, filter for 5X or 10X coverage, and prepare data for further analyses in R. 
* DNA Methylation Characterization
    * This workflow describes how to characterize DNA methylation in various genomic features
    * The HPC workflow is linked [here](https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/Past_Genomic_Feature_Analysis_20221014.md)
    * The R workflow is linked [here]((https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/Methylation_Characterization.md))
* SNP extraction
    * 
* Plotting, Statsitics, WGCNA, Functional enrichment
    * This R script is linked [here](https://github.com/kevinhwong1/Thermal_Transplant_Molecular/blob/main/scripts/WGBS_GM.Rmd)