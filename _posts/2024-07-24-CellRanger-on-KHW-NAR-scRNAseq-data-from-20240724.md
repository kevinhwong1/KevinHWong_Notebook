---
layout: post
title: CellRanger on KHW NAR scRNAseq data from 20240724
date: '2024-07-24'
categories: Analysis
tags: scRNAseq
---

# Sample information

| Sample # | Customer Sample Name |       BaseSpace Sample ID       | Flow Cell Type |    Read Type   | Total number of clusters (reads) per sample |
|:--------:|:--------------------:|:-------------------------------:|:--------------:|:--------------:|:-------------------------------------------:|
|     1    |     gfas_control     | AndradeRodriguez-19876-001_GEX3 |     10B-300    | paired-end,151 |                 581,421,834                 |
|     2    |      gfas_bleach     | AndradeRodriguez-19876-002_GEX3 |     10B-300    | paired-end,151 |                 636,439,053                 |


# Upload files
scp -r /Users/cnidarianimmunity/Desktop/KHW/Andrade_Rodriguez-01-07-19876-425120288/2024* kxw755@pegasus.ccs.miami.edu:/nethome/kxw755
