---
layout: post
title: 20210415 WGBS PicoMethyl-seq library prep for Thermal Transplant Trial 2
date: '2021-04-15'
categories: Processing, Protocols
tags: Porites astreoides, WGBS, PicoMethyl-seq
---


**Project:** Thermal Transplant Molecular

### Goal

Second trial of WGBS library prep using the Zymo Pico Methyl-seq library prep kit and [this protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resourses/DNA_RNA-protocols/WGBS-PicoMethyl-Protocol.md). The modification I made today was reducing the number of PCR cycles in the second amplification from 12 cycles to 10.


Sample dilution and mix calculations can be found on [this google sheet](https://docs.google.com/spreadsheets/d/1kthTxfiwn0RAWAQLLW3-pWBg5MBleQaFaEdjgEvHr58/edit#gid=0)


### Protocol

I am using a modified version of [this protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resourses/DNA_RNA-protocols/WGBS-PicoMethyl-Protocol.md).

The other modification I made to this protocol was changing the input DNA from 1ng to 10ng. This worked in previous library preps performed by E. Strand and M. Schedl.


#### Samples and Indices

| Sample | Coral.ID | Life Stage | i5 index # | i7 index # |
|:------:|:--------:|:----------:|:----------:|:----------:|
|  18-9  |   P-4-B  |    ADULT   |      4     |      4     |
|  18-91 |  R-14-B  |    ADULT   |      5     |      5     |
|  L-728 |   P-6-A  |   LARVAE   |      6     |      6     |


#### Qubit Results
To test DNA quantity: [Qubit](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2019-05-31-Qubit-Protocol.md)  

I used the DNA broad range assay.

|    Sample   | DNA 1 (ng/ul) | DNA 2 (ng/ul) |
|:-----------:|:-------------:|:-------------:|
|  Standard 1 |     182.09    |               |
|  Standard 2 |    19121.60   |               |
|  18-9 (#4)  |      17.1     |      16.8     |
|  18-91 (#5) |      15.8     |      15.3     |
|  L-728 (#6) |      13.4     |      13.2     |


#### Tapestation Results

[D5000 Tapestation protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/DNA-Tapestation/)

[Full Tapestation Report](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/Tapestation_Results/2021-04-15_tapestation.pdf)

KW #4 (18-9):
![18-9]({{ site.baseurl}}/images/20210415_18-9.png "Summary")

KW #5 (18-91):
![18-91]({{ site.baseurl}}/images/20210415_18-91.png "Summary")

KW #6 (L-728):
![L-728]({{ site.baseurl}}/images/20210415_L-728.png "Summary")


#### Summary

IT WORKED BEAUTIFULLY :)
