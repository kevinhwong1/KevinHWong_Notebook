---
layout: post
title: 20210401 WGBS PicoMethyl-seq library prep for Thermal Transplant Trial 1
date: '2021-04-01'
categories: Processing, Protocols
tags: Porites astreoides, WGBS, PicoMethyl-seq
---

**Project:** Thermal Transplant Molecular

### Goal

First trial of WGBS library prep using the Zymo Pico Methyl-seq library prep kit and [this protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resourses/DNA_RNA-protocols/WGBS-PicoMethyl-Protocol.md).


Sample dilution and mix calculations can be found on [this google sheet](https://docs.google.com/spreadsheets/d/1kthTxfiwn0RAWAQLLW3-pWBg5MBleQaFaEdjgEvHr58/edit#gid=0)


### Protocol

I am using a modified version of [this protocol](https://github.com/Putnam-Lab/Lab_Management/blob/master/Lab_Resourses/DNA_RNA-protocols/WGBS-PicoMethyl-Protocol.md).

The only modification I made to this protocol was changing the input DNA from 1ng to 10ng. This worked in previous library preps performed by E. Strand and M. Schedl.


#### Samples and Indices

| Sample | Coral.ID | Life Stage | i5 index # | i7 index # |
|:------:|:--------:|:----------:|:----------:|:----------:|
| 18-358 |  P-10-A  |    ADULT   |      1     |      1     |
|  18-20 |   R-7-B  |    ADULT   |      2     |      2     |
|  L-933 |  P-10-A  |   LARVAE   |      3     |      3     |


#### Qubit Results
To test DNA quantity: [Qubit](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2019-05-31-Qubit-Protocol.md)  

I used the DNA broad range assay.

|    Sample   | DNA 1 (ng/ul) | DNA 2 (ng/ul) |
|:-----------:|:-------------:|:-------------:|
|  Standard 1 |     184.99    |               |
|  Standard 2 |    19664.63   |               |
| 18-358 (#1) |      21.6     |      21.0     |
|  18-20 (#2) |      12.2     |      11.8     |
|  L-933 (#3) |      5.80     |      5.64     |


#### Tapestation Results

[D5000 Tapestation protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/DNA-Tapestation/)

[Full Tapestation Report](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/Tapestation_Results/2021-04-01_tapestation.pdf)

KW #1 (18-358):
![18-358]({{ site.baseurl}}/images/20210401_18-358.png "Summary")

KW #2 (18-20):
![18-358]({{ site.baseurl}}/images/20210401_18-20.png "Summary")

KW #3 (L-933):
![L-933]({{ site.baseurl}}/images/20210401_L-933.png "Summary")


#### Summary

Overall, this protocol worked, however there is PCR over-amplification. We will need to optimize the number number of PCR cycles to reduce this over-amplification.
