---
layout: post
title: 20210506 WGBS PicoMethyl-seq library prep for Thermal Transplant Trial 7
date: '2021-05-06'
categories: Processing
tags: Porites astreoides, WGBS, PicoMethyl-seq
---

**Project:** Thermal Transplant Molecular

### Goal

Seventh trial of WGBS library prep using the Zymo Pico Methyl-seq library prep kit and [this protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/). I finished the first prep of every sample and re-did one sample that previously failed (i.e. too low to read on the Qubit, L-862).

Sample dilution and mix calculations can be found on [this google sheet](https://docs.google.com/spreadsheets/d/1kthTxfiwn0RAWAQLLW3-pWBg5MBleQaFaEdjgEvHr58/edit#gid=0)

### Protocol

I am using a modified version of [this protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/).

For all samples, I used the following modifications:
- Use the **combined i5 and i7 index primers (an input of 1uL)**
- Reduce the Library Amp Master Mix (2x) to **13uL** (was originally 14uL)


#### Samples and Indices

| Sample | Coral ID | Life Stage | i5 Index # | i7 Index # |
|:------:|:--------:|:----------:|:----------:|:----------:|
| 18-346 |  P-6-B   |    Adult   |     41     |     41     |
| L-571  |  R-11-A  |    Larvae  |     42     |     42     |
| 18-250 |  R-11-B  |    Adult   |     43     |     43     |
| 18-418 |  R-9-A   |    Adult   |     44     |     44     |
| 18-454 |  P-12-A  |    Adult   |     45     |     45     |
|  18-44 |  R-19-B  |    Adult   |     46     |     46     |
| 18-466 |  R-15-B  |    Adult   |     47     |     47     |
| L-862  |  P-14-A  |    Larvae  |     48     |     48     |

**L-862** was the sample I re-did.


#### Qubit Results
To test DNA quantity: [Qubit](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2019-05-31-Qubit-Protocol.md)  

I used the DNA broad range assay.

|    Sample    | DNA 1 (ng/ul) | DNA 2 (ng/ul) |
|:------------:|:-------------:|:-------------:|
|  Standard 1  |     184.29    |               |
|  Standard 2  |    20541.44   |               |
|  41          |      7.28     |      6.94     |
|  42          |      10.1     |      10.2     |
|  43          |      8.22     |      8.18     |
|  44          |      21.2     |      21.0     |
|  45          |      5.42     |      5.38     |
|  46          |      21.6     |      21.4     |
|  47          |      19.1     |      19.0     |
|  48          |      24.4     |      24.2     |


#### Tapestation Results

[D5000 Tapestation protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/DNA-Tapestation/)

[Full Tapestation Report](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/Tapestation_Results/2021-05-06_tapestation.pdf)

![#41]({{ site.baseurl}}/images/20210506_41.png "Summary")

![#42]({{ site.baseurl}}/images/20210506_42.png "Summary")

![#43]({{ site.baseurl}}/images/20210506_43.png "Summary")

![#44]({{ site.baseurl}}/images/20210506_44.png "Summary")

![#45]({{ site.baseurl}}/images/20210506_45.png "Summary")

![#46]({{ site.baseurl}}/images/20210506_46.png "Summary")

![#47]({{ site.baseurl}}/images/20210506_47.png "Summary")

![#48]({{ site.baseurl}}/images/20210506_48.png "Summary")
