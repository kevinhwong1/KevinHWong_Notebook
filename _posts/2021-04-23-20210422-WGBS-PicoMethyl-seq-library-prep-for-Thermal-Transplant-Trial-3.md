---
layout: post
title: 20210422 WGBS PicoMethyl-seq library prep for Thermal Transplant Trial 3
date: '2021-04-23'
categories: Processing
tags: Porites astreoides, WGBS, PicoMethyl-seq
---

**Project:** Thermal Transplant Molecular

### Goal

Third trial of WGBS library prep using the Zymo Pico Methyl-seq library prep kit and [this protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/).

Today, I wanted to test 2 different modifications during the second amplification:

*Modification #1:*
- Use the **combined i5 and i7 index primers (an input of 1uL)**
- Reduce the Library Amp Master Mix (2x) to **13uL** (was originally 14uL)

*Modification #2:*
- Use the **KAPA HiFi HotStart ReadyMix (2X)** instead of the Library Amp Master Mix (2x)

Sample dilution and mix calculations can be found on [this google sheet](https://docs.google.com/spreadsheets/d/1kthTxfiwn0RAWAQLLW3-pWBg5MBleQaFaEdjgEvHr58/edit#gid=0)

### Protocol

I am using a modified version of [this protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/).

#### Samples and Indices

| Sample | Coral ID | Life Stage | i5 Index # | i7 Index # | Modification # |
|:------:|:--------:|:----------:|:----------:|:----------:|:--------------:|
| 18-118 |   P-6-A  |    Adult   |      7     |      7     |        1       |
| 18-106 |  P-19-A  |    Adult   |      8     |      8     |        1       |
| 18-334 |  P-19-B  |    Adult   |      9     |      9     |        2       |
| 18-442 |   R-9-B  |    Adult   |     10     |     10     |        2       |
| 18-167 |  R-19-A  |    Adult   |     11     |     11     |       NA       |
| L-1053 |  P-12-A  |   Larvae   |     12     |     12     |       NA       |
| L-1093 |  R-19-A  |   Larvae   |     13     |     13     |       NA       |
|  L-704 |  P-15-A  |   Larvae   |     14     |     14     |       NA       |

#### Qubit Results
To test DNA quantity: [Qubit](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2019-05-31-Qubit-Protocol.md)  

I used the DNA broad range assay.

I ran the qubit and tapestation on the following day (20210423). After the bead clean up step I stored the libraries in the fridge at 4C for ~16 hours.

|    Sample    | DNA 1 (ng/ul) | DNA 2 (ng/ul) |
|:------------:|:-------------:|:-------------:|
|  Standard 1  |     176.98    |               |
|  Standard 2  |    19057.22   |               |
|  18-118 (#7) |      25.2     |      24.8     |
|  18-106 (#8) |      25.2     |      24.8     |
|  18-334 (#9) |      28.0     |      27.8     |
| 18-442 (#10) |      73.2     |      72.0     |
| 18-167 (#11) |      29.4     |      28.8     |
| L-1053 (#12) |      8.90     |      8.50     |
| L-1093 (#13) |      24.6     |      24.0     |
|  L-704 (#14) |      7.08     |      6.92     |


#### Tapestation Results

[D5000 Tapestation protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/DNA-Tapestation/)

[Full Tapestation Report](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/Tapestation_Results/2021-04-23_tapestation.pdf)


![#7]({{ site.baseurl}}/images/20210423_7.png "Summary")

![#8]({{ site.baseurl}}/images/20210423_8.png "Summary")

![#9]({{ site.baseurl}}/images/20210423_9.png "Summary")

![#10]({{ site.baseurl}}/images/20210423_10.png "Summary")

![#11]({{ site.baseurl}}/images/20210423_11.png "Summary")

![#12]({{ site.baseurl}}/images/20210423_12.png "Summary")

![#13]({{ site.baseurl}}/images/20210423_13.png "Summary")

![#14]({{ site.baseurl}}/images/20210423_14.png "Summary")


#### Summary

Both modifications produced good libraries. The KAPA HiFi modification increased the yield. This also occurred in previous library preps.
