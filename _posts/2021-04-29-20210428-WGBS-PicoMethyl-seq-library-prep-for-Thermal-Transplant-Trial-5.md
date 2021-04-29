---
layout: post
title: 20210428 WGBS PicoMethyl-seq library prep for Thermal Transplant Trial 5
date: '2021-04-29'
categories: Processing
tags: Porites astreoides, WGBS, PicoMethyl-seq
---

**Project:** Thermal Transplant Molecular

### Goal

Fifth trial of WGBS library prep using the Zymo Pico Methyl-seq library prep kit and [this protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/).

Sample dilution and mix calculations can be found on [this google sheet](https://docs.google.com/spreadsheets/d/1kthTxfiwn0RAWAQLLW3-pWBg5MBleQaFaEdjgEvHr58/edit#gid=0)

### Protocol

I am using a modified version of [this protocol](https://kevinhwong1.github.io/KevinHWong_Notebook/Thermal-Transplant-WGBS-PicoMethyl-Protocol/).

For all samples, I used the following modifications:
- Use the **combined i5 and i7 index primers (an input of 1uL)**
- Reduce the Library Amp Master Mix (2x) to **13uL** (was originally 14uL)


#### Samples and Indices

| Sample | Coral ID | Life Stage | i5 Index # | i7 Index # |
|:------:|:--------:|:----------:|:----------:|:----------:|
| 18-406 |  R-15-A  |    Adult   |     23     |     23     |
|  18-32 |  P-12-B  |    Adult   |     24     |     24     |
| 18-262 |  P-15-B  |    Adult   |     25     |     25     |
| 18-322 |   P-9-A  |    Adult   |     26     |     26     |
|  18-79 |  R-14-A  |    Adult   |     27     |     27     |
|  L-661 |   R-9-A  |   Larvae   |     28     |     28     |
| L-1029 |   P-4-A  |   Larvae   |     29     |     29     |
|  L-862 |  P-14-A  |   Larvae   |     30     |     30     |
| L-1038 |   R-8-A  |   Larvae   |     31     |     31     |
| L-1257 |  P-19-A  |   Larvae   |     32     |     32     |

#### Qubit Results
To test DNA quantity: [Qubit](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2019-05-31-Qubit-Protocol.md)  

I used the DNA broad range assay.

I ran the qubit and tapestation on the following day (20210427). After the bead clean up step I stored the libraries in the fridge at 4C for ~16 hours.

|    Sample    | DNA 1 (ng/ul) | DNA 2 (ng/ul) |
|:------------:|:-------------:|:-------------:|
|  Standard 1  |     177.01    |               |
|  Standard 2  |    19683.34   |               |
|  23          |      24.0     |      23.6     |
|  24          |      13.1     |      12.5     |
|  25          |      15.2     |      14.8     |
|  26          |      3.68     |      3.40     |
|  27          |      5.64     |      5.24     |
|  28          |      3.16     |      2.56     |
|  29          |      7.84     |      7.54     |
|  30          |      LOW      |      LOW      |
|  31          |      14.2     |      13.9     |
|  32          |      5.78     |      5.36     |

#### Tapestation Results

[D5000 Tapestation protocol](https://meschedl.github.io/MESPutnam_Open_Lab_Notebook/DNA-Tapestation/)

[Full Tapestation Report](https://github.com/kevinhwong1/KevinHWong_Notebook/blob/master/images/Tapestation_Results/2021-04-28_tapestation.pdf)

![#23]({{ site.baseurl}}/images/20210428_23.png "Summary")

![#24]({{ site.baseurl}}/images/20210428_24.png "Summary")

![#25]({{ site.baseurl}}/images/20210428_25.png "Summary")

![#27]({{ site.baseurl}}/images/20210428_27.png "Summary")

![#29]({{ site.baseurl}}/images/20210428_29.png "Summary")

![#31]({{ site.baseurl}}/images/20210428_31.png "Summary")

![#32]({{ site.baseurl}}/images/20210428_32.png "Summary")
