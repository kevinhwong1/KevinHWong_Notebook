---
layout: post
title: Touchdown PCR for ITS2
date: '2022-03-16'
categories: Protocols
tags: ITS2, Porites astreoides
---

# Touchdown PCR for ITS2 processing

The goal of touchdown PCR ([Korbie and Mattick 2008](https://www.nature.com/articles/nprot.2008.133))is to increase specificity of the primers during PCR. I am currently having troubles with our [original ITS2 protocol](https://github.com/emmastrand/EmmaStrand_Notebook/blob/master/_posts/2020-01-31-ITS2-Sequencing-Protocol.md) on coral samples with very low Symbiodiniaceae densities. Our previous protocol is sequencing coral DNA when the Symbiodiniaceae densities are very low.

![]({{ site.baseurl}}/images/Relative_Abundance_Barplot.png)

**Previous gel results:**

G3 (Vial 43) is higher than the other samples.

![]({{ site.baseurl}}/images/20211109_ITS2_GelA.png)

Samples 74, 41, 40 had very faint bands.

![]({{ site.baseurl}}/images/20211116_ITS2_Gel.png)


# 20220315

This protocol I used from [Korbie and Mattick 2008](https://www.nature.com/articles/nprot.2008.133) and used a Tm of 52°C. In short, this protocol starts 10°C+Tm and decreases 1°C every cycle until it reaches Tm.

### Samples

I am using DNA that is already diluted to 4ng/uL in [this post](https://kevinhwong1.github.io/KevinHWong_Notebook/20211104-ITS2-Test-set-for-KW-AH-ES-samples/)

| Fragment_ID |       Group       | Vial | Well_dilution_plate |
|:-----------:|:-----------------:|:----:|:-------------------:|
|    R8-52    |      Bleached     |  43  |          G3         |
|    R37-52   |      Bleached     |  74  |          A6         |
|    R35-52   | Partial-Mortality |  41  |          F3         |
|    R36-52   | Partial-Mortality |  40  |          E3         |
|    R7-52    |      Control      |  34  |          H2         |


### Master Mix Calculations

**Determining number of wells:**


| 20220315 test set |    |
|:-----------------:|:--:|
|     # Samples     |  4 |
|  Positive Control |  1 |
|  Negative Control |  1 |
|       Error       |  1 |
|  Total # Samples  |  7 |
|   Total # Wells   | 14 |


**Calculating Master mix:**

|                    | Volume (uL) |  Total # Wells  | Master Mix Vol (uL) |
|:------------------:|:-----------:|:---------------:|:-------------------:|
| Phusion Master Mix |     12.5    |        14       |         175         |
|   F Primer (10µM)  |     0.5     |        14       |          7          |
|   R Primer (10µM)  |     0.5     |        14       |          7          |
|   Ultra Pure H2O   |     10.5    |        14       |         147         |
|                    |             |                 |                     |
|                    |             | Total Vol of MM |         336         |
|                    |             |    24 uL/Well   |          14         |

### PCR settings

| Phase | Step Number | Number of cycles | Step     | Temperature | Time  |
|-------|-------------|------------------|----------|-------------|-------|
| 1     | 1           | 1                | Denature | 95          | 3 min |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 62          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 61          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 60          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 59          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 58          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 57          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 56          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 55          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 54          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 53          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 52          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 2     | 1           | 25               | Denature | 95          | 30 s  |
| 2     | 2           | 25               | Anneal   | 52          | 45 s  |
| 2     | 3           | 25               | Elongate | 72          | 60 s  |
| 3     | 1           | 1                | Elongate | 72          | 5 min |
| 3     | 2           | 1                | Halt     | 4           | Hold  |


### Gel Results

Today I ran a 1% gel on a 100 well gel.

![Gel]({{ site.baseurl}}/images/20220315_ITS2_Gel.png "Gel")

This protocol did not seem to work based off of the gel results.

# 20220317

Today I ran 2 different variations of the thermocycler settings:

1. TD_2: Starting 10°C + Tm, going down a degree every cycle until Tm-5°C and continuing the remaining cycles at Tm-5°C.

2. TD_3: A protocol modified from H. Putnam, 7 min at 95°C, 20 cycles of 94°C for 30s, 62°C for 30s (decreased by 0.5°C at each cycle), 72°C for 1 min, followed by 12 additional cycles with an annealing temperature set at 52°C, final extension of 72°C for 7 min.

### TD_2 thermocycler settings

| Phase | Step Number | Number of cycles | Step     | Temperature | Time  |
|-------|-------------|------------------|----------|-------------|-------|
| 1     | 1           | 1                | Denature | 95          | 3 min |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 62          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 61          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 60          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 59          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 58          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 57          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 56          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 55          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 54          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 53          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 52          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 51          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 50          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 49          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 48          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 95          | 30 s  |
| 1     | 3           | 1                | Anneal   | 47          | 45 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 2     | 1           | 20               | Denature | 95          | 30 s  |
| 2     | 2           | 20               | Anneal   | 47          | 45 s  |
| 2     | 3           | 20               | Elongate | 72          | 60 s  |
| 3     | 1           | 1                | Elongate | 72          | 5 min |
| 3     | 2           | 1                | Halt     | 4           | Hold  |

### TD_3 thermocycler settings

| Phase | Step Number | Number of cycles | Step     | Temperature | Time  |
|-------|-------------|------------------|----------|-------------|-------|
| 1     | 1           | 1                | Denature | 95          | 7 min |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 62          | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 61.5        | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 60          | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 59.5        | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 59          | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 58.5        | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 58          | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 57.5        | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 57          | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 56.5        | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 56          | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 55.5        | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 55          | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 54.5        | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 54          | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 53.5        | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 53          | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 52.5        | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 1                | Denature | 94          | 30 s  |
| 1     | 3           | 1                | Anneal   | 52          | 30 s  |
| 1     | 4           | 1                | Elongate | 72          | 60 s  |
| 1     | 2           | 12               | Denature | 94          | 30 s  |
| 1     | 3           | 12               | Anneal   | 52          | 30 s  |
| 1     | 4           | 12               | Elongate | 72          | 60 s  |
| 3     | 1           | 1                | Elongate | 72          | 7 min |
| 3     | 2           | 1                | Halt     | 4           | Hold  |


### Samples Master Mix calculations

I used the same samples and master mix calculations for each thermocycler setting.

### Gel Results


![Gel]({{ site.baseurl}}/images/20220317_ITS2_Gel.png "Gel")
