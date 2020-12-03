---
layout: post
title: 20201203 Tag-Seq prep for Porites July Bleaching Experiment
date: '2020-12-03'
categories: Processing, Protocols
tags: RNA, Porites astreoides, Porites July Bleaching, Tag-seq
---

**Project:** [Porites Bleaching 2019](https://github.com/kevinhwong1/Porites_Rim_Bleaching_2019)

### Goal

Dilute extracted RNA to 10 ug/uL in 50 uL and consolidate into one plate for Tag-seq. Molecular QC information can be found [here](https://docs.google.com/spreadsheets/d/1bLsWHy7dzJcp06hSLESLgp66e4YscQf4IkLxQdHO2Ss/edit#gid=1053192266).

### Samples

| Fragment ID 	| Timepoint (Bag) 	| Group     	| Vial # 	| Extraction.Date 	| RNA Conc (ng/uL) 	| Sample Input for 10ng/ul (uL) 	| Diluent Vol (uL)  	| Total Vol (uL) 	| Well 	|
|-------------	|-----------------	|-----------	|--------	|-----------------	|------------------	|-------------------------------	|-------------------	|----------------	|------	|
| R29         	| A1 (MOLEC)      	| Bleached  	| 44     	| 20201117        	| 13.2             	| 37.88                         	| 12.12             	| 50             	| A1   	|
| R20         	| A4 (Whole frag) 	| Bleached  	| 23     	| 20201106        	| 15.6             	| 32.05                         	| 17.95             	| 50             	| A2   	|
| R19         	| A4 (Whole frag) 	| Bleached  	| 49     	| 20201119        	| 21.8             	| 22.94                         	| 27.06             	| 50             	| A3   	|
| R7          	| A2 (MOLEC)      	| Control   	| 37     	| 20201117        	| 24               	| 20.83                         	| 29.17             	| 50             	| A4   	|
| R7          	| A4 (Whole frag) 	| Control   	| 34     	| 20201111        	| 25.4             	| 19.69                         	| 30.31             	| 50             	| A5   	|
| R40         	| A1 (MOLEC)      	| Control   	| 88     	| 20201202        	| 33.8             	| 14.79                         	| 35.21             	| 50             	| A6   	|
| R37         	| A4 (Whole frag) 	| Bleached  	| 74     	| 20201125        	| 19.2             	| 26.04                         	| 23.96             	| 50             	| A7   	|
| R26         	| A4 (Whole frag) 	| Mortality 	| 87     	| 20201202        	| 13.2             	| 37.88                         	| 12.12             	| 50             	| A8   	|
| R29         	| A4 (Whole frag) 	| Bleached  	| 65     	| 20201120        	| 15.6             	| 32.05                         	| 17.95             	| 50             	| A9   	|
| R26         	| A1 (MOLEC)      	| Mortality 	| 46     	| 20201117        	| 16.4             	| 30.49                         	| 19.51             	| 50             	| A10  	|
| R40         	| A2 (MOLEC)      	| Control   	| 14     	| 20201029        	| 18.4             	| 27.17                         	| 22.83             	| 50             	| A11  	|
| R19         	| A2 (MOLEC)      	| Bleached  	| 39     	| 20201117        	| 14.8             	| 33.78                         	| 16.22             	| 50             	| A12  	|
| R28         	| A4 (Whole frag) 	| Mortality 	| 24     	| 20201106        	| 28               	| 17.86                         	| 32.14             	| 50             	| B1   	|
| R35         	| A2 (MOLEC)      	| Mortality 	| 13     	| 20201029        	| 15.2             	| 32.89                         	| 17.11             	| 50             	| B2   	|
| R28         	| A2 (MOLEC)      	| Mortality 	| 36     	| 20201111        	| 10.2             	| 49.02                         	| 0.98              	| 50             	| B3   	|
| R36         	| A4 (Whole frag) 	| Mortality 	| 40     	| 20201117        	| 14.2             	| 35.21                         	| 14.79             	| 50             	| B4   	|
| R20         	| A1 (MOLEC)      	| Bleached  	| 30     	| 20201110        	| 25.6             	| 19.53                         	| 30.47             	| 50             	| B5   	|
| R17         	| A1 (MOLEC)      	| Control   	| 15     	| 20201029        	| 18.8             	| 26.6                          	| 23.4              	| 50             	| B6   	|
| R40         	| A4 (Whole frag) 	| Control   	| 80     	| 20201126        	| 20.4             	| 24.51                         	| 25.49             	| 50             	| B7   	|
| R32         	| A4 (Whole frag) 	| Control   	| 33     	| 20201111        	| 13.4             	| 37.31                         	| 12.69             	| 50             	| B8   	|
| R37         	| A2 (MOLEC)      	| Bleached  	| 53     	| 20201119        	| 12.6             	| 39.68                         	| 10.32             	| 50             	| B9   	|
| R28         	| A1 (MOLEC)      	| Mortality 	| 71     	| 20201125        	| 20.6             	| 24.27                         	| 25.73             	| 50             	| B10  	|
| R11         	| A2 (MOLEC)      	| Mortality 	| 22     	| 20201106        	| 31.4             	| 15.92                         	| 34.08             	| 50             	| B11  	|
| R8          	| A2 (MOLEC)      	| Bleached  	| 58     	| 20201120        	| 14.2             	| 35.21                         	| 14.79             	| 50             	| B12  	|
| R17         	| A2 (MOLEC)      	| Control   	| 51     	| 20201119        	| 14               	| 35.71                         	| 14.29             	| 50             	| C1   	|
| R11         	| A4 (Whole frag) 	| Mortality 	| 63     	| 20201120        	| 20.8             	| 24.04                         	| 25.96             	| 50             	| C2   	|
| R11         	| A1 (MOLEC)      	| Mortality 	| 32     	| 20201111        	| 38.8             	| 12.89                         	| 37.11             	| 50             	| C3   	|
| R23         	| A2 (MOLEC)      	| Control   	| 11     	| 20201027        	| 15.8             	| 31.65                         	| 18.35             	| 50             	| C4   	|
| R7          	| A1 (MOLEC)      	| Control   	| 56     	| 20201119        	| 25.2             	| 19.84                         	| 30.16             	| 50             	| C5   	|
| R35         	| A4 (Whole frag) 	| Mortality 	| 41     	| 20201117        	| 10.8             	| 46.3                          	| 3.7               	| 50             	| C6   	|
| R23         	| A4 (Whole frag) 	| Control   	| 73     	| 20201125        	| 13.4             	| 37.31                         	| 12.69             	| 50             	| C7   	|
| R32         	| A2 (MOLEC)      	| Control   	| 35     	| 20201111        	| 17.8             	| 28.09                         	| 21.91             	| 50             	| C8   	|
| R8          	| A4 (Whole frag) 	| Bleached  	| 43     	| 20201117        	| 25               	| 20                            	| 30                	| 50             	| C9   	|
| R23         	| A1 (MOLEC)      	| Control   	| 67     	| 20201125        	| 15               	| 33.33                         	| 16.67             	| 50             	| C10  	|
| R32         	| A1 (MOLEC)      	| Control   	| 54     	| 20201119        	| 27               	| 18.52                         	| 31.48             	| 50             	| C11  	|
| R29         	| A2 (MOLEC)      	| Bleached  	| 21     	| 20201106        	| 29.2             	| 17.12                         	| 32.88             	| 50             	| C12  	|
| R17         	| A4 (Whole frag) 	| Control   	| 89     	| 20201202        	| 25.4             	| 19.69                         	| 30.31             	| 50             	| D1   	|
| R37         	| A1 (MOLEC)      	| Bleached  	| 10     	| 20201027        	| 24.6             	| 20.33                         	| 29.67             	| 50             	| D2   	|
| R35         	| A1 (MOLEC)      	| Mortality 	| 68     	| 20201125        	| 15.8             	| 31.65                         	| 18.35             	| 50             	| D3   	|
| R8          	| A1 (MOLEC)      	| Bleached  	| 31     	| 20201111        	| 13               	| 38.46                         	| 11.54             	| 50             	| D4   	|
| R20         	| A2 (MOLEC)      	| Bleached  	| 12     	| 20201027        	| 27.4             	| 18.25                         	| 31.75             	| 50             	| D5   	|
| R26         	| A2 (MOLEC)      	| Mortality 	| 57     	| 20201120        	| 23.2             	| 21.55                         	| 28.45             	| 50             	| D6   	|
| R36         	| A2 (MOLEC)      	| Mortality 	| 62     	| 20201120        	| 12.4             	| 40.32                         	| 9.68              	| 50             	| D7   	|
| R36         	| A1 (MOLEC)      	| Mortality 	| 55     	| 20201119        	| 16.2             	| 30.86                         	| 19.14             	| 50             	| D8   	|
| R19         	| A1 (MOLEC)      	| Bleached  	| 29     	| 20201110        	| 36.4             	| 13.74                         	| 36.26             	| 50             	| D9   	|


### Plate Map

![Plate Map]({{ site.baseurl}}/images/20201203_PlateMap.png "Plate Map")

### Protocol

1) Remove samples from the -80&deg;C freezer and place on a tube rack in an ice bucket to thaw.

2) Fill each well of the 96-well PCR plate with the corresponding diluent volume (I used Ultra Pure Water).

3) Place plate onto a block frozen at -80&deg;C.

4) Once samples are thawed, vortex and spin down the sample before adding the appropriate volume into the corresponding well.

5) Seal plate with aluminum foil tape and seal each well. Write on the label and the side of the plate with the following information:
  - Date
  - Name
  - Project
  - Job Number

![Plate label]({{ site.baseurl}}/images/20201203_PlateLabel.jpg "Plate label")

6) Place remaining samples and plate back into the -80&deg;C freezer.
