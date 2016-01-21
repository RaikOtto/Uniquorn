---
title: "Overview of the Uniquorn package"
author: "Raik Otto"
date: "19. Januar 2016"
output: html_document
---


# Uniquorn R package

Package to identify cancer cell lines (CL)s based on their weighted mutational fingerprint.

# How to make it work: 

Start R

# 1 Preparation 

optional if you have the packages installed and loaded

`install.packages("devtools")`

`library("devtools")`

`source("https://bioconductor.org/biocLite.R")`

# 2 Install Uniquorn

`install_github("RaikOtto/Uniquorn")`

Note that some systems require the prior command `options(unzip = 'internal')` to install from github 

# 3 Test run

Here the NCI-60 exome sequenced HT29 Cancer Cell line, reference genome GRCh37/ HG19


```r
library("Uniquorn")
```


```r
HT29_vcf_file = system.file("extdata/HT29.vcf.gz", package="Uniquorn")
```















