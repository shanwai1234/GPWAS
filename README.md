# GPWAS [![](https://img.shields.io/badge/Release-v1.0.1-blue.svg)](https://github.com/shanwai1234/GPWAS/commits/master) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
## Genome-Phenome Wide Association Study

### Authors:
> Zhikai Liang, Yumou Qiu and James C. Schnable.

### Contact:
> [zliang@huskers.unl.edu](Zhikai Liang)

---
# Installation and Package loading
```r
# install "devtools" package in your R environment
> devtools::install_github("shanwai1234/GPWAS")
> library(GPWAS)
```
# Data Preparation
## Genotype Data
| Gene | SNP | Sample1 | Sample2 | Sample3 |
| :---: | :---: |:---: |:---: | :---: |
|Zm00001d000501|S3_27396852| 0 | 1 | 0 |
|Zm00001d000501|S3_27397388| 0 | 2 | 1 |
|Zm00001d000632|S3_156234383| 0 | 1 | 1 |
|Zm00001d000632|S3_156234434| 1 | 2 | 1 |
|Zm00001d000654|S3_178548901| 2 | 1 | 1 |
## Phenotype Data
## PC covariates

