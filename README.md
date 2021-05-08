# GPWAS [![](https://img.shields.io/badge/Release-v1.0.1-blue.svg)](https://github.com/shanwai1234/GPWAS/commits/master) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
## Genome-Phenome Wide Association Study

<p align="center">
<img src="GPWAS-logo2.png" height="300px" width="200px">
</a>
</p>

### Authors:
> Zhikai Liang, Yumou Qiu and James C. Schnable.

### Contact:
> [zliang@huskers.unl.edu] or [liang795@umn.edu] (Zhikai Liang)
> Or raising your questions at ```issues```

### Citation:
> Liang, Z., Qiu, Y. and Schnable, JC (2020). Genome-phenome wide association in maize and Arabidopsis identifies a common molecular and evolutionary signature. [Molecular Plant](https://www.cell.com/molecular-plant/fulltext/S1674-2052(20)30065-4)
---
# Installation and Package loading
```r
# R version is required >= 3.4.4
# When the first time to use the package, please make sure MASS and leaps packages are installed under your R environment, if not, please use commands below to install
> install.packages("leaps")
> install.packages("MASS")
# install "devtools" package in your R environment
> devtools::install_github("shanwai1234/GPWAS")
> library(GPWAS)
```
# Data Preparation

All of input data are required to be organized in following format.

## Genotype Data

**Gene**: Keep the title of this item as "Gene" and do not change it. All of gene names should be kept below it.

**SNP**: Keep the title of this item as "SNP" and do not change it. The format of SNP should be written as "S"+chromosome+"\_"+SNP position.

**Sample**: Individual sample name.

**Note**: All items should be split by \tab. Making sure there is no missing data in your genotype file. For each SNP per gene, they should not be completely identical. We do recommend you to filter SNP based upon MAF (Minor Allele Frequency).

| Gene | SNP | Sample1 | Sample2 | Sample3 |
| :---: | :---: |:---: |:---: | :---: |
|Zm00001d000501|S3_27396852| 0 | 1 | 0 |
|Zm00001d000501|S3_27397388| 0 | 2 | 1 |
|Zm00001d000632|S3_156234383| 0 | 1 | 1 |
|Zm00001d000632|S3_156234434| 1 | 2 | 1 |
|Zm00001d000654|S3_178548901| 2 | 1 | 1 |

## Phenotype Data

**Pheno**: Name of phenotype name.

**Sample**: Individual sample name.

**Note**: All items should be split by space. Make sure there is no missing data in your phenotype file or any row containing missing data will be removed for following analysis. The number of analyzed phenotype is better not exceed the number of individuals in the studied population. If you have extremely high-dimensional phenotypes, it is suggested to reduce ones that are too similar with others.

| Taxa | Pheno1 | Pheno2 | Pheno3 |
| :---: | :---: |:---: |:---: |
| Sample1 | 23 | 50 | 55 |
| Sample2 | 10 | 12 | 8 |
| Sample3 | 120 | 150 | 133 |

## PC covariates
Version 1.0.1 GPWAS package controls population structure using PC scores generated by PCA. For each individual gene, population structure is controlled by PC scores calculated by the rest of other chromosomes. But generating PC scores were not included in GPWAS package. You are suggested to calculate it using function like [prcomp](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html) or [TASSEL](https://www.maizegenetics.net/tassel)

|  | "PC1" | "PC2" | "PC3" |
| :---: | :---: |:---: |:---: |
| "Sample1" | -11 | 8 | 31 |
| "Sample2" | -36 | 35 | 138 |
| "Sample3" | 15 | 10 | -7 |

You need to prepare separate PC covariate files excluding each individual chromosome. If you have 10 chromosomes, you need to prepare 10 separate PC covariate files.

**Note**: All items should be split by space. The order of samples should be identical to the order of samples in both genotype and phenotype file.

**Example**: When you want to exclude chromosome 1, you need to make the file name such as "exclude-chr1.txt". Then store all of these files to a folder.

## Imputation suggestion
GPWAS requires both genotype and phenotype matrice to be complete. That being said, there is no missing data points allowed. For genotype data, imputation software [Beagle](https://faculty.washington.edu/browning/beagle/beagle.html) could work for this purpose. For phenotype data, it depends on your missing data rate. If your missing data point is in a low rate, you could just use the mean value of other observed data points for a trait to replace missing data point for this trait. It just needs a simple code to implement it, an example could be found in [here](https://stackoverflow.com/questions/25835643/replace-missing-values-with-column-mean). If missing data point in your data is common, you could either consider use [Phenix](https://mathgen.stats.ox.ac.uk/genetics_software/phenix/phenix.html) or [softImpute](https://cran.r-project.org/web/packages/softImpute/index.html) to implement the phenotype impputation.  

# How to use
```r
> gpwas(ingeno, inpheno, inpc, gp, gv, R = num)
```
 **ingeno**: Input genotype file name/directory. It is recommended to split big genotype file into multiple in order to reduce memory load.

 **inpheno**: Input phenotype file name/directory.

 **inpc**: Input folder with PCA parsed population structure covariance. If n number of chromosomes, n number of separate files should be included, as SNPs on each chromosome is excluded for performing PCA once.

 **gp**: Output file name/directory for selected phenotypes with every gene as well as p value of each selected phenotypes (both gene names and values are just examples, PC[,1] shows the p value for the first principal component, and same for other selected principal components, Pheno1 and Pheno2 are incorporated phenotypes for this Gene1. P-values in this gp file are generally not considered. One the gene in gv file is considered as signficant or in high rank of p-value, you could check in gp file for that specific gene and see which traits were selected for it).

| "Gene" | "Predictor" | "Predictor p-value" |
| :---: | :---: |:---: |
| "Gene1" | PC[,1] | 1.0e-4 |
| "Gene1" | PC[,2] | 1.2e-2 |
| "Gene1" | PC[,3] | 4.3e-2 |
| "Gene1" | Pheno1 | 3.7e-4 |
| "Gene1" | Pheno2 | 1.4e-3 |

 **gv**: Output file name/directory of terminated p value for each gene (both gene names and values are just examples).

| "Gene" | "GPWAS p-value" |
| :---: | :---: |
| "Gene1" | 2.2e-6 |
| "Gene2" | 1.3e-8 |
| "Gene3" | 4.5e-10 |

 **R**: Number of iteration for scanning all of input phenotypes with one specific gene. Too big number will be redundancy and computationally cost. Suggested ranging from 10-50.

 **Note**: Once you continueously seeing "No new add-in" and "No leave-out" means the model is stable.

 ```r
 # Customizing more
 > gpwas(ingeno, inpheno, inpc, g, gp, gv, R = num, pc = 3, selectIn = 0.01, selectOut = 0.01)
 ```
**g**: A list of specific gene that needs to analyze. By default the model will run for all of genes detected in the input genotype file.

Example as below:

| Candidate |
| :---: |
| Gene1 |
| Gene2 |
| Gene3 |

**pc**: Number of PCs that needs to be included to control the population structure.

**selectIn**: p value threshold to determine if a phenotype could be selected in the model.

**selectOut**: p value threshold to determine if a phenotype could be dropped out from the model.

```r
# Run the demo data
# Demo data was stored in Data/ directory of GPWAS package
> gpwas(ingeno='GPWAS-demo.geno', inpheno='GPWAS-demo.pheno', pc=3, inpc = 'population-structure-demo', gp='output-geno-phenotypes.txt', gv='output-geno-pvalue.txt', R=5)
```

# GPWAS genes selection

After running GPWAS model for collected phenotype and genotype data in a given population, you will obtain two files **gp** and **gv** (depends on your provided file names) with stored selected phenotypes per gene and p value per gene. GPWAS genes could be selected upon significant level per gene. However, if you are willing to set a threshold for GPWAS genes selection, we recommend you to shuffle each phenotype across all genotypes in your original phenotype data matrix for N times, then re-run GPWAS on each shuffled phenotype matrix and merge p-value per gene in your N **gv** outputs. Finally, you will get two list of p-value per gene: real data **gv** and N permutation data **gv**. When setting a p-value threshold **P**, you will know **R** genes in real data **gv** above **P** and **M** genes in permutation data above **P**. FDR could be determined using
```
FDR = (M/N)/R
```
A better estimation on FDR can be achieved through increasing permutation time. P-value threshold can be determined based on a permutation-based FDR approach using the R code deposited under ```utilis/```.
# Speed up the computation of GPWAS

Depend on the size of actual data matrix you have, we recommend you to split your genotype matrix into multiple subsets if you have too many phenotypes or/and too dense SNP per gene or/and too many individuals in the given population. Then submitting jobs in parallel to a computing cluster would shorten computing time efficiently.

# Common error reports solution

```
Error in anova.mlm(fit1) : residuals have rank 3 < 4
Calls: gpwas -> anova -> anova.mlm
Execution halted
```
Solution: This is likely multiple SNPs for a certain gene contain exact same information. You should make sure there are no duplicated SNPs inside one gene, even though their SNP name is different

# phenotype data shuffling suggestion
This is one suggested methods for shuffling the phenotype data.
```
# this is for one-round shuffling of your phenotype data.
before = data.frame(matrix(c(1:90), nrow = 9)) # your original phenotype data
after = data.frame(row.names=1:9)
for (x in c(1:ncol(before))){
  set.seed(100) # the number indicates one round of phentype shulling. If you want to get multiple shuffled phenotypes, you HAVE to change the number!
  after = cbind(after, sample(before[,x]))
}
colnames(after) = colnames(before)
after # your shuffled phenotype data
```
