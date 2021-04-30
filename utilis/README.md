### Determining P-value threshold from GPWAS based on the permutation-based FDR approach
Selecting a p-value based upon FDR to determine your p-value threshold
Using the command
```
Rscript FDR-GPWAS.R real-pvalue.txt merged-permutation-pvalue.txt 5
```
Note: ```real-pvalue.txt``` is the gene p-value output from GPWAS using real phenotype Data; ```merged-permutation-pvalue.txt``` is the gene p-value output from GPWAS using permuted phenotype data; ```5``` changes this to the number of permutation you run. 
