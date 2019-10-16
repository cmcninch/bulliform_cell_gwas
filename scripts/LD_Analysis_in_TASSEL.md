#Compute pairwise LD (R^2) genomewide in TASSEL
---
####Colton McNinch (August 12, 2019)
---
##Step 1. Create separate `.hmp.txt` files for each chromosoms in TASSEL. Follow steps in `Create_Numeric_SNP_Matrix.md` file.
* Use the taxa, MAF (<0.05) and missing rate (<0.05) hapmap filtered file.
* In TASSEL do: Data -> Separate
* Save each of the 10 chromosome files in hapmap format in the following notation: widiv_chrom_XX, where XX is the chromosome name.

## Step 2. Used tassel from the command line to compute pairwise ld between each SNP and the 500 SNPs upstream and downstream of each SNP. Chromosome 01 is shown as an example.
```
./run_pipeline.pl -Xms512m -Xmx10g -fork1 -h chromosome_01.hmp.txt -sortPositions -ld -ldType SlidingWindow -ldWinSize 500 -export chrom_01_ld.txt
```

## Step 3. Only obtain the information needed for the resulting files. Chromosome 01 is shown as an example. Remaining of the analyses are done in R.
```
cut -f2,8,13,14 chrom_01_ld.txt > chrom_01_ld_final.txt
```
