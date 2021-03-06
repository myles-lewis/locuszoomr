# locuszoomr

This is a pure R implementation of locuszoom for plotting genetic data at 
genomic loci accompanied by gene annotations. It uses the LDlink API to query 
linkage disequilibrium data from 1000 Genomes.

# Installation

Bioconductor packages `ensembldb` and an Ensembl database package are required 
before installation.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v75")
```

Install from Github
```
devtools::install_github("myles-lewis/locuszoomr")
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
```

`locuszoomr` automatically leverages the `LDlinkR` package to query 1000 genomes 
for linkage disequilibrium (LD) across SNPs. In order to make use of this API
function you will need a personal access token (see the `LDlinkR` vignette), 
available from the LDlink website https://ldlink.nci.nih.gov/?tab=apiaccess.

Requests to LDlink are cached using the `memoise` package, to reduce API 
requests. This is helpful when modifying plots for aesthetic reasons.

# Example locus plot

```
# Locus plot using SLE GWAS data from Bentham et al 2015
# FTP download full summary statistics from
# https://www.ebi.ac.uk/gwas/studies/GCST003156
library(data.table)
SLE_gwas <- fread('../bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv')

loc <- locus(SLE_gwas, gene = 'UBE2L3', flank = 1e5, LDtoken = "..")
plot(loc)
```
