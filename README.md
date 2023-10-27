# locuszoomr

This is a pure R implementation of locuszoom for plotting genetic data at 
genomic loci accompanied by gene annotations. It uses the LDlink API to query 
linkage disequilibrium data from 1000 Genomes.

# Installation

Bioconductor packages `ensembldb` and an Ensembl database installed either as a
package or obtained through Bioconductor packages `AnnotationHub` are required
before installation.

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v75")
```

Install from CRAN
```
install.packages("locuszoomr")
```

Install from Github
```
devtools::install_github("myles-lewis/locuszoomr")
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
# Using subset of data embedded in the package
library(locuszoomr)
data(SLE_gwas_sub)  ## limited subset of data from SLE GWAS

library(EnsDb.Hsapiens.v75)
loc <- locus(SLE_gwas_sub, gene = 'UBE2L3', flank = 1e5)
summary(loc)
locus_plot(loc)

# Or FTP download the full summary statistics from
# https://www.ebi.ac.uk/gwas/studies/GCST003156
library(data.table)
SLE_gwas <- fread('../bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv')

loc <- locus(SLE_gwas, gene = 'UBE2L3', flank = 1e5, LDtoken = "..")
plot(loc)
```
