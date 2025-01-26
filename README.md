# locuszoomr

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/locuszoomr)](https://CRAN.R-project.org/package=locuszoomr)
[![Downloads](https://cranlogs.r-pkg.org/badges/locuszoomr)](https://CRAN.R-project.org/package=locuszoomr)
<!-- badges: end -->

This is a pure R implementation of locuszoom for plotting genetic data at
genomic loci accompanied by gene annotations. Plots can be produced in base
graphics, ggplot2 or plotly. Plots can be stacked or laid out with multiple
plots per page, or the gene track can be plotted separately and added to your
own plots.

The LDlink API can be queried to obtain linkage disequilibrium data from 1000
Genomes. Recombination rate can also be shown by querying UCSC genome browser.

See the detailed vignette for code examples.

# Installation

Bioconductor package `ensembldb` and an Ensembl database installed either as a
package or accessed through Bioconductor package `AnnotationHub` are required
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

`locuszoomr` can leverage the `LDlinkR` package to query the 1000 Genomes
Project for linkage disequilibrium (LD) across SNPs. In order to make use of
this API function you will need a personal access token, available from the 
[LDlink website](https://ldlink.nih.gov/?tab=apiaccess).

We recommend that users who want to add recombination rate lines to multiple
plots download the recombination rate track from UCSC and use it as described in
the 'Add recombination rate' section in the vignette.

# Example locus plot

```
# Locus plot using SLE GWAS data from Bentham et al 2015
# Using subset of data embedded in the package
library(locuszoomr)
data(SLE_gwas_sub)

library(EnsDb.Hsapiens.v75)
loc <- locus(gene = 'UBE2L3', SLE_gwas_sub, flank = 1e5,
             ens_db = "EnsDb.Hsapiens.v75")
summary(loc)
locus_plot(loc)

# Or FTP download the full summary statistics from
# https://www.ebi.ac.uk/gwas/studies/GCST003156
library(data.table)
SLE_gwas <- fread('../bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv')

loc <- locus(gene = 'UBE2L3', SLE_gwas, flank = 1e5,
             ens_db = "EnsDb.Hsapiens.v75")
locus_plot(loc)
```

# Example layered plot shown in the paper

```
library(locuszoomr)
library(EnsDb.Hsapiens.v75)

data(SLE_gwas_sub)
loc <- locus(data = SLE_gwas_sub, gene = 'IRF5', flank = c(1e5, 2e5),
             ens_db = "EnsDb.Hsapiens.v75")

# add recombination rate
loc <- link_recomb(loc, genome = "hg19")

# add LD and eQTL data
# users must obtain an API token from https://ldlink.nih.gov/?tab=apiaccess
loc <- link_LD(loc, token = "your_API_token")
loc <- link_eqtl(loc, token = "your_API_token")

# set up layered plot with 2 plots & a gene track
pdf("locuszoomr_demo.pdf", width = 4.5, height = 7)
oldpar <- set_layers(2)
scatter_plot(loc, xticks = FALSE, labels = c("index", "rs113708239"),
             label_x = c(-4, 4))
eqtl_plot(loc, xlab = "")
genetracks(loc, highlight = "IRF5")
par(oldpar)
dev.off()
```

# Citation

If you use this package please cite as:

Lewis MJ, Wang S. (2025) locuszoomr: an R package for visualising
publication-ready regional gene locus plots. *Bioinformatics Advances* 2025;
vbaf006, [doi:10.1093/bioadv/vbaf006](https://doi.org/10.1093/bioadv/vbaf006)
