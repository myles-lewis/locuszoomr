News
=====

# locuszoomr 0.2.0
###### 21/12/2023

### New features
* Improved ggplot2 gene track plotting via `gg_genetracks()` to enable easy 
layering of several ggplot2 plots above a row of gene tracks (thanks to nickhir 
for the suggestion).
* For those that only want the gene tracks for their own plots, this is now 
easier by simply not specifying `data` (or setting it to `NULL`) when calling 
`locus()`.
* Added function `quick_peak()` for quickly finding peaks in GWAS datasets.
* Added function `link_recomb()` for retrieving recombination data from UCSC.
* Recombination rate is shown on a secondary y axis by `locus_plot()` and 
`locus_ggplot()`.
* Added `...` to `link_LD()` and `link_eqtl()` to allow passing of additional 
arguments such as `genome_build` to `LDlinkR` queries.

### Changes
* Argument `LDtoken` in `link_LD()` and `link_eqtl()` has been renamed `token` 
to be consistent with `LDlinkR`.

### Bugfixes
* Fixed bug when plotting LD with absent levels in `locus_ggplot()` and 
`locus_plotly()`.
* Fixed plots with no gene tracks (thanks to Tom Willis).
* Genes with missing gene symbols now display the ensembl gene ID.

# locuszoomr 0.1.3
###### 03/11/2023

* Added arrows to the gene tracks in `locus_plotly()`
* Fixed bug relating `yzero` argument in scatter plots
* Improved labelling
* Fixed CRAN ERROR relating to package EnsDb.Hsapiens.v75 in Suggests

# locuszoomr 0.1.2
###### 02/11/2023

* This is the initial build of locuszoomr
