News
=====

# locuszoomr 0.3.2
###### 18/08/2024

* Fix for SNPs with chromosome coordinate format in `link_LD()` (only works with 
`LDproxy` method).
* Fix for non-human ensembl databases e.g. mouse in `locus()`.
* Record ensembl version, organism and genome in locus objects.
* Bugfix: give an error message if gene is not found in ensembl database in 
`locus()`.

# locuszoomr 0.3.1
###### 28/06/2024

* Add toggle for using webGL in `scatter_plotly()`.
* Add height control in plotly functions.

# locuszoomr 0.3.0
###### 16/04/2024

* Allow `index_snp` to be a vector to highlight more than 1 SNP per region 
(suggested by Luke Pilling).
* Altered default colour scheme.
* Multiple improvements to plotly version.
* Added option to use the much faster `LDproxy` in `link_LD()`. This is now the
default option.
* Added support for plotting loci with eQTL data to show multiple genes in 
different colours.
* Added ability to overlay up/down pointing triangles to show sign of beta 
coefficient for significant SNPs.
* Added highlighting of selected genes with individual colours in the gene 
tracks in `locus_plot()`, `locus_ggplot()`, `genetracks()` and 
`gg_genetracks()`.
* Enable use of downloadable recombination rate track files from UCSC in 
`link_recomb()`, which is much faster when plotting multiple loci.

# locuszoomr 0.2.1
###### 17/02/2024

* Added labels to `locus_ggplot()` and `gg_scatter()` (thanks to Tom Willis).
* Improved error handling in `link_recomb()`
* Ensure index SNP is plotted on top in `locus_plot()` and `locus_ggplot()`.
* In `scatter_plot()` arguments `chromCol` and `sigCol` are replaced by `scheme` 
which now allows setting of the index SNP colour.

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
