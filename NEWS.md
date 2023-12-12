News
=====

# locuszoomr 0.1.5
###### 12/12/2023

* Improved ggplot2 gene track plotting via `gg_genetracks()` to enable easy 
layering of several ggplot2 plots above a row of gene tracks
* For those that only want the gene tracks for their own plots, this is now 
easier by simply not specifying `data` (or setting it to `NULL`) when calling 
`locus()`.

# locuszoomr 0.1.3
###### 03/11/2023

* Added arrows to the gene tracks in `locus_plotly()`
* Fixed bug relating `yzero` argument in scatter plots
* Improved labelling
* Fixed CRAN ERROR relating to package EnsDb.Hsapiens.v75 in Suggests

# locuszoomr 0.1.2
###### 02/11/2023

* This is the initial build of locuszoomr
