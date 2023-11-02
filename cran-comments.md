## R CMD check results

0 errors | 0 warnings | 0 note

Version 0.1.3
* Fix CRAN NOTE: in vignette added conditional `if(require())` for 
Bioconductor package 'EnsDb.Hsapiens.v75', so that a vignette can still be made
even if EnsDb.Hsapiens.v75 package is missing.

* 'EnsDb.Hsapiens.v75' is still listed under Suggests as it is only needed for
examples and the vignette. But I suspect that will still result in a NOTE 
"Package suggested but not available for checking". Am wondering whether to move
package to Imports as then there will be no CRAN NOTE on checking?


Version 0.1.2
* Fixed CRAN ERROR: in examples added conditional `if(require())` for packages
in Suggests

* Fixed DESCRIPTION spaces and weblink for API
* Fixed return value for locus_ggplot.Rd

* This is a new release.
