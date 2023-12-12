test_that("Error is raised when locus has no SNPs", {
  load(file = system.file('data', 'SLE_gwas_sub.rdata', package='locuszoomr', mustWork = T))
  
  library(EnsDb.Hsapiens.v75)

  expect_error(locus(data = SLE_gwas_sub, ens_db = 'EnsDb.Hsapiens.v75', labs = 'rsid', seqname = 3, xrange = c(2e6, 100+2e6)), 'Locus contains no SNPs/datapoints')
})

test_that("Locus returns xrange of appropriate dimension when more than one gene matches specified name", {
  #16929816:16972118
  
  dummy_daf <- data.frame(chrom = 17, pos = 16929816, rsid = 'rs1', p = 1)
  
  library(EnsDb.Hsapiens.v86)
  
  expect_warning(locus(data = dummy_daf, ens_db = 'EnsDb.Hsapiens.v86', labs = 'rsid', gene = 'TNFRSF13B', flank = 1e5),  'Identified 2 genes matching name \'TNFRSF13B\', taking first\n')
 
  suppressWarnings({actual <- locus(data = dummy_daf, ens_db = 'EnsDb.Hsapiens.v86', labs = 'rsid', gene = 'TNFRSF13B', flank = 1e5)}) 
  
  expect_equal(actual$xrange, c(16829816, 17072118))
})

test_that("locus function exits successfully when no transcripts are found", {
  load(file = system.file('data', 'no_tx_data.rdata', package='locuszoomr', mustWork = T))
  
  library(EnsDb.Hsapiens.v86)
  
  loc <- locus(data = no_tx_data, chrom = 'chrom', pos = 'pos', p = 'p', labs = 'snp', seqname = 6, xrange = c(29833417, 29841260), ens_db = 'EnsDb.Hsapiens.v86')
 
  # Really just wanting the locus function to exit successfully, but need an assertion 
  expect_equal(class(loc), 'locus')  
})