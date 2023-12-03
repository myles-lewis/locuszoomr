test_that("Error is raised when locus has no SNPs", {
  load(file = system.file('data', 'SLE_gwas_sub.rdata', package='locuszoomr', mustWork = T))
  
  library(EnsDb.Hsapiens.v75)

  expect_error(locus(data = SLE_gwas_sub, ens_db = 'EnsDb.Hsapiens.v75', labs = 'rsid', seqname = 3, xrange = c(2e6, 100+2e6)), 'Locus contains no SNPs/datapoints')
})
