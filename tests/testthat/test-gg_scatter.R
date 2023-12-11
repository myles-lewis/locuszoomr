test_that("gg_scatter can handle reduced no. of r2 intervals", {
  load(file = system.file('data', 'SLE_gwas_sub.rdata', package='locuszoomr', mustWork = T))

  library(EnsDb.Hsapiens.v75)

  filtered_r2_sle <- with(SLE_gwas_sub, SLE_gwas_sub[r2 > 0.0 & r2 < 0.2,])

  no_na_sle_loc <- locus(data = filtered_r2_sle, ens_db =  "EnsDb.Hsapiens.v75", gene = 'UBE2L3', flank = 1e5, labs = 'rsid', LD = 'r2')

  pl <- gg_scatter(no_na_sle_loc)

  expect_equal(class(pl), c('gg', 'ggplot'))
})
