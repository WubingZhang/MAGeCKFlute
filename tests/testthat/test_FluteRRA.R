# test_FluteRRA
test_that("FluteRRA pipeline - only gene summary file", {
  expect_message(
    FluteRRA(gene_summary = gzfile("../testdata/Pan2018_B16_Pmel1_rra.gene_summary.txt.gz"),
             proj = "testFluteRRA", organism = "mmu"))
})
test_that("FluteRRA pipeline - both gene summary and sgRNA summary", {
  expect_message(
    FluteRRA(gene_summary = gzfile("../testdata/Pan2018_B16_Pmel1_rra.gene_summary.txt.gz"),
             sgrna_summary = gzfile("../testdata/Pan2018_B16_Pmel1_rra.sgrna_summary.txt.gz"),
             proj = "testFluteRRA", organism = "mmu"))
})
