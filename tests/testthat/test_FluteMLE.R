# test_FluteMLE
test_that("FluteMLE pipeline - single sample", {
  expect_message(
    FluteMLE(gzfile("../testdata/Pan2018_B16_Pmel1_mle.gene_summary.txt.gz"),
             treatname = "Pmel-1_IFNg", ctrlname = "Pmel-1_OT1_Ctrl_IFNg",
             proj = "testFluteMLE", organism = "mmu", norm_method = "cell_cycle"))
  })
test_that("FluteMLE pipeline - multiple samples", {
  expect_message(
    FluteMLE(gzfile("../testdata/Pan2018_B16_Pmel1_mle.gene_summary.txt.gz"),
             treatname = "Pmel-1_IFNg",
             ctrlname = c("Pmel-1_OT1_Ctrl_IFNg", "Pmel-1_Input_IFNg"),
             proj = "testFluteMLE", organism = "mmu"))
})
test_that("FluteMLE pipeline - Depmap", {
  expect_message(
    FluteMLE(gzfile("../testdata/Pan2018_B16_Pmel1_mle.gene_summary.txt.gz"),
             treatname = "Pmel-1_IFNg", ctrlname = "Depmap",
             incorporateDepmap = TRUE, proj = "testFluteMLE", organism = "mmu"))
})
