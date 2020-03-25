# test TransGeneID
test_that("Transform gene symbol to gene entrez ID", {
  expect_equal(TransGeneID("HLA-A"), c("HLA-A"="3105"))
})
test_that("Transform gene entrez ID to gene symbol", {
  expect_equal(TransGeneID("3105", fromType = "Entrez",
                           toType = "Symbol"), c("3105"="HLA-A"))
})
test_that("Transform gene symbols between organisms", {
  expect_equal(TransGeneID("Pvr", fromType = "Symbol",
                           toType = "Symbol", fromOrg = "mmu",
                           toOrg = "hsa"), c("Pvr"="PVR"))
})
