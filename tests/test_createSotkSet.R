library(SOTK)
library(testthat)

context("Checking the initionalization function.")

test_that("Invalid/exclusive arguments", {
        pvNmfObj <- readRDS("../inst/testdata/PV.RDS")
        dataL <- list(PV = pvNmfObj)
        rankL <- list(PV = c(2:6))

        expect_error(SOSet(NMFobjL = pvNmfObj, NMFrankL = rankL, dataCol =  c("PV" = "#7570B3"), corMet = "spearman"),
                "ERROR::Check the names in NMFrankL whether they match with the names in NMFobjL.")

        expect_error(SOSet(NMFobjL = dataL, NMFrankL = rankL, dataCol =  c("PV" = "#7570B3"), corMet = "concordance_index"),
                "ERROR::Provide either pearson, spearman, or kendall.")
})