context("Test user input")

data(cpgtoy)
data(GSlisttoy)

test_that("Invalid cpg pvalue input", {
    df = data.frame(cpgID = names(cpg.pval), pval = cpg.pval)

    expect_error(methylglm(cpg.pval = df))
    expect_error(methylgometh(cpg.pval = df, sig.cut = 0.001))
    expect_error(methylRRA(cpg.pval = df))
    expect_error(methylRRA(cpg.pval = runif(1000)))

    cpg.pval[1:10] = 0
    expect_error(methylgometh(cpg.pval = cpg.pval, sig.cut = 0.001))
    expect_error(methylRRA(cpg.pval = cpg.pval))
    expect_error(methylRRA(cpg.pval = cpg.pval))

})

test_that("Invalid GS.list input", {
    df = data.frame(
        GS = c(rep("GS1", length(GS.list[[1]])),
               rep("GS2", length(GS.list[[2]])),
               rep("GS3", length(GS.list[[3]]))
               ),
        gene = c(GS.list[[1]], GS.list[[2]], GS.list[[3]]))

    expect_error(methylglm(cpg.pval = cpg.pval, GS.list = df))
    expect_error(methylRRA(cpg.pval = cpg.pval, GS.list = df))
    expect_error(methylgometh(cpg.pval = cpg.pval,
                              sig.cut = 0.01, GS.list = df))
})


test_that("Invalid sig.cut", {
    expect_error(methylgometh(cpg.pval = cpg.pval, sig.cut = 1))
    expect_error(methylgometh(cpg.pval = cpg.pval, sig.cut = 0))
    expect_error(methylgometh(cpg.pval = cpg.pval, sig.cut = 10))
    expect_error(methylgometh(cpg.pval = cpg.pval, sig.cut = -1))
})



