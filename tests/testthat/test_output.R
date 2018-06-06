context("Test function output")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(cpgtoy)
data(GSlisttoy)
GS.list = GS.list[1:10]

test_that("check for valid output", {
    res1 = methylglm(cpg.pval = cpg.pval, GS.list = GS.list,
    GS.idtype = "SYMBOL", minsize = 100, maxsize = 300)
    expect_is(res1, 'data.frame')
    expect_equal(dim(res1)[2], 4)
    expect_true(all(res1$pvalue>=0 & res1$pvalue<=1))
    expect_true(all(res1$padj>=0 & res1$padj<=1))
    expect_true(all(colnames(res1) %in% c("ID", "size", "pvalue", "padj")))

    res2 = methylRRA(cpg.pval = cpg.pval, method = "ORA", GS.list = GS.list)
    expect_is(res2, 'data.frame')
    expect_equal(dim(res2)[2], 4)
    expect_true(all(res2$pvalue>=0 & res2$pvalue<=1))
    expect_true(all(res2$padj>=0 & res2$padj<=1))
    expect_true(all(colnames(res2) %in% c("ID", "size", "pvalue", "padj")))

    res3 = methylRRA(cpg.pval = cpg.pval, method = "GSEA", GS.list = GS.list)
    expect_is(res3, 'data.frame')
    expect_equal(dim(res3)[2], 7)
    expect_true(all(res3$pvalue>=0 & res3$pvalue<=1))
    expect_true(all(res3$p.adjust>=0 & res3$p.adjust<=1))
    expect_true(all(colnames(res3) %in%
                        c("ID","setSize","enrichmentScore",
                          "NES","pvalue","p.adjust","leading_edge")))

    res4 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.001,
                        GS.type = "KEGG", minsize = 200, maxsize = 205)
    expect_is(res4, 'data.frame')
    expect_true(all(c("N", "DE", "P.DE", "FDR")%in%colnames(res4)))
    expect_true(all(res4$P.DE>=0 & res4$P.DE<=1))
    expect_true(all(res4$FDR>=0 & res4$FDR<=1))

    #cpg.pval = 0.99
    #expect_warning(methylRRA(cpg.pval = cpg.pval,
                             #method = "ORA", GS.list = GS.list))
})

