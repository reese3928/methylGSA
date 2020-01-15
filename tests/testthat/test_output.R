context("Test function output")

data(CpG2Genetoy)
data(cpgtoy)
data(GSlisttoy)
GS.list = GS.list[1:10]
FullAnnot = prepareAnnot(CpG2Gene)

test_that("check for valid output", {
    res1 = methylglm(cpg.pval = cpg.pval, FullAnnot = FullAnnot,
                     GS.list = GS.list, GS.idtype = "SYMBOL",
                     minsize = 100, maxsize = 300)
    expect_is(res1, 'data.frame')
    expect_equal(dim(res1)[2], 4)
    expect_true(all(res1$pvalue>=0 & res1$pvalue<=1))
    expect_true(all(res1$padj>=0 & res1$padj<=1))
    expect_true(all(colnames(res1) %in% c("ID", "Size", "pvalue", "padj")))

    res2 = methylRRA(cpg.pval = cpg.pval, FullAnnot = FullAnnot,
                     method = "ORA", GS.list = GS.list)
    expect_is(res2, 'data.frame')
    expect_equal(dim(res2)[2], 5)
    expect_true(all(res2$pvalue>=0 & res2$pvalue<=1))
    expect_true(all(res2$padj>=0 & res2$padj<=1))
    expect_true(all(colnames(res2) %in% c("ID", "Count", 
                                          "Size", "pvalue", "padj")))

    res3 = methylRRA(cpg.pval = cpg.pval, FullAnnot = FullAnnot,
                     method = "GSEA", GS.list = GS.list)
    expect_is(res3, 'data.frame')
    expect_equal(dim(res3)[2], 8)
    expect_true("core_enrichment"%in%colnames(res3))
    expect_true(all(res3$pvalue>=0 & res3$pvalue<=1))
    expect_true(all(res3$p.adjust>=0 & res3$p.adjust<=1))
    expect_true(all(colnames(res3) %in%
                        c("ID","Size","enrichmentScore",
                          "NES","pvalue","padj","leading_edge",
                          "core_enrichment")))
})

