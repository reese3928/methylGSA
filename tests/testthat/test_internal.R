context("Test internal functions")

library(org.Hs.eg.db)
library(reactome.db)
data(GSlisttoy)

test_that("check getGS", {
    geneids = GS.list[[1]]
    GOs = suppressMessages(
        na.omit(unique(
            select(org.Hs.eg.db, geneids, "GO",keytype = "SYMBOL")$GO)))
    GOdf = suppressMessages(
        select(org.Hs.eg.db, GOs, "SYMBOL", keytype = "GO"))
    genesymbol = GOdf$SYMBOL
    names(genesymbol) = GOdf$GO
    GO.list = split(genesymbol, names(genesymbol))

    expect_identical(getGS(geneids, "GO"), GO.list)

    KEGGs = suppressMessages(
        na.omit(unique(
            select(org.Hs.eg.db, geneids, "PATH",keytype = "SYMBOL")$PATH)))
    KEGGdf = suppressMessages(
        select(org.Hs.eg.db, KEGGs, "SYMBOL", keytype = "PATH"))
    genesymbol = KEGGdf$SYMBOL
    names(genesymbol) = KEGGdf$PATH
    KEGG.list = split(genesymbol, names(genesymbol))

    expect_identical(getGS(geneids, "KEGG"), KEGG.list)


    gene.entrez = suppressMessages(
        select(org.Hs.eg.db, geneids,
               columns = "ENTREZID",keytype = "SYMBOL")$ENTREZID)
    reactome.df = suppressMessages(
        select(reactome.db, gene.entrez,
               columns = "REACTOMEID", keytype = "ENTREZID"))
    reactom2symbol = suppressMessages(
        select(org.Hs.eg.db, reactome.df$ENTREZID,
               columns = "SYMBOL", keytype = "ENTREZID")$SYMBOL)
    names(reactom2symbol) = reactome.df$REACTOMEID
    reactome.list = split(reactom2symbol, names(reactom2symbol))

    expect_identical(getGS(geneids, "Reactome"), reactome.list)

})


