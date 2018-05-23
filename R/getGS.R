#' @title Get Gene Sets
#'
#' @description This function gets gene sets information.
#' @param geneids A vector contains all gene ids of interest.
#' @param GS.type A string. "GO", "KEGG", or "Reactome".
#' @import org.Hs.eg.db
#' @import reactome.db
#' @import AnnotationDbi
#' @return A list contains all gene sets of
#' interest and their corresponding genes.
#' @references Carlson M (2017). org.Hs.eg.db:
#' Genome wide annotation for Human. R package version 3.5.0.
#' @references Ligtenberg W (2017). reactome.db:
#' A set of annotation maps for reactome. R package version 1.62.0.

getGS = function(geneids, GS.type){
    cat(paste("retrieving", GS.type, "sets...\n"))
    if(GS.type == "GO"){
        GOs = suppressMessages(
            na.omit(unique(
                select(org.Hs.eg.db, geneids, "GO",keytype = "SYMBOL")$GO)))
        GOdf = suppressMessages(
            select(org.Hs.eg.db, GOs, "SYMBOL", keytype = "GO"))
        genesymbol = GOdf$SYMBOL
        names(genesymbol) = GOdf$GO
        GO.list = split(genesymbol, names(genesymbol))
        return(GO.list)
    }

    if(GS.type == "KEGG"){
        KEGGs = suppressMessages(
            na.omit(unique(
                select(org.Hs.eg.db, geneids, "PATH",keytype = "SYMBOL")$PATH)))
        KEGGdf = suppressMessages(
            select(org.Hs.eg.db, KEGGs, "SYMBOL", keytype = "PATH"))
        genesymbol = KEGGdf$SYMBOL
        names(genesymbol) = KEGGdf$PATH
        KEGG.list = split(genesymbol, names(genesymbol))
        return(KEGG.list)
    }

    if(GS.type == "Reactome"){
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
        return(reactome.list)
    }
}




