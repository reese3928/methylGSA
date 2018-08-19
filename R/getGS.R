#' @title Get Gene Sets
#'
#' @description This function gets gene sets information.
#' @param geneids A vector contains all gene ids of interest. Gene ids should
#' be gene symbol.
#' @param GS.type A string. "GO", "KEGG", or "Reactome".
#' @export
#' @import org.Hs.eg.db
#' @import reactome.db
#' @importFrom AnnotationDbi select
#' @return A list contains all gene sets of
#' interest and their corresponding genes.
#' @references Carlson M (2017). org.Hs.eg.db:
#' Genome wide annotation for Human. R package version 3.5.0.
#' @references Ligtenberg W (2017). reactome.db:
#' A set of annotation maps for reactome. R package version 1.62.0.
#' @examples
#' geneids = c("FKBP5", "NDUFA1", "STAT5B")
#' GO.list = getGS(geneids, "KEGG")
#' head(GO.list)

getGS = function(geneids, GS.type){
    message("retrieving ", GS.type, " sets...")
    if(GS.type == "KEGG")
        GS.type = "PATH"
    if(GS.type == "GO")
        GS.type = "GOALL"
    if(GS.type == "Reactome"){
        ## first convert id to entrezid to use reactome.db
        gene.entrez = suppressMessages(
            select(org.Hs.eg.db, geneids,
                        columns = "ENTREZID",keytype = "SYMBOL")$ENTREZID)
        GOdf = suppressMessages(
            select(reactome.db, gene.entrez,
                        columns = "REACTOMEID", keytype = "ENTREZID"))
        genesymbol = suppressMessages(
            select(org.Hs.eg.db, GOdf$ENTREZID,
                        columns = "SYMBOL", keytype = "ENTREZID")$SYMBOL)
        GS.type = "REACTOMEID"

    }
    else{
        GOs = suppressMessages(
            na.omit(unique(
                select(org.Hs.eg.db, geneids,
                            GS.type,keytype = "SYMBOL")[,GS.type])))
        GOdf = suppressMessages(
            select(org.Hs.eg.db, GOs, "SYMBOL", keytype = GS.type))
        genesymbol = GOdf$SYMBOL
    }
    names(genesymbol) = GOdf[,GS.type]
    GO.list = split(genesymbol, names(genesymbol))
    return(GO.list)
}


