
#' @title Get gene set description
#'
#' @description This function gets description of gene sets.
#' @param GSids A vector contains gene set IDs.
#' @param GS.type A string. "GO", "KEGG", or "Reactome".
#' @export
#' @import GO.db
#' @importFrom clusterProfiler download_KEGG
#' @import reactome.db
#' @import AnnotationDbi
#' 
#' @return A vector contains gene sets description.
#' @references Carlson M (2018). GO.db: A set of annotation maps describing 
#' the entire Gene Ontology. R package version 3.6.0.
#' @references Yu G, Wang L, Han Y, He Q (2012). clusterProfiler: an R package 
#' for comparing biological themes among gene clusters. OMICS: A Journal 
#' of Integrative Biology, 16(5), 284-287.
#' @references Ligtenberg W (2017). reactome.db:
#' A set of annotation maps for reactome. R package version 1.62.0.
#' @examples
#' GSids = c("GO:0007389", "GO:0000978", "GO:0043062")
#' Description = getDescription(GSids, "GO")
#' head(Description)

getDescription <- function(GSids, GS.type){
    if(GS.type=="GO"){
        goterms = unlist(Term(GOTERM))
        temp = as.character(goterms[GSids])
        return(temp)
    }
    if(GS.type=="KEGG"){
        KEGGID2NAME = 
            download_KEGG("hsa", keggType="KEGG", 
                          keyType="kegg")$KEGGPATHID2NAME
        rownames(KEGGID2NAME) = KEGGID2NAME$from
        temp = as.character(KEGGID2NAME[paste0("hsa",GSids), "to"])
        return(temp)
    }
    if(GS.type=="Reactome"){
        Reactome2NAME = unlist(as.list(reactomePATHID2NAME))
        temp = as.character(Reactome2NAME[GSids])
        return(temp)
    }
    
}

