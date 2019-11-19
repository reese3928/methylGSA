#' @title Prepare user-supplied mapping between CpGs and genes.
#'
#' @description This function prepares CpG to gene mapping which will be
#' used by methylRRA and methylglm.
#' @param CpG2Gene A matrix, or a data frame or a list contains CpG to gene
#' mapping. For a matrix or data frame, 1st column should be CpG ID and 2nd
#' column should be gene name. For a list, entry names should be gene names,
#' and elements correpond to CpG IDs.
#' @param geneidtype A string. "SYMBOL", "ENSEMBL", "ENTREZID" or "REFSEQ".
#' Default is "SYMBOL".
#' @export
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi select
#' @return A data frame contains ready to use CpG to gene mapping.
#' @references Carlson M (2017). org.Hs.eg.db:
#' Genome wide annotation for Human. R package version 3.5.0.
#' @examples
#' data(CpG2Genetoy)
#' FullAnnot = prepareAnnot(CpG2Gene)
#' head(FullAnnot)

prepareAnnot <- function(CpG2Gene, geneidtype = "SYMBOL"){
    geneidtype = match.arg(
        geneidtype,c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"))

    if(is.matrix(CpG2Gene)|is.data.frame(CpG2Gene)){
        if(!is.character(CpG2Gene[,1])){
            stop("CpG ID should be characters")
        }
        if(ncol(CpG2Gene)!=2){
            stop("CpG2Gene should contain two columns")
        }
        FullAnnot = data.frame(CpG2Gene)
    }else if(is.list(CpG2Gene)){
        FullAnnot = data.frame(
            CpG = unlist(CpG2Gene),
            gene = rep(names(CpG2Gene),vapply(CpG2Gene, length, FUN.VALUE = 0)))
    }else{
        stop("CpG2Gene should be a matrix or a data frame or a list.")
    }
    
    colnames(FullAnnot) = c("Name", "UCSC_RefGene_Name")

    if(geneidtype!="SYMBOL"){
        temp = suppressMessages(
            select(org.Hs.eg.db, FullAnnot$UCSC_RefGene_Name,
                        columns = "SYMBOL",keytype = geneidtype))
        FullAnnot$UCSC_RefGene_Name = temp$SYMBOL
    }
    rownames(FullAnnot) = FullAnnot$Name
    return(FullAnnot)
}


