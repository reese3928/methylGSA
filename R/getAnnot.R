#' @title Get CpG annotation
#'
#' @description This function gets CpG IDs and their corresponding gene symbols.
#' @param array.type A string. Either "450K" or "EPIC". Default is "450K".
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#' @importFrom stringr str_length
#' @importFrom minfi getAnnotation
#' @details The implementation of the function is modified
#' from .flattenAnn function in missMethyl package.
#' @return A data frame contains CpG IDs and gene symbols.
#' @references Hansen KD (2016). IlluminaHumanMethylation450kanno.ilmn12.hg19:
#' Annotation for Illumina's 450k methylation arrays. R package version 0.6.0.
#' @references Hansen KD (2016). IlluminaHumanMethylationEPICanno.ilm10b2.hg19:
#' Annotation for Illumina's EPIC methylation arrays. R package version 0.6.0,
#' https://bitbucket.com/kasperdanielhansen/Illumina_EPIC.
#' @references Phipson B, Maksimovic J and Oshlack A (2015).
#' “missMethyl: an R package for analysing methylation data from
#' Illuminas HumanMethylation450 platform.” Bioinformatics, pp. btv560.

getAnnot = function(array.type){
    if(array.type=="450K")
        FullAnnot = getAnnotation(
            IlluminaHumanMethylation450kanno.ilmn12.hg19
            ::IlluminaHumanMethylation450kanno.ilmn12.hg19)
    else
        FullAnnot = getAnnotation(
            IlluminaHumanMethylationEPICanno.ilm10b2.hg19
            ::IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

    FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name")]
    FullAnnot = FullAnnot[str_length(rownames(FullAnnot))==10,]
    FullAnnot = FullAnnot[!FullAnnot$UCSC_RefGene_Name=="",]
    temp = sapply(strsplit(FullAnnot$UCSC_RefGene_Name,split=";"), '[', 1)
    ## get the first gene in each USCS_RefGene_Name
    FullAnnot$UCSC_RefGene_Name = temp

    return(FullAnnot)
}
