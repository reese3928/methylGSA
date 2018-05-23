#' @title Adjusting number of probes in gene set testing
#' using gometh or gsameth in missMethyl
#'
#' @description This function calls gometh or gsameth function
#' in missMethyl package to adjust number of probes in gene set testing
#' @param cpg.pval A named vector containing p-values of
#' differential methylation test. Names should be CpG IDs.
#' @param sig.cut A numeric value indicating cut-off value for significant CpG.
#' @param array.type A string. Either "450K" or "EPIC". Default is "450K".
#' @param GS.list A list. Default is NULL. If there is no input list,
#' Gene Ontology is used. Entry names are gene sets names, and elements
#' correpond to genes that gene sets contain.
#' @param GS.idtype A string. "SYMBOL", "ENSEMBL", "ENTREZID" or "REFSEQ".
#' Default is "SYMBOL"
#' @param GS.type A string. "GO", "KEGG", or "Reactome"
#' @param minsize An integer. If the number of genes in a gene set
#' is less than this integer, this gene set is not tested. Default is 100.
#' @param maxsize An integer. If the number of genes in a gene set
#' is greater than this integer, this gene set is not tested. Default is 500.
#' @export
#' @import stats
#' @import missMethyl
#' @import minfi
#' @import reactome.db
#' @import AnnotationDbi
#' @import org.Hs.eg.db
#' @import stringr
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import IlluminaHumanMethylationEPICanno.ilm10b3.hg19
#' @return A data frame contains gene set tests results.
#' @references Phipson, B., Maksimovic, J., and Oshlack, A. (2015).
#' missMethyl: an R package for analysing methylation data from Illuminas
#' HumanMethylation450 platform. Bioinformatics, btv560.
#' @references Ligtenberg W (2017). reactome.db: A set of annotation maps
#' for reactome. R package version 1.62.0.
#' @references Carlson M (2017). org.Hs.eg.db: Genome wide annotation
#' for Human. R package version 3.5.0.
#' @references Hansen KD (2016). IlluminaHumanMethylation450kanno.ilmn12.hg19:
#' Annotation for Illumina's 450k methylation arrays. R package version 0.6.0.
#' @references Hansen KD (2017). IlluminaHumanMethylationEPICanno.ilm10b3.hg19:
#' Annotation for Illumina's EPIC methylation arrays. R package version 0.6.0,
#'  https://bitbucket.com/kasperdanielhansen/Illumina_EPIC.
#' @examples
#' library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' data(cpgtoy)
#' res = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.001,
#' minsize = 200, maxsize = 210)
#' head(res, 15)

methylgometh <- function(cpg.pval, sig.cut, array.type = "450K",
                                GS.list=NULL, GS.idtype = "SYMBOL",
                                GS.type = "GO",
                                minsize = 100, maxsize = 500){
    if(!is.vector(cpg.pval) | !is.numeric(cpg.pval) | is.null(names(cpg.pval)))
        stop("Input CpG pvalues should be a named vector\n")
    if(sum(cpg.pval==0)>0)
        stop("Input CpG pvalues should not contain 0\n")
    if(!is.list(GS.list)&!is.null(GS.list))
        stop("Input gene sets should be a list\n")
    if(!is.numeric(sig.cut) | sig.cut>=1 | sig.cut<=0)
        stop("sig.cut should be a number between 0 and 1\n")
    GS.type = match.arg(GS.type, c("GO", "KEGG", "Reactome"))

    sig.cpg = names(cpg.pval)[cpg.pval < sig.cut]

    if(is.null(GS.list) & (GS.type=="GO"|GS.type=="KEGG")){
        cat(paste(GS.type, "are being tested...\n"))
        res = gometh(sig.cpg = sig.cpg, all.cpg = names(cpg.pval),
                            collection = GS.type, array.type=array.type,
                            plot.bias = FALSE, prior.prob = TRUE)
        res = res[res$N>=minsize & res$N<=maxsize,]
        res = res[order(res$P.DE),]
        cat("Done!\n")
        return(res)
    }

    if(is.null(GS.list) & GS.type=="Reactome"){
        cat(paste(GS.type, "are being tested...\n"))
        if(array.type=="450K")
            FullAnnot = getAnnotation(
                IlluminaHumanMethylation450kanno.ilmn12.hg19::
                    IlluminaHumanMethylation450kanno.ilmn12.hg19)
        else
            FullAnnot = getAnnotation(
                IlluminaHumanMethylationEPICanno.ilm10b3.hg19::
                    IlluminaHumanMethylationEPICanno.ilm10b3.hg19)

        FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name")]
        FullAnnot = FullAnnot[str_length(rownames(FullAnnot))==10,]
        FullAnnot = FullAnnot[!FullAnnot$UCSC_RefGene_Name=="",]
        temp = sapply(strsplit(FullAnnot$UCSC_RefGene_Name,split=";"), '[', 1)

        gene.entrez = suppressMessages(
            select(org.Hs.eg.db, temp, columns = "ENTREZID",
                        keytype = "SYMBOL")$ENTREZID)
        reactome.df = suppressMessages(
            select(reactome.db, gene.entrez, columns = "REACTOMEID",
                        keytype = "ENTREZID"))
        reactom2entrez = reactome.df$ENTREZID
        names(reactom2entrez) = reactome.df$REACTOMEID
        reactome.list = split(reactom2entrez, names(reactom2entrez))
        reactome.list.sizes = unlist(lapply(reactome.list, length))
        reactome.list.sub =
            reactome.list[
                reactome.list.sizes>=minsize & reactome.list.sizes<=maxsize]

        res = gsameth(sig.cpg = sig.cpg, all.cpg = names(cpg.pval),
                            collection = reactome.list.sub,
                            array.type=array.type,
                            plot.bias = FALSE, prior.prob = TRUE)
        res = res[res[,1]>=minsize & res[,1]<=maxsize,]
        res = res[order(res[,3]),]
        cat("Done!\n")
        return(res)
    }

    else{
        if(!is.null(GS.list) & GS.idtype!="ENTREZID"){
            cat("Converting gene IDs...\n")
        GS.list.entrezid = suppressMessages( lapply(GS.list, function(x)
            return( select(org.Hs.eg.db, x, columns = "ENTREZID",
                                keytype = GS.idtype)$ENTREZID )) )
        }

        else
            GS.list.entrezid = GS.list

        GS.list.entrezid = lapply(GS.list.entrezid, na.omit)
        GS.sizes = unlist(lapply(GS.list.entrezid, length))
        GS.list.sub = GS.list.entrezid[GS.sizes>=minsize & GS.sizes<=maxsize]
        ## filter gene sets by their sizes
        cat(length(GS.list.sub), "gene sets are being tested...\n")
        res = gsameth(sig.cpg = sig.cpg, all.cpg = names(cpg.pval),
                            collection = GS.list.sub, array.type=array.type,
                            plot.bias = FALSE, prior.prob = TRUE)
        res = res[order(res[,3]),]
        cat("Done!\n")
        return(res)
    }
}

