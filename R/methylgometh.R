#' @title Adjusting number of probes in gene set testing
#' using gometh or gsameth in missMethyl
#'
#' @description This function calls gometh or gsameth function
#' in missMethyl package to adjust number of probes in gene set testing
#' @param cpg.pval A named vector containing p-values of
#' differential methylation test. Names should be CpG IDs.
#' @param sig.cut A numeric value indicating cut-off value for significant CpG.
#' Default is 0.001. This argument will be ignored if topDE is provided. 
#' @param topDE An integer. The top number of CpGs to be declared as 
#' significant.
#' @param array.type A string. Either "450K" or "EPIC". Default is "450K".
#' @param GS.list A list. Default is NULL. If there is no input list,
#' Gene Ontology is used. Entry names are gene sets names, and elements
#' correpond to genes that gene sets contain.
#' @param GS.idtype A string. "SYMBOL", "ENSEMBL", "ENTREZID" or "REFSEQ".
#' Default is "SYMBOL".
#' @param GS.type A string. "GO", "KEGG", or "Reactome"
#' @param minsize An integer. If the number of genes in a gene set
#' is less than this integer, this gene set is not tested. Default is 100.
#' @param maxsize An integer. If the number of genes in a gene set
#' is greater than this integer, this gene set is not tested. Default is 500.
#' @export
#' @import stats
#' @importFrom missMethyl gometh
#' @importFrom missMethyl gsameth
#' @import reactome.db
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#' @importFrom stringr str_length
#' @return A data frame contains gene set tests results.
#' @references Phipson, B., Maksimovic, J., and Oshlack, A. (2015).
#' missMethyl: an R package for analysing methylation data from Illuminas
#' HumanMethylation450 platform. Bioinformatics, btv560.
#' @references Ligtenberg W (2017). reactome.db: A set of annotation maps
#' for reactome. R package version 1.62.0.
#' @references Carlson M (2017). org.Hs.eg.db: Genome wide annotation
#' for Human. R package version 3.5.0.
#' @examples
#' \dontrun{
#' library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' data(cpgtoy)
#' res = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.001, GS.type = "KEGG",
#' minsize = 200, maxsize = 205)
#' head(res)
#' }


methylgometh <- function(cpg.pval, sig.cut = 0.001, topDE = NULL,
                                array.type = "450K",
                                GS.list=NULL, GS.idtype = "SYMBOL",
                                GS.type = "GO",
                                minsize = 100, maxsize = 500){
    ## check input
    if(!is.vector(cpg.pval)|!is.numeric(cpg.pval)|is.null(names(cpg.pval))){
        stop("Input CpG pvalues should be a named vector")
    }
    if(sum(cpg.pval==0)>0){
        stop("Input CpG pvalues should not contain 0")
    }
    if(!is.list(GS.list)&!is.null(GS.list)){
        stop("Input gene sets should be a list")
    }
        
    stopifnot(length(sig.cut)==1)
    if(!is.numeric(sig.cut) | sig.cut>=1 | sig.cut<=0){
        stop("sig.cut should be a number between 0 and 1")
    }
    if(!is.null(topDE)){
        if(!is.numeric(topDE)|floor(topDE)<=0){
            stop("topDE should be a positive integer")
        }
    }
    
    stopifnot(length(array.type)==1)
    if(array.type!="450K" & array.type!="EPIC"){
        stop("Input array type should be either 450K or EPIC")
    }
        
    GS.type = match.arg(GS.type, c("GO", "KEGG", "Reactome"))
    
    stopifnot(length(minsize)==1)
    if(!is.numeric(minsize) | minsize<0){
        stop("minsize should be a positive number")
    }
    stopifnot(length(maxsize)==1)
    if(!is.numeric(maxsize) | maxsize<0){
        stop("maxsize should be a positive number")
    }
    if(maxsize<minsize){
        stop("maxsize should be greater than minsize")
    }

    if(is.null(topDE)){
        ## if topDE is not provided, declare significant genes by sig.cut
        sig.cpg = names(cpg.pval)[cpg.pval < sig.cut]
        if(length(sig.cpg)==0){
            warning("No CpG is significant under cut-off ", sig.cut, 
                    ". Use a higher cut-off or \"topDE\" argument.")
        }else{
            message(length(sig.cpg), 
                    " CpGs are significant under cut-off ", sig.cut)
        }
    }else{
        stopifnot(length(topDE)==1)
        if(!is.numeric(topDE)|floor(topDE)<=0){
            stop("topDE should be a positive integer")
        }
        cpg.pval = cpg.pval[order(cpg.pval)]
        sig.cpg = names(cpg.pval)[seq_len(floor(topDE))]
    }
    
    if(is.null(GS.list) & (GS.type=="GO"|GS.type=="KEGG")){
        message(GS.type, " are being tested...")
        res = gometh(sig.cpg = sig.cpg, all.cpg = names(cpg.pval),
            collection = GS.type, array.type=array.type,
            plot.bias = FALSE, prior.prob = TRUE)
        res = res[res$N>=minsize & res$N<=maxsize,]
        res = res[order(res$P.DE),]
        res$ID = rownames(res)
        if(GS.type=="GO"){
            colnames(res) = c("Description", "Ont", "Size", "Count", 
                "pvalue", "padj", "ID")
        }else{
            colnames(res) = c("Description", "Size", "Count", 
                "pvalue", "padj", "ID")
        }
        ## re-calculate adjusted p-value
        res$padj = p.adjust(res$pvalue,method = "BH")
        message("Done!")
        return(res)
    }

    if(is.null(GS.list) & GS.type=="Reactome"){
        message(GS.type, " are being tested...")
        if(array.type=="450K"){
            tempAnnot = getAnnot("450K")
        }else{
            tempAnnot = getAnnot("EPIC")
        }
            
        temp = unique(tempAnnot$UCSC_RefGene_Name)

        gene.entrez = suppressMessages(
            select(org.Hs.eg.db, temp, columns = "ENTREZID",
                        keytype = "SYMBOL")$ENTREZID)
        reactome.df = suppressMessages(
            select(reactome.db, gene.entrez, columns = "REACTOMEID",
                keytype = "ENTREZID"))
        reactome.df = na.omit(reactome.df)
        reactom2entrez = reactome.df$ENTREZID
        names(reactom2entrez) = reactome.df$REACTOMEID
        reactome.list = split(reactom2entrez, names(reactom2entrez))
        reactome.list.sizes = vapply(reactome.list, length, FUN.VALUE = 1)
        reactome.list.sub =
            reactome.list[
                reactome.list.sizes>=minsize & reactome.list.sizes<=maxsize]

        res = gsameth(sig.cpg = sig.cpg, all.cpg = names(cpg.pval),
                            collection = reactome.list.sub,
                            array.type=array.type,
                            plot.bias = FALSE, prior.prob = TRUE)
        res = data.frame(res)
        res = res[res$N>=minsize & res$N<=maxsize,]
        ID = rownames(res)
        Description = getDescription(rownames(res), "Reactome")
        res = cbind(ID, Description, res)
        res = res[order(res[,"P.DE"]),]
        colnames(res) = c("ID", "Description", "Size", 
            "Count", "pvalue", "padj")
        ## re-calculate adjusted p-value
        res$padj = p.adjust(res$pvalue,method = "BH")
        message("Done!")
        return(res)
    }else{
        if(!is.null(GS.list) & GS.idtype!="ENTREZID"){
            message("Converting gene IDs...")
        GS.list.entrezid = suppressMessages( lapply(GS.list, function(x)
            select(org.Hs.eg.db, x, columns = "ENTREZID",
                        keytype = GS.idtype)$ENTREZID))
        }else{
            GS.list.entrezid = GS.list
        }
            
        GS.list.entrezid = lapply(GS.list.entrezid, na.omit)
        GS.sizes = vapply(GS.list.entrezid, length, FUN.VALUE = 1)
        GS.list.sub = GS.list.entrezid[GS.sizes>=minsize & GS.sizes<=maxsize]
        ## filter gene sets by their sizes
        message(length(GS.list.sub), " gene sets are being tested...")
        res = gsameth(sig.cpg = sig.cpg, all.cpg = names(cpg.pval),
                            collection = GS.list.sub, array.type=array.type,
                            plot.bias = FALSE, prior.prob = TRUE)
        res = data.frame(res)
        ID = rownames(res)
        res = cbind(ID, res)
        res = res[order(res[,"P.DE"]),]
        colnames(res) = c("ID", "Size", "Count", "pvalue", "padj")
        message("Done!")
        return(res)
    }
}


