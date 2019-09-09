#' @title Implement logistic regression adjusting
#' for number of probes in enrichment analysis
#'
#' @description This function implements logistic regression adjusting
#' for number of probes in enrichment analysis.
#' @param cpg.pval A named vector containing p-values of differential
#' methylation test. Names should be CpG IDs.
#' @param array.type A string. Either "450K" or "EPIC". Default is "450K".
#' This argument will be ignored if FullAnnot is provided.
#' @param FullAnnot A data frame provided by prepareAnnot function.
#' Default is NULL.
#' @param group A string. "all", "body" or "promoter". Default is "all".
#' @param GS.list A list. Default is NULL. If there is no input list,
#' Gene Ontology is used. Entry names are gene sets names, and elements
#' correpond to genes that gene sets contain.
#' @param GS.idtype A string. "SYMBOL", "ENSEMBL", "ENTREZID" or "REFSEQ".
#' Default is "SYMBOL"
#' @param GS.type A string. "GO", "KEGG", or "Reactome". Default is "GO"
#' @param minsize An integer. If the number of genes in a gene set is
#' less than this integer, this gene set is not tested. Default is 100.
#' @param maxsize An integer. If the number of genes in a gene set is greater
#' than this integer, this gene set is not tested. Default is 500.
#' @param parallel either TRUE or FALSE indicating whether parallel should be
#' used. Default is FALSE
#' @param BPPARAM an argument provided to \code{\link{bplapply}}. See
#' \code{\link[BiocParallel]{register}} for details.
#' @details The implementation of this function is modified from goglm
#' function in GOglm package.
#' @export
#' @import stats
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi select
#' @importFrom BiocParallel bplapply bpparam
#' @return A data frame contains gene set tests results.
#' @references Mi G, Di Y, Emerson S, Cumbie JS and Chang JH (2012)
#' Length bias correction in Gene Ontology enrichment analysis using
#' logistic regression. PLOS ONE, 7(10): e46128
#' @references Phipson, B., Maksimovic, J., and Oshlack, A. (2015).
#' missMethyl: an R package for analysing methylation data from Illuminas
#' HumanMethylation450 platform. Bioinformatics, btv560.
#' @references Carlson M (2017). org.Hs.eg.db: Genome wide annotation for
#' Human. R package version 3.5.0.
#' @examples
#' data(CpG2Genetoy)
#' data(cpgtoy)
#' data(GSlisttoy)
#' GS.list = GS.list[1:10]
#' FullAnnot = prepareAnnot(CpG2Gene)
#' res = methylglm(cpg.pval = cpg.pval, FullAnnot = FullAnnot,
#' GS.list = GS.list, GS.idtype = "SYMBOL")
#' head(res)

methylglm <- function(cpg.pval, array.type = "450K", FullAnnot = NULL,
                            group = "all", GS.list=NULL, GS.idtype = "SYMBOL", 
                            GS.type = "GO", minsize = 100, maxsize = 500,
                            parallel = FALSE, BPPARAM = bpparam()){
    if(!is.vector(cpg.pval) | !is.numeric(cpg.pval) | is.null(names(cpg.pval)))
        stop("Input CpG pvalues should be a named vector")
    if(sum(cpg.pval==0)>0)
        stop("Input CpG pvalues should not contain 0")
    if(!is.list(GS.list)&!is.null(GS.list))
        stop("Input gene sets should be a list")
    
    GS.idtype = match.arg(
        GS.idtype,c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"))
    if(!is.null(GS.list) & GS.idtype!="SYMBOL")
        GS.list = suppressMessages(
            lapply(GS.list, function(x)
                select(org.Hs.eg.db, x, columns = "SYMBOL",
                            keytype = GS.idtype)$SYMBOL))
    GS.type = match.arg(GS.type, c("GO", "KEGG", "Reactome"))
    group = match.arg(group, c("all", "body", "promoter"))
    
    stopifnot(length(minsize)==1)
    if(!is.numeric(minsize) | minsize<0)
        stop("minsize should be a positive number")
    stopifnot(length(maxsize)==1)
    if(!is.numeric(maxsize) | maxsize<0)
        stop("maxsize should be a positive number")
    if(maxsize<minsize){
        stop("maxsize should be greater than minsize")
    }
    
    stopifnot(length(parallel)==1)
    if(!is.logical(parallel)){
        stop("parallel should be either TRUE or FALSE")
    }

    if(is.null(FullAnnot)){
        stopifnot(length(array.type)==1)
        if(array.type!="450K" & array.type!="EPIC")
            stop("Input array type should be either 450K or EPIC")
        if(array.type=="450K")
            FullAnnot = getAnnot("450K", group)
        else
            FullAnnot = getAnnot("EPIC", group)
    }

    cpg.intersect = intersect(names(cpg.pval), rownames(FullAnnot))
    cpg.pval = cpg.pval[cpg.intersect]
    FullAnnot.sub = FullAnnot[names(cpg.pval), ]
    ## match user input CpG to our FullAnnot database
    names(cpg.pval) = FullAnnot.sub$UCSC_RefGene_Name
    ## change CpG ids to gene symbols
    geneID.list = split(cpg.pval, names(cpg.pval))
    ## convert cpg.pval to a list, each element of this
    #list is a gene and it's corresponding cpg pvalue
    gene.pval = vapply(geneID.list, min, FUN.VALUE = 1)
    ## for each gene, take the minimum p-value
    probes = vapply(geneID.list, length, FUN.VALUE = 1)
    ## get number of probes for each gene
    
    flag = 0
    if(is.null(GS.list)){
        GS.list = getGS(names(geneID.list), GS.type = GS.type)
        flag = 1
    }
    
    GS.list = lapply(GS.list, na.omit)

    GS.sizes = vapply(GS.list, length, FUN.VALUE = 1)
    GS.list.sub = GS.list[GS.sizes>=minsize & GS.sizes<=maxsize]
    ## filter gene sets by their sizes
    
    if(!parallel){
        message(length(GS.list.sub), " gene sets are being tested...")
    }else{
        numworkers = BPPARAM$workers
        if(!numworkers[[1]]){
            numworkers = 1
        }
        message(length(GS.list.sub), " gene sets are being tested... using ",
            numworkers, " workers.")
    }
    
    glmFit = function(i){
        gs = GS.list.sub[[i]]
        y = as.numeric(names(gene.pval)%in%gs)
        df = data.frame(NegLogP = -log(gene.pval), probes = log(probes), y = y)
        glm.fit = glm(y ~ NegLogP + probes, family = "quasibinomial",
            data = df, control = list(maxit = 25))
        sumry = summary(glm.fit)
        sign(coefficients(glm.fit)[[2]])*sumry$coef[ ,"Pr(>|t|)"][[2]]
    }
    if(parallel){
        temp = bplapply(seq_along(GS.list.sub), glmFit, BPPARAM = BPPARAM)
        gs.pval = vapply(temp, function(x) x[1], 1)
    }else{
        gs.pval = vapply(seq_along(GS.list.sub), glmFit, 1)
    }
    gs.pval[gs.pval<=0] = 1
    ID = names(GS.list.sub)
    size = vapply(GS.list.sub, length, FUN.VALUE = 1)

    gs.padj = p.adjust(gs.pval, method = "BH")
    if(flag==1){
        des = getDescription(GSids = ID, GS.type = GS.type)
        res = data.frame(ID = ID, Description = des, Size = size,
            pvalue = gs.pval, padj = gs.padj)
    }
    else
        res = data.frame(ID = ID, Size = size, pvalue = gs.pval, padj = gs.padj)
    rownames(res) = ID
    res = res[order(res$pvalue), ]
    message("Done!")
    return(res)
}

