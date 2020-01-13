#' @title Enrichment analysis after adjusting multiple p-values of
#' each gene by Robust Rank Aggregation
#'
#' @description This function implements enrichment after adjusting
#' multiple p-values of each gene by Robust Rank Aggregation.
#' @param cpg.pval A named vector containing p-values of differential
#' methylation test. Names should be CpG IDs.
#' @param array.type A string. Either "450K" or "EPIC". Default is "450K".
#' This argument will be ignored if FullAnnot is provided.
#' @param FullAnnot A data frame provided by prepareAnnot function.
#' Default is NULL.
#' @param group A string. "all", "body", "promoter1" or "promoter2". 
#' Default is "all". If group = "body", only CpGs on gene body will be 
#' considered in methylRRA. If group = "promoter1" or group = "promoter2", 
#' only CpGs on promoters will be considered. Here is the definition of "body", 
#' "promoter1" and "promoter2" according to the annotation in 
#' IlluminaHumanMethylation450kanno.ilmn12.hg19 or 
#' IlluminaHumanMethylationEPICanno.ilm10b4.hg19. 
#' \itemize{
#'   \item body: CpGs whose gene group correspond to "Body" or "1stExon" 
#'   \item promoter1: CpGs whose gene group correspond to "TSS1500" or "TSS200"
#'   \item promoter2: CpGs whose gene group correspond to "TSS1500", "TSS200", 
#'   "1stExon", or "5'UTR". 
#' }
#' If group = "all", all CpGs are considered regardless of their gene group.
#' @param method A string. "ORA" or "GSEA". Default is "ORA"
#' @param sig.cut A numeric value indicating FDR cut-off for significant gene
#' in ORA. Default is 0.05. This argument will be ignored if topDE is provided
#' or method = "GSEA" is used.
#' @param topDE An integer. The top number of genes to be declared as 
#' significant after robust rank aggregation. This argument will be ignored
#' if method = "GSEA" is used.
#' @param GS.list A list. Default is NULL. If there is no input list,
#' Gene Ontology is used. Entry names are gene sets names, and elements
#' correpond to genes that gene sets contain.
#' @param GS.idtype A string. "SYMBOL", "ENSEMBL", "ENTREZID" or
#' "REFSEQ". Default is "SYMBOL".
#' @param GS.type A string. "GO", "KEGG", or "Reactome". Default is "GO"
#' @param minsize An integer. If the number of genes in a gene set is
#' less than this integer, this gene set is not tested. Default is 100.
#' @param maxsize An integer. If the number of genes in a gene set is
#' greater than this integer, this gene set is not tested. Default is 500.
#' @export
#' @import stats
#' @import RobustRankAggreg
#' @importFrom clusterProfiler GSEA
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi select
#' @return A data frame contains gene set tests results.
#' @references Kolde, Raivo, et al. Robust rank aggregation for gene
#' list integration and meta-analysis. Bioinformatics 28.4 (2012): 573-580.
#' @references Phipson, B., Maksimovic, J., and Oshlack, A. (2015).
#' missMethyl: an R package for analysing methylation data from Illuminas
#' HumanMethylation450 platform. Bioinformatics, btv560.
#' @references Yu, Guangchuang, et al. clusterProfiler: an R package for
#' comparing biological themes among gene clusters. Omics: a journal of
#' integrative biology 16.5 (2012): 284-287.
#' @references Carlson M (2017). org.Hs.eg.db: Genome wide annotation for
#' Human. R package version 3.5.0.
#' @examples
#' data(CpG2Genetoy)
#' data(cpgtoy)
#' data(GSlisttoy)
#' GS.list = GS.list[1:10]
#' FullAnnot = prepareAnnot(CpG2Gene)
#' res1 = methylRRA(cpg.pval = cpg.pval, FullAnnot = FullAnnot,
#' method = "ORA", GS.list = GS.list)
#' head(res1)

methylRRA <- function(cpg.pval, array.type = "450K", FullAnnot = NULL, 
                            group = "all", method = "ORA", sig.cut = 0.05, 
                            topDE = NULL, GS.list=NULL, 
                            GS.idtype = "SYMBOL", GS.type = "GO", 
                            minsize = 100, maxsize = 500){
    ##  check input
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
        stopifnot(length(topDE)==1)
        if(!is.numeric(topDE)|floor(topDE)<=0){
            stop("topDE should be a positive integer")
        }
    }
    
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
    
    ## match argument and make gene id to be gene symbol
    method = match.arg(method, c("ORA", "GSEA"))
    GS.idtype = match.arg(GS.idtype,
                                c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"))
    if(!is.null(GS.list) & GS.idtype!="SYMBOL"){
        GS.list = suppressMessages(lapply(GS.list, function(x)
            select(org.Hs.eg.db, x, columns = "SYMBOL",
                   keytype = GS.idtype)$SYMBOL))
    }
    GS.type = match.arg(GS.type, c("GO", "KEGG", "Reactome"))
    group = match.arg(group, c("all", "body", "promoter1", "promoter2"))

    ##  get annotation
    if(is.null(FullAnnot)){
        stopifnot(length(array.type)==1)
        if(array.type!="450K" & array.type!="EPIC"){
            stop("Input array type should be either 450K or EPIC")
        }
        if(array.type=="450K"){
            FullAnnot = getAnnot("450K", group)
        }else{
            FullAnnot = getAnnot("EPIC", group)
        }
            
    }

    cpg.intersect = intersect(names(cpg.pval), rownames(FullAnnot))
    cpg.pval = cpg.pval[cpg.intersect]
    ## match user input CpG to our FullAnnot database
    FullAnnot.sub = FullAnnot[names(cpg.pval), ]
    ## change CpG ids to gene symbols
    names(cpg.pval) = FullAnnot.sub$UCSC_RefGene_Name
    ## convert cpg.pval to a list, each element of this
    # list is a gene and it's corresponding cpg pvalue
    geneID.list = split(cpg.pval, names(cpg.pval))
    ## for each gene, compute its rho score
    rho = vapply(geneID.list, rhoScores, FUN.VALUE = 1)
    ## for each gene, compute its beta value
    #minbeta = vapply(geneID.list, function(x)
    #    min(betaScores(x)), FUN.VALUE = 1)
    
    flag = 0
    if(is.null(GS.list)){
        GS.list = getGS(names(geneID.list), GS.type = GS.type)
        flag = 1
    }
    GS.list = lapply(GS.list, na.omit)

    ## calculate gene set size
    GS.sizes = vapply(GS.list, length, FUN.VALUE = 0)
    ## filter gene sets by their sizes
    GS.list.sub = GS.list[GS.sizes>=minsize & GS.sizes<=maxsize]
    message(length(GS.list.sub), " gene sets are being tested...")
    size = vapply(GS.list.sub, length, FUN.VALUE = 0)
    ID = names(GS.list.sub)
    gs.pval = rep(NA, length(GS.list.sub))
    Count = rep(NA, length(GS.list.sub))

    if(method=="ORA"){
        if(is.null(topDE)){ 
            ## if topDE is not provided, declare significant genes by sig.cut
            rhoadj = p.adjust(rho, method = "BH")
            DEgenes = names(rhoadj)[rhoadj<sig.cut]
            if(length(DEgenes)==0){
                warning("No gene is significant under cut-off ", sig.cut, 
                        ". Use a higher cut-off or \"topDE\" argument.")
            }else{
                message(length(DEgenes), 
                        " genes are significant under cut-off ", sig.cut)
            } 
        }else{
            ## if topDE is provided, use topDE
            rho = rho[order(rho)]
            DEgenes = names(rho)[seq_len(floor(topDE))]
        }
        
        ##  total number of genes
        N = length(unique(FullAnnot.sub$UCSC_RefGene_Name))
        ##  number of DE genes
        m = length(DEgenes)
        if(m==0){
            gs.pval = 1
            Count = 0 
        }else{
            overlap.genes = lapply(GS.list.sub, 
                function(x) DEgenes[DEgenes%in%x])
            Q = vapply(overlap.genes, length, 0)
            Count = Q
            gs.pval = phyper(q = Q, m = m, n = N-m, k = size, lower.tail=FALSE)
            overlap.genes.vec = 
                vapply(overlap.genes, function(x) paste(x, collapse = ","), "")
        }
        gs.padj = p.adjust(gs.pval, method = "BH")
        
        if(flag==1){
            des = getDescription(GSids = ID, GS.type = GS.type)
            if(all(Count==0)){
                res = data.frame(ID = ID, Description = des, Count = Count, 
                    Size = size, pvalue = gs.pval, padj = gs.padj)
            }else{
                res = data.frame(ID = ID, Description = des, Count = Count, 
                    overlap = overlap.genes.vec, Size = size,
                    pvalue = gs.pval, padj = gs.padj)
            }
            
        }else{
            if(all(Count==0)){
                res = data.frame(ID = ID, Count = Count, Size = size,
                    pvalue = gs.pval, padj = gs.padj)
            }else{
                res = data.frame(ID = ID, Count = Count, 
                    overlap = overlap.genes.vec, Size = size, 
                    pvalue = gs.pval, padj = gs.padj)
            }
        }
        rownames(res) = ID
        res = res[order(res$pvalue), ]
        message("Done!")
        return(res)
    }

    if(method == "GSEA"){
        zscore = qnorm(rho/2, lower.tail = FALSE)
        zscore = zscore[order(zscore, decreasing = TRUE)]
        zscore[is.infinite(zscore)] = suppressWarnings(
            max(zscore[-which(is.infinite(zscore))]))
        
        ## drop the gene sets that has all zero z-scores
        allzero = vapply(GS.list, function(x){
            x = intersect(x, names(zscore))
            ifelse(all(na.omit(zscore[x])==0), TRUE, FALSE)
            }, FALSE) 
        GS.list = GS.list[!allzero]
        
        GS2gene = data.frame(
            ont = rep(names(GS.list), vapply(GS.list, length, FUN.VALUE = 0)),
                gene = unlist(GS.list))
        
        GSEAres = GSEA(geneList = zscore, minGSSize = minsize,
                            maxGSSize = maxsize, pvalueCutoff = 1,
                            TERM2GENE = GS2gene)
        GSEAres = GSEAres[,c("ID","setSize","enrichmentScore", "NES","pvalue",
            "p.adjust", "leading_edge","core_enrichment")]
        GSEAres = data.frame(GSEAres)
        if(flag==1){
            des = getDescription(GSids = GSEAres$ID, GS.type = GS.type)
            GSEAres$Description = des 
            GSEAres = GSEAres[,c("ID","Description","setSize","enrichmentScore",
                "NES","pvalue","p.adjust","leading_edge","core_enrichment")]
        }
        colnames(GSEAres)[which(colnames(GSEAres)=="setSize")] = "Size"
        colnames(GSEAres)[which(colnames(GSEAres)=="p.adjust")] = "padj"
        message("Done!")
        return(GSEAres)
    }
}
