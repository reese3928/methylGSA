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
#' @param group A string. "all", "body" or "promoter". Default is "all".
#' @param method A string. "ORA" or "GSEA". Default is "ORA"
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
                            group = "all", method = "ORA", GS.list=NULL, 
                            GS.idtype = "SYMBOL", GS.type = "GO", 
                            minsize = 100, maxsize = 500){
    if(!is.vector(cpg.pval) | !is.numeric(cpg.pval) | is.null(names(cpg.pval)))
        stop("Input CpG pvalues should be a named vector")
    if(sum(cpg.pval==0)>0)
        stop("Input CpG pvalues should not contain 0")
    if(!is.list(GS.list)&!is.null(GS.list))
        stop("Input gene sets should be a list")
    method = match.arg(method, c("ORA", "GSEA"))
    GS.idtype = match.arg(GS.idtype,
                                c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"))
    if(!is.null(GS.list) & GS.idtype!="SYMBOL")
        GS.list = suppressMessages(lapply(GS.list, function(x)
            select(org.Hs.eg.db, x, columns = "SYMBOL",
                        keytype = GS.idtype)$SYMBOL))
    GS.type = match.arg(GS.type, c("GO", "KEGG", "Reactome"))
    group = match.arg(group, c("all", "body", "promoter"))

    if(is.null(FullAnnot)){
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
    # list is a gene and it's corresponding cpg pvalue
    rho = vapply(geneID.list, rhoScores, FUN.VALUE = 1)
    ## for each gene, compute its rho score
    minbeta = vapply(geneID.list, function(x)
        min(betaScores(x)), FUN.VALUE = 1)
    ## for each gene, compute its beta value

    flag = 0
    if(is.null(GS.list)){
        GS.list = getGS(names(geneID.list), GS.type = GS.type)
        flag = 1
    }
    GS.list = lapply(GS.list, na.omit)

    GS.sizes = vapply(GS.list, length, FUN.VALUE = 0)
    GS.list.sub = GS.list[GS.sizes>=minsize & GS.sizes<=maxsize]
    ## filter gene sets by their sizes
    message(length(GS.list.sub), " gene sets are being tested...")
    size = vapply(GS.list.sub, length, FUN.VALUE = 0)
    ID = names(GS.list.sub)
    gs.pval = rep(NA, length(GS.list.sub))
    Count = rep(NA, length(GS.list.sub))

    if(method=="ORA"){
        rhoadj = p.adjust(rho, method = "BH")
        DEgenes = names(rhoadj)[rhoadj<0.05]

        N = length(unique(FullAnnot.sub$UCSC_RefGene_Name))
        m = length(DEgenes)
        if(m==0){
            gs.pval = 1
            Count = 0 
            warning("No gene is significant under FDR 0.05.")
        }
        else{
            for(i in seq_along(GS.list.sub)){
                q = sum(DEgenes %in% GS.list.sub[[i]])
                Count[i] = q
                gs.pval[i] = phyper(q = q, m = m, n = N-m,
                                        k = size[i], lower.tail=FALSE)
            }
        }
        gs.padj = p.adjust(gs.pval, method = "BH")
        
        if(flag==1){
            des = getDescription(GSids = ID, GS.type = GS.type)
            res = data.frame(ID = ID, Description = des, Count = Count, 
                             Size = size, pvalue = gs.pval, padj = gs.padj)
        }
        else
            res = data.frame(ID = ID, Count = Count, Size = size,
                             pvalue = gs.pval, padj = gs.padj)
        rownames(res) = ID
        res = res[order(res$pvalue), ]
        message("Done!")
        return(res)
    }

    if(method == "GSEA"){
        GS2gene = data.frame(
            ont = rep(names(GS.list), vapply(GS.list, length, FUN.VALUE = 0)),
                gene = unlist(GS.list))
        zscore = qnorm(minbeta/2, lower.tail = FALSE)
        zscore = zscore[order(zscore, decreasing = TRUE)]
        zscore[is.infinite(zscore)] = suppressWarnings(
            max(zscore[-which(is.infinite(zscore))]))
        GSEAres = GSEA(geneList = zscore, minGSSize = minsize,
                            maxGSSize = maxsize, pvalueCutoff = 1,
                            TERM2GENE = GS2gene)
        GSEAres = GSEAres[,c("ID","setSize","enrichmentScore",
                                "NES","pvalue","p.adjust","leading_edge")]
        GSEAres = data.frame(GSEAres)
        if(flag==1){
            des = getDescription(GSids = GSEAres$ID, GS.type = GS.type)
            GSEAres$Description = des 
            GSEAres = GSEAres[,c("ID","Description","setSize","enrichmentScore",
                                 "NES","pvalue","p.adjust","leading_edge")]
        }
        colnames(GSEAres)[which(colnames(GSEAres)=="setSize")] = "Size"
        colnames(GSEAres)[which(colnames(GSEAres)=="p.adjust")] = "padj"
        message("Done!")
        return(GSEAres)
    }
}
