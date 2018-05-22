#' @title Enrichment analysis after adjusting multiple p-values of each gene by Robust Rank Aggregation
#'
#' @description This function implements enrichment after adjusting multiple p-values of each gene by Robust Rank Aggregation.
#' @param cpg.pval A named vector containing p-values of differential methylation test. Names should be CpG IDs.
#' @param array.type A string. Either "450K" or "EPIC". Default is "450K".
#' @param method A string. "ORA" or "GSEA". Default is "ORA"
#' @param GS.list A list. Default is NULL. If there is no input list, Gene Ontology is used. Entry names are gene sets names, and elements correpond to genes that gene sets contain.
#' @param GS.idtype A string. "SYMBOL", "ENSEMBL", "ENTREZID" or "REFSEQ". Default is "SYMBOL"
#' @param GS.type A string. "GO", "KEGG", or "Reactome". Default is "GO"
#' @param minsize An integer. If the number of genes in a gene set is less than this integer, this gene set is not tested. Default is 100.
#' @param maxsize An integer. If the number of genes in a gene set is greater than this integer, this gene set is not tested. Default is 500.
#' @export
#' @import stats
#' @import RobustRankAggreg
#' @import clusterProfiler
#' @import org.Hs.eg.db
#' @import AnnotationDbi
#' @return A data frame contains gene set tests results.
#' @references Kolde, Raivo, et al. Robust rank aggregation for gene list integration and meta-analysis. Bioinformatics 28.4 (2012): 573-580.
#' @references Phipson, B., Maksimovic, J., and Oshlack, A. (2015). missMethyl: an R package for analysing methylation data from Illuminas HumanMethylation450 platform. Bioinformatics, btv560.
#' @references Yu, Guangchuang, et al. clusterProfiler: an R package for comparing biological themes among gene clusters. Omics: a journal of integrative biology 16.5 (2012): 284-287.
#' @references Carlson M (2017). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.5.0.
#' @examples
#' library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' res1 = methylRRA(cpg.pval = cpg.pval, method = "ORA", minsize = 200, maxsize = 220)
#' head(res1, 15)
#' res2 = methylRRA(cpg.pval = cpg.pval, method = "GSEA", minsize = 200, maxsize = 220)
#' head(res2, 10)

methylRRA <- function(cpg.pval, array.type = "450K", method = "ORA", GS.list=NULL, GS.idtype = "SYMBOL", GS.type = "GO", minsize = 100, maxsize = 500){
  if(!is.vector(cpg.pval) | !is.numeric(cpg.pval) | is.null(names(cpg.pval)) )
    stop("Input CpG pvalues should be a named vector\n")
  if(sum(cpg.pval==0)>0)
    stop("Input CpG pvalues should not contain 0\n")
  if(!is.list(GS.list)&!is.null(GS.list))
    stop("Input gene sets should be a list\n")
  method = match.arg(method, c("ORA", "GSEA"))
  GS.idtype = match.arg(GS.idtype, c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"))
  if(!is.null(GS.list) & GS.idtype!="SYMBOL")
    GS.list = suppressMessages(lapply(GS.list, function(x)
      return( select(org.Hs.eg.db, x, columns = "SYMBOL", keytype = GS.idtype)$SYMBOL)) )
  GS.type = match.arg(GS.type, c("GO", "KEGG", "Reactome"))

  if(array.type!="450K" & array.type!="EPIC")
    stop("Input array type should be either 450K or EPIC\n")
  if(array.type=="450K")
    FullAnnot = getAnnot("450K")
  else
    FullAnnot = getAnnot("EPIC")

  cpg.intersect = intersect(names(cpg.pval), rownames(FullAnnot))
  cpg.pval = cpg.pval[cpg.intersect]
  FullAnnot.sub = FullAnnot[names(cpg.pval), ]  ## match user input CpG to our FullAnnot database
  names(cpg.pval) = FullAnnot.sub$UCSC_RefGene_Name   ## change CpG ids to gene symbols
  geneID.list = split(cpg.pval, names(cpg.pval))   ## convert cpg.pval to a list, each element of this list is a gene and it's corresponding cpg pvalue
  rho = unlist(lapply(geneID.list, rhoScores))   ## for each gene, compute its rho score
  minbeta = unlist( lapply(geneID.list, function(x){ return(min(betaScores(x))) } ) )  ## for each gene, compute its beta value

  if(is.null(GS.list))
    GS.list = getGS(names(geneID.list), GS.type = GS.type)
  GS.list = lapply(GS.list, na.omit)

  GS.sizes = unlist(lapply(GS.list, length))
  GS.list.sub = GS.list[GS.sizes>=minsize & GS.sizes<=maxsize]   ## filter gene sets by their sizes
  cat(length(GS.list.sub), "gene sets are being tested...\n")
  size = unlist(lapply(GS.list.sub, length))
  ID = names(GS.list.sub)
  gs.pval = rep(NA, length(GS.list.sub))

  if(method=="ORA"){
    rhoadj = p.adjust(rho, method = "BH")
    DEgenes = names(rhoadj)[rhoadj<0.05]

    N = length(unique(FullAnnot.sub$UCSC_RefGene_Name))
    m = length(DEgenes)
    if(m==0){
      gs.pval=1
      warning("No gene is significant under FDR 0.05.")
    }
    else{
      for(i in 1:length(GS.list.sub)){
        q = sum(DEgenes %in% GS.list.sub[[i]])
        gs.pval[i] = phyper(q = q, m = m, n = N-m, k = size[i], lower.tail=FALSE)
      }
    }
    gs.padj = p.adjust(gs.pval, method = "BH")
    res = data.frame(ID = ID, size = size, pvalue = gs.pval, padj = gs.padj)
    rownames(res) = ID
    res = res[order(res$pvalue), ]
    cat("Done!\n")
    return(res)
  }

  if(method == "GSEA"){
    GS2gene = data.frame(ont = rep(names(GS.list), sapply(GS.list, length)), gene = unlist(GS.list))
    zscore = qnorm(rho/2, lower.tail = FALSE)
    zscore = zscore[order(zscore, decreasing = TRUE)]
    zscore[is.infinite(zscore)] = suppressWarnings(max(zscore[-which(is.infinite(zscore))]))
    GSEAres = GSEA(geneList = zscore, minGSSize = minsize, maxGSSize = maxsize, pvalueCutoff = 1, TERM2GENE = GS2gene)
    GSEAres = GSEAres[,c(1,3,4,5,6,7,10)]
    cat("Done!\n")
    return(GSEAres)
  }
}
