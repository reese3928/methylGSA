## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(comment = "", message=FALSE, warning = FALSE)

## ------------------------------------------------------------------------
library(methylGO)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b3.hg19)

## ------------------------------------------------------------------------
#data(cpg.pval)
head(cpg.pval, 20)

## ------------------------------------------------------------------------
res1 = methylglm(cpg.pval = cpg.pval, minsize = 200, maxsize = 220)
head(res1, 15)

## ------------------------------------------------------------------------
res2 = methylRRA(cpg.pval = cpg.pval, method = "ORA", minsize = 200, maxsize = 220)
head(res2, 15)

## ------------------------------------------------------------------------
res3 = methylRRA(cpg.pval = cpg.pval, method = "GSEA", minsize = 200, maxsize = 220)
head(res3, 10)

## ------------------------------------------------------------------------
res4 = methylgometh(cpg.pval = cpg.pval, sig.cut = 0.001, minsize = 200, maxsize = 220)
head(res4, 15)

## ------------------------------------------------------------------------
#data(GS.list)
## to make the display compact, only a proportion of each gene set is shown
head(lapply(GS.list, function(x) return(x[1:30])), 3)   

## ------------------------------------------------------------------------
res5 = methylRRA(cpg.pval = cpg.pval, array.type = "450K", method = "ORA", GS.list = GS.list, GS.idtype = "SYMBOL", minsize = 100, maxsize = 300)
head(res5, 10)

## ------------------------------------------------------------------------
res6 = methylglm(cpg.pval = cpg.pval, array.type = "450K", GS.type = "Reactome", minsize = 100, maxsize = 200)
head(res6, 10)

