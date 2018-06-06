

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(stringr)
FullAnnot = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name")]
FullAnnot = FullAnnot[str_length(rownames(FullAnnot))==10,]
FullAnnot = FullAnnot[!FullAnnot$UCSC_RefGene_Name=="",]

set.seed(123)
cpg.pval = runif(dim(FullAnnot)[1])
names(cpg.pval) = FullAnnot$Name
cpg.pval = round(cpg.pval, 4)
##  make sure pvalue is not 0
cpg.pval[cpg.pval==0] = 0.001
##  make some cpg to be significant
ind = sample(dim(FullAnnot)[1], 1000)
cpg.pval[ind] = 0.001
save(cpg.pval, file = 'data/cpgtoy.RData', compress = 'xz')





