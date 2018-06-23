

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(stringr)
FullAnnot = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name")]
FullAnnot = FullAnnot[str_length(rownames(FullAnnot))==10,]
FullAnnot = FullAnnot[!FullAnnot$UCSC_RefGene_Name=="",]
temp = vapply(strsplit(FullAnnot$UCSC_RefGene_Name,split=";"),
              '[', 1, FUN.VALUE=character(1))
FullAnnot$UCSC_RefGene_Name = temp
colnames(FullAnnot) = c("CpG", "Gene")
rownames(FullAnnot) = NULL
CpG2Gene = FullAnnot

save(CpG2Gene, file = 'data/CpG2Genetoy.RData', compress = 'xz')




