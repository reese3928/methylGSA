
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(stringr)

FullAnnot = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name")]
FullAnnot = FullAnnot[str_length(rownames(FullAnnot))==10,]
FullAnnot = FullAnnot[!FullAnnot$UCSC_RefGene_Name=="",]
temp = vapply(strsplit(FullAnnot$UCSC_RefGene_Name,split=";"),
              '[', 1, FUN.VALUE=character(1))
FullAnnot$UCSC_RefGene_Name = temp
gene2cpg = split(FullAnnot$Name, FullAnnot$UCSC_RefGene_Name)
probes = unlist(lapply(gene2cpg, length))
probes_order = probes[order(probes)]

set.seed(123)
GS.sizes = sample(100:300, 200)  ## sample gene sizes
GS.sizes = GS.sizes[cumsum(GS.sizes)<length(probes_order)]
GS.sizes[length(GS.sizes)] =
    length(probes_order) - sum(GS.sizes[1:(length(GS.sizes)-1)])
GS.sizescum = cumsum(GS.sizes)
GS.list = vector("list", length(GS.sizes))
GS.list[[1]] = names(probes_order)[1:GS.sizes[1]]
names(GS.list)[1] = "GS1"
for(i in 2:length(GS.sizes)){
    GS.list[[i]] = names(probes_order)[ (GS.sizescum[i-1]+1):GS.sizescum[i] ]
    names(GS.list)[i] = paste0("GS", i)
}
save(GS.list, file = 'data/GSlisttoy.RData', compress = 'xz')


