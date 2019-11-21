methylGSA
===========

`methylGSA` is a Bioconductor package and Shiny app for DNA methylation data length bias adjustment in gene set testing. 

The Bioconductor package can be found [here](https://bioconductor.org/packages/release/bioc/html/methylGSA.html).    
The Bioconductor package vignette can be found [here](https://bioconductor.org/packages/release/bioc/vignettes/methylGSA/inst/doc/methylGSA-vignette.html).   
The `methylGSA` paper can be found [here](https://doi.org/10.1093/bioinformatics/bty892).

Shiny app installation
------------
The following packages are required to be installed before launching the app.    
Packages from CRAN:
```{r}    
install.packages("DT")    
install.packages("ggplot2")       
install.packages("shinycssloaders")     
```

Packages from Bioconductor:    
If analyzing 450K array, `IlluminaHumanMethylation450kanno.ilmn12.hg19` needs to be installed.
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
```
If analyzing EPIC array, `IlluminaHumanMethylationEPICanno.ilm10b4.hg19` needs to be installed.
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
```

Launch the app
------------
After installation, the app can be launched with the following code:
```{r}
library(methylGSA)
methylGSA::runExample()
```

Step-by-step instructions
------------
A step-by-step instructions on the workflow of the app can be found [here](http://www.ams.sunysb.edu/%7Epfkuan/methylGSA/ShinyAppInstruction.html).

Citation
------------
Ren, X., Kuan, P. (2018). methylGSA: a Bioconductor package and Shiny app for DNA methylation data length bias adjustment in gene set testing. Bioinformatics, bty892.    

@article{ren2018methylgsa,    
title={methylGSA: a Bioconductor package and Shiny app for DNA methylation data length bias adjustment in gene set testing},    
author={Ren, Xu and Kuan, Pei Fen},    
journal={Bioinformatics},    
year={2018}    
}





