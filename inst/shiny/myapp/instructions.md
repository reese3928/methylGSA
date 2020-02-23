## Step by step instructions for methylGSA shiny app    
#### Xu Ren and Pei Fen Kuan  
#### 2020-02-20

This app is a web interface for gene set testing using the outcome of differential methylation. 

### Step 1
Click on "Browse" to upload differential methylation result. 

<img src="https://github.com/reese3928/methylGO-RShiny-app/raw/master/step1.png" style="width:40.0%">


Users are expected to upload a txt file with 1st column CpG IDs, and 2nd column p-values correspond to the CpGs. For example:

<hr />

cg13869341 &nbsp;&nbsp; 0.307766 <br>
cg14008030 &nbsp;&nbsp; 0.257672 <br>
cg12045430 &nbsp;&nbsp; 0.552322 <br>
cg20826792 &nbsp;&nbsp; 0.056383 <br>
cg00381604 &nbsp;&nbsp; 0.468549 <br>
cg20253340 &nbsp;&nbsp; 0.483770 <br>
cg21870274 &nbsp;&nbsp; 0.812402 <br>
... <br>
... <br>
... <br>
cg21106100 &nbsp;&nbsp; 0.079276 <br>
cg08265308 &nbsp;&nbsp; 0.748265 <br>
cg14273923 &nbsp;&nbsp; 0.553923 <br>

<hr />

Files should be no more than 100MB. Please check [**cpg_file_example.txt**](https://github.com/reese3928/methylGO-RShiny-app/raw/master/cpy_file_example.txt) for a toy dataset on Illumina 450 K Beadchip. 




### Step 2
Choose array type and group.

<img src="https://github.com/reese3928/methylGO-RShiny-app/raw/master/step2_1.png" style="width:40.0%">
<img src="https://github.com/reese3928/methylGO-RShiny-app/raw/master/step2_2.png" style="width:40.0%">


Array type:  
`450K`: Illumina 450 K Beadchip   
`EPIC`: Illumina EPIC Beadchip    

Group:  
Group is the type of CpG to be considered in methylRRA or methylglm. By default, group is `all`, which means all CpGs are considered regardless of their gene group. If group is `body`, only CpGs on gene body will be considered. If group is `promoter1` or `promoter2`, only CpGs on promoters will be considered. Based on the annotation in IlluminaHumanMethylation450kanno.ilmn12.hg19 and IlluminaHumanMethylationEPICanno.ilm10b4.hg19, `body`, `promoter1` and `promoter2` are defined as:    
`body`: CpGs whose gene group correspond to “Body” or “1stExon”    
`promoter1`: CpGs whose gene group correspond to “TSS1500” or “TSS200”    
`promoter2`: CpGs whose gene group correspond to “TSS1500”, “TSS200”, “1stExon”, or “5’UTR”


### Step 3
Select gene sets to test.

<img src="https://github.com/reese3928/methylGO-RShiny-app/raw/master/step3.png" style="width:40.0%">


`Gene Ontology`: http://www.geneontology.org  
`KEGG` (Kyoto Encyclopedia of Genes and Genomes): https://www.genome.jp/kegg/  
`Reactome`: https://reactome.org


### Step 4
Select gene set sizes.

<img src="https://github.com/reese3928/methylGO-RShiny-app/raw/master/step4.png" style="width:40.0%">

Gene sets outsize this range will not be tested. 


### Step 5
Select gene set testing method.

<img src="https://github.com/reese3928/methylGO-RShiny-app/raw/master/step5.png" style="width:40.0%">

Check [this page](https://bioconductor.org/packages/devel/bioc/vignettes/methylGSA/inst/doc/methylGSA-vignette.html) for a description of methods. 


### Step 6
Hit "GO!".

<img src="https://github.com/reese3928/methylGO-RShiny-app/raw/master/step6.png" style="width:40.0%">

Gene set testing is running in the background. Depending on the number of gene sets tested, it may take several seconds to several minutes. Once it is done, the result is going to show up on the right hand panel. 

<br>

### Interpreting the results

The result of gene set testing is presented in two formats, namely table and box plot.

Hit "Download as csv" or "Download as txt" to save the results table to a csv or txt file:

<img src="https://github.com/reese3928/methylGO-RShiny-app/raw/master/restable.png" style="width:70.0%">

<br>

Various options are provided for users to customize the box plot:

<img src="https://github.com/reese3928/methylGO-RShiny-app/raw/master/resplot1.png" style="width:70.0%">

`Number of gene sets to display`: The number of genes sets to display on the plot.  
`x-axis`: Either the number of significant genes in each gene set or the total number of genes. The number of significant genes is not available for methylglm and methylRRA(GSEA) because they are FCS methods.   
`Color by`: Color the barplot by either raw p-value or adjusted p-value. 


Hit "Download as pdf" or "Download as png" to save the boxplot to a pdf or png file:

<img src="https://github.com/reese3928/methylGO-RShiny-app/raw/master/resplot2.png" style="width:70.0%">




