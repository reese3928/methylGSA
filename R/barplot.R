#' @title Barplot for methylGSA analysis result
#'
#' @description This function visualizes methylGSA analysis result by barplot.
#' 
#' @param res A data frame which contains methylGSA analysis result. 
#' @param xaxis A string which specify the x-axis in the barplot. 
#' Either "Size" (number of genes in gene set) or "Count" 
#' (number of significant genes in gene set). Default is "Size". "Count" option
#' is not available for methylglm and methylRRA(GSEA) result.
#' @param num An integer. Number of gene sets to display on the barplot. 
#' Default is 5.
#' @param colorby A string. Either "pvalue" or "padj". Default is "padj".
#' @param title A string. Barplot title. Default is NULL.
#' @import ggplot2
#' 
#' @details The implementation of the function is adapted from barplot function
#' in enrichplot package.
#' @export
#' @return ggplot object
#' @references Yu G (2018). enrichplot: Visualization of Functional Enrichment 
#' Result. R package version 1.0.2, https://github.com/GuangchuangYu/enrichplot.
#' @examples
#' res = data.frame(ID = c("04144", "04510", "04740", "04810", "05200"),
#'                  Description = c("Endocytosis", "Focal adhesion", 
#'                  "Olfactory transduction", 
#'                  "Regulation of actin cytoskeleton", "Pathways in cancer"),
#'                  Size = c(201, 200, 388, 213, 326),
#'                  pvalue = c(0.481, 0.696, 1, 1, 1),
#'                  padj = 1
#'                  )
#' barplot(res)

barplot <- function(res, xaxis = "Size", num = 5, 
                    colorby = "padj", title = ""){
    
    # if Description is not provided, use ID as text
    txt = ifelse("Description"%in%colnames(res), "Description", "ID")
    
    xaxis = match.arg(xaxis, c("Size", "Count"))
    colorby = match.arg(colorby, c("pvalue", "padj"))
    
    if(!is.numeric(num))
        stop("num should be an integer.")
    
    if(!"Count"%in%colnames(res)&xaxis=="Count"){
        warning("\"Count\" option not available. \"Size\" is used.")
        xaxis = "Size"
    }
    
    if("Description"%in%colnames(res)){
        res = res[!is.na(res$Description),]
        res$Description = factor(res$Description, 
            levels=rev(unique(res$Description)))
    }
    else{
        res = res[!is.na(res$ID),]
        res$ID = factor(res$ID, levels=rev(unique(res$ID)))
    }
    
    if(nrow(res)>num)
        res = res[seq_len(num),]
    
    ggplot(res, aes_string(x = txt, y = xaxis, fill = colorby)) +
        scale_fill_continuous(low="red", high="blue", name = colorby, 
            guide=guide_colorbar(reverse=TRUE))+
        geom_bar(stat = "identity") + coord_flip() +
        ggtitle(title) + xlab(NULL) + ylab(xaxis)
}


