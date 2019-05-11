

library(DT)
library(ggplot2)

server <- function(input, output) {
    
    options(shiny.maxRequestSize=100*1024^2)
    
    output$restext <- renderUI({
        HTML(paste("<b>", input$test.method, "analysis result:", "</b>")) })
    
    tableInput <- reactive({
        inFile <- eventReactive(input$go, input$cpg.pval)
        if(is.null(inFile()))
            return(NULL)
        
        temp = read.table(inFile()$datapath)
        cpg.pval1 = temp[,2]
        names(cpg.pval1) = temp[,1]
        
        inputmethod = eventReactive(input$go, input$test.method)
        
        if(input$array.type=="450K"){
            suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
        }else{
            suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))
        }
        
        if(inputmethod() == "methylglm"){
            res <- eventReactive(input$go, {
                methylglm(cpg.pval = cpg.pval1, array.type = input$array.type,
                          group = input$group, GS.list=NULL, 
                          GS.idtype = "SYMBOL", GS.type = input$GS.list, 
                          minsize = input$minsize, maxsize = input$maxsize)
            })
        }
        
        if(inputmethod() == "RRA(ORA)"){
            res <- eventReactive(input$go,{
                methylRRA(cpg.pval = cpg.pval1, array.type = input$array.type, 
                          group = input$group, method = "ORA", 
                          GS.list=NULL, GS.idtype = "SYMBOL", GS.type = input$GS.list, 
                          minsize = input$minsize, maxsize = input$maxsize)
            }) 
        }
        
        if(inputmethod() == "RRA(GSEA)"){
            res <- eventReactive(input$go,{
                methylRRA(cpg.pval = cpg.pval1, array.type = input$array.type, 
                          group = input$group, method = "GSEA", 
                          GS.list=NULL, GS.idtype = "SYMBOL", GS.type = input$GS.list, 
                          minsize = input$minsize, maxsize = input$maxsize)
            })
        }
        
        if(inputmethod() == "gometh"){
            res <- eventReactive(input$go,{
                methylgometh(cpg.pval = cpg.pval1, sig.cut = 0.001, array.type = input$array.type, 
                             GS.list=NULL, GS.idtype = "SYMBOL", GS.type = input$GS.list, 
                             minsize = input$minsize, maxsize = input$maxsize)
            })
        }
        
        res
        
    })
    
    output$resTable <- DT::renderDataTable({
        df = tableInput()
        DT::datatable(data=df()) %>% formatSignif(c("pvalue", "padj"), digits = 8)
    })
    
    plotInput <- reactive({
        inputmethod = eventReactive(input$go, input$test.method)
        
        res = tableInput()
        
        if(inputmethod() %in% c("methylglm","RRA(GSEA)")){
            validate(
                need(input$xaxis!="Count", 'Number of significant genes is not available for methlglm')
            )
        }
        
        barplot(res(), xaxis = input$xaxis, num = input$ngs, colorby = input$colorby)
        
    })
    
    output$resPlot = renderPlot({
        print(plotInput())
    })
    
    output$download1 <- downloadHandler(filename = function(){ paste0(input$test.method,"_result.csv") }, 
                                        content = function(fname){ write.csv(tableInput()(), fname)})
    
    output$download2 <- downloadHandler(filename = function(){ paste0(input$test.method,"_result.txt") }, 
                                        content = function(fname){ write.table(tableInput()(), fname, quote = FALSE, row.names = FALSE)})
    
    output$downloadPlot1 <- downloadHandler(filename = function() { paste0(input$test.method, '_plot.pdf') },
                                            content = function(fname){ ggsave(fname, plotInput(), device = "pdf", width = 30, height = 30, units = "cm")})
    output$downloadPlot2 <- downloadHandler(filename = function() { paste0(input$test.method, '_plot.png') },
                                            content = function(fname){ ggsave(fname, plotInput(), device = "png", width = 30, height = 30, units = "cm")})
    
}


