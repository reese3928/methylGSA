library(shinycssloaders)
ui <- navbarPage("methylGSA",
                 tabPanel("Main",
                          # Sidebar with a slider input for number of bins 
                          sidebarLayout(
                              sidebarPanel(
                                  HTML("Please upload CpG IDs and their p-values. See help page for instructions." 
                                  ),
                                  
                                  fileInput(inputId = "cpg.pval", 
                                            label = ""),
                                  
                                  fluidRow(
                                      column(6, 
                                             selectInput(inputId = "array.type", 
                                                         label = "Array type:",
                                                         choices = c("450K", "EPIC"), 
                                                         selected = "450K")
                                      ),
                                      column(6,
                                             selectInput(inputId = "group", 
                                                         label = "Group:",
                                                         choices = c("all", "body", "promoter1", "promoter2"), 
                                                         selected = "all")
                                      )
                                  ),
                                  
                                  selectInput(inputId = "GS.list", 
                                              label = "Gene sets tested:",
                                              choices = c("Gene Ontology" = "GO", "KEGG", "Reactome"), 
                                              selected = "GO"),
                                  
                                  sliderInput(inputId = "minsize", 
                                              label = "Minimum gene set size:",
                                              min = 0, max = 1000, value = 100, step = 10), 
                                  
                                  sliderInput(inputId = "maxsize", 
                                              label = "Maximum gene set size:",
                                              min = 0, max = 1000, value = 500, step = 10), 
                                  
                                  selectInput(inputId = "test.method", 
                                              label = "Test method:",
                                              choices = c("methylglm", "gometh", "RRA(ORA)", "RRA(GSEA)"), 
                                              selected = "methylglm"),
                                  
                                  actionButton("go","GO!")
                                  
                                  #HTML("Please note, it may take some time for the results to show up.")
                                  
                              ),
                              
                              
                              mainPanel(
                                  htmlOutput("restext"),
                                  tabsetPanel(type = "tabs",
                                              tabPanel("Table",
                                                       withSpinner(DT::dataTableOutput("resTable"), type = 8),
                                                       downloadButton('download1',"Download as csv"),
                                                       downloadButton('download2',"Download as txt"),
                                                       tags$head(tags$style("#restext{font-size: 25px;}"))),
                                              tabPanel("Plot",
                                                       fluidRow(
                                                           column(4, 
                                                                  numericInput(inputId = "ngs",
                                                                               label = "Number of gene sets to display:",
                                                                               value = 5,
                                                                               min = 1,
                                                                               max = NA)
                                                           ),
                                                           column(4,
                                                                  selectInput(inputId = "xaxis", 
                                                                              label = "x-axis:",
                                                                              choices = c("Number of significant genes" = "Count", "Total genes" = "Size"), 
                                                                              selected = "Size")
                                                           ),
                                                           column(4,
                                                                  selectInput(inputId = "colorby", 
                                                                              label = "Color by:",
                                                                              choices = c("Adjusted p-value" = "padj", "Raw p-value" = "pvalue"), 
                                                                              selected = "padj")
                                                                  
                                                           )
                                                         
                                                       ),
                                                       plotOutput("resPlot"),
                                                       downloadButton('downloadPlot1',"Download as pdf"),
                                                       downloadButton('downloadPlot2',"Download as png")
                                                       
                                              )
                                  )
                                  
                              )
                              
                              
                          )
                 ),
                 
                 tabPanel("Help",
                          fluidRow(
                              column(8,
                                     includeMarkdown("instructions.md")
                              )
                          )
                          
                 ),
                 
                 tabPanel("About",
                          fluidRow(
                              column(6,
                                     includeMarkdown("About.md")
                              )
                              
                          )
                          
                 )
                 
)


