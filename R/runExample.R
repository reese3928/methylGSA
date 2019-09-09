#' @title methylGSA shiny app
#' 
#' @description This is an interface for Bioconductor package methylGSA.
#' @note In order to run the app, the following R/Bioconductor packages needs
#' to be installed properly: shinycssloaders, DT, ggplot2, 
#' IlluminaHumanMethylation450kanno.ilmn12.hg19 (if analyzing 450K array)
#' IlluminaHumanMethylationEPICanno.ilm10b4.hg19 (if analyzing EPIC array)
#' @param run Run the app or not. Default is TRUE
#' @importFrom shiny runApp
#' @return The shiny app will be opened in a web browser.
#' @export
#' @examples 
#' ## Please note: in this example, the argument run is set to be FALSE in 
#' ## order to pass R CMD check. However, when using the app, users are 
#' ## expected to launch the app by runExample()
#' runExample(FALSE)
runExample = function(run=TRUE) {
    if(run){
        appDir = system.file("shiny", "myapp", package = "methylGSA")
        if (appDir == "") {
            stop("Could not find example directory. 
                 Try re-installing `methylGSA`.", call. = FALSE)
        }
        runApp(appDir, display.mode = "normal")
        
    }
}

