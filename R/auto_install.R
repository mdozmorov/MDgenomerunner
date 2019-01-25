#' Automatically install CRAN, Bioconductor, and/or GitHub packages
#' 
#' @param .cran_packages a vector of CRAN package names to be installed. Default - NULL, nothing to install.
#' @param .bioc_packages a vector of Bioconductor package names to be installed. Default - NULL, nothing to install.
#' @param .github_packages a vector of GitHub package names to be installed. Default - NULL, nothing to install.
#'
#' @return Nothing, just installs the packages
#' @export
#' @examples
#' 
#' \dontrun{
#' .cran_packages  <-  c("knitr", "phyloseqGraphTest", "phyloseq", "shiny",
#'                       "miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest",
#'                       "vegan", "plyr", "dplyr", "ggrepel", "nlme",
#'                       "reshape2", "devtools", "PMA","structSSI","ade4",
#'                       "igraph", "ggnetwork", "intergraph", "scales")
#' .bioc_packages <- c("phyloseq", "genefilter", "impute")
#' .github_packages <- c("jfukuyama/phyloseqGraphTest")
#' auto_install(.cran_packages, .bioc_packages, .github_packages)
#' }
#' 
#' @note Idea adapted from \link[https://github.com/spholmes/F1000_workflow/blob/master/src/analysis-setup.R]{F1000 microbiome workflow}
##

auto_install <- function(.cran_packages = NULL, .bioc_packages = NULL, .github_packages = NULL) {
  
  if (!is.null(.cran_packages)) {
    .inst <- .cran_packages %in% installed.packages()
    if (any(!.inst)) {
      install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
    }
  }
  
  if (!is.null(.bioc_packages)) {
    .inst <- .bioc_packages %in% installed.packages()
    if (any(!.inst)) {
      source("http://bioconductor.org/biocLite.R")
      biocLite(.bioc_packages[!.inst])
    }
  }
  
  if (!is.null(.github_packages)) {
    .inst <- .github_packages %in% installed.packages()
    if (any(!.inst)) {
      devtools::install_github(.github_packages[!.inst])
    }
  }
  
}
