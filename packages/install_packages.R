##########################################
#  http://www.majorbio.com/
#  Copyright (C) in 2019 Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

#' @ sourceDir(path = "F:\\MAJORBIO\\R Course\\Intro\\packages") in windows OS
#' @ sourceDir(path = "/mnt/ilustre/centos7users/meng.luo/project/") in Linux OS
#'


sourceInstall <- function(path, win = TRUE,linux = FALSE, packages = TRUE, ...) {
  print(paste("R code running in ", R.version$platform, sep = " "))
  print(paste("Current", R.version$version.string, sep = " "))
  if (path & win & packages) {
    for (nm in list.files(path)) {
      install.packages(paste(path, nm, sep = "\\"))
    }
  } else if (win & packages) {
    install.packages("ggplot2")
    install.packages("RColorBrewer")
    install.packages("circlize")
    install.packages("Rcpp")
    install.packages("VennDiagram")
    install.packages("ggsci")
    if (!requireNamespace('BiocManager', quietly = TRUE))
      install.packages('BiocManager')
    BiocManager::install('ComplexHeatmap')
    BiocManager::install('ComplexHeatmap')
    BiocManager::install('EnhancedVolcano')
    BiocManager::install("DESeq2")
    #devtools::install_github("jokergoo/circlize")
  } else if (linux & packages) {
    install.packages("ggplot2")
    install.packages("RColorBrewer")
    install.packages("circlize")
    install.packages("Rcpp")
    install.packages("VennDiagram")
    install.packages("ggsci")
    if (!requireNamespace('BiocManager', quietly = TRUE))
      install.packages('BiocManager')
    BiocManager::install('ComplexHeatmap')
    BiocManager::install('ComplexHeatmap')
    BiocManager::install('EnhancedVolcano')
    BiocManager::install("DESeq2")
  } else {
    install.packages(list.files())
  }
}
