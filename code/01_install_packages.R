##########################################
#  http://www.majorbio.com/
#  Copyright (C) 2019 in Majorbio workshop
#  contact: meng.luo@majorbio.com
##########################################

# install R packages

install.packages("name")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("name")

devtools::install_github("name") 


# load packages

library(name) 

require(name) 


# examples 

install.packages('nVennR')
install.packages('VennDiagram')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

devtools::install_github("jokergoo/ComplexHeatmap")

library(ComplexHeatmap)
library(nVennR)
library(VennDiagram)

require(ComplexHeatmap)
require(nVennR)
require(VennDiagram)






