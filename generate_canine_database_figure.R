setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/canine_database_figure")
getwd()

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

library(Biobase)

# BiocManager::install(c("GEOquery"))
library(GEOquery)
# install.packages("meta")
library(meta)

gse <- getGEO("GSE27217", destdir=".", GSEMatrix = FALSE)
Meta(gse)
gse_list <- getGEO("GSE27217", destdir=".", GSEMatrix = TRUE)
gse_set <- getGEO(filename='GSE27217_series_matrix.txt.gz')
Biobase::pData(gse_set)

data = getGEO("GSE27217")
datExpr = exprs (data[[1]])



