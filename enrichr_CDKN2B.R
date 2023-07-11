library(dplyr)
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/CDKN2B_enrichr_plot")
getwd()

library("enrichR")
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

gene_list_file <- read.csv(file = "DEG.csv")

db <- c("GO_Biological_Process_2021")

gene_list_file <- gene_list_file$Gene.name
if (websiteLive) {
  enriched  <- enrichr(gene_list_file, db)
}

if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "Combined.Score")
}

