setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/canine_database_figure")
getwd()
library(dplyr)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")

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

pData(gse_set)$data_processing[1]

data = getGEO("GSE27217")
datExpr = as.data.frame(exprs (data[[1]]))

library(tibble)
datExpr2 <- tibble::rownames_to_column(datExpr, "ID")

Meta(gse)$platform_id
GPL3738 <- getGEO('GPL3738', destdir=".")
Meta(GPL3738)$title
colnames(Table(GPL3738))
Table(GPL3738)[1:10,1:4]
Table(GPL3738)[1:10,c("ID","GB_ACC","Gene Title","Gene Symbol","ENTREZ_GENE_ID")]
annotation_data <- GPL3738@dataTable@table
annotated_expression_data <- full_join(annotation_data, 
                            datExpr2, 
                            by = c("ID" = "ID"))
write.csv(annotated_expression_data,
          file = "annotated_expression_data.csv")

annotated_expression_data <- 
  read.csv("annotated_expression_data.csv")
annotated_expression_data_with_CV <- 
  read.csv("annotated_expression_data_with_CV.csv")

at_least_2_tumor_list <- 
  read.csv("at_least_2_tumor_list.csv")

annotated_expression_data_with_CV_filtered = 
  annotated_expression_data_with_CV %>% 
  filter(Gene.Symbol %in% at_least_2_tumor_list[,1] == TRUE) %>% 
  filter(CV >= 0.5) %>% 
  select()



write.csv(annotated_expression_data_with_CV_filtered,
          file = "annotated_expression_data_with_CV_filtered.csv")
