#TCGAbiolinks(TCGA RNAseq data download)-----------------------------------

#BiocManager::install("TCGAbiolinks")
#browseVignettes("TCGAbiolinks")
#BiocManager::install("SummarizedExperiment")
#install.packages("DT")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data") 
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

setwd("D:/OS")
getwd()
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
library(DT)

query <- GDCquery(project = "TARGET-OS",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  experimental.strategy = "RNA-Seq",
                  legacy = FALSE)
GDCdownload(query, method = "api", files.per.chunk = 10)
RNAseq_data <- GDCprepare(query, save = TRUE, save.filename = "TARGET-OS.rda")
RNAseq_data_matrix= as.data.frame(assay(RNAseq_data))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EDASeq")
library(GeneStructureTools)
RNAseq_data_matrix <- 
  tibble::rownames_to_column(RNAseq_data_matrix, "row_names")
RNAseq_data_matrix[,1] <- 
  removeVersion(RNAseq_data_matrix[,1])
RNAseq_data_matrix <- RNAseq_data_matrix %>%
  group_by(row_names) %>%
  filter(n() == 1)
RNAseq_data_matrix <- 
  tibble::column_to_rownames(RNAseq_data_matrix, "row_names")

library(EDASeq)
dataNorm <- TCGAanalyze_Normalization(RNAseq_data_matrix, 
                                      geneInfoHT, 
                          method = "geneLength")
write.csv(dataNorm, 
          file="TARGET-OS_dataNorm.csv")
write.csv(RNAseq_data_matrix, 
          file="TARGET-OS_data_matrix.csv")

clin.os <- GDCquery_clinic("TARGET-OS", "clinical")


clinical.data <- 
  read.csv(file = "TARGET_OS_ClinicalData_Discovery_20210520.csv")
YTHDF2_expression <- 
  read.csv(file = "YTHDF2_expression.csv")
YTHDF2_expression_normalized <- 
  read.csv(file = "YTHDF2_expression_normalized.csv")

joined_data <- left_join(YTHDF2_expression_normalized, 
          clinical.data, 
          by = c("Sample"= "TARGET.USI"))
joined_data <- joined_data %>% 
  arrange(ENSG00000198492)
write.csv(joined_data, file = "joined_data.csv")
joined_data_50_50 <- 
  read.csv(file = "joined_data_50_50_for_survival.csv")
install.packages("ranger")
install.packages("ggfortify")
library(survival)
library(ranger)
library(ggplot2)
library(ggfortify)



km_trt_fit <- survfit(Surv(Time.to.death.in.days, status) ~ expression, 
                      data=joined_data_50_50)
autoplot(km_trt_fit)
cox <- coxph(Surv(Time.to.death.in.days, status) ~ expression, 
             data = joined_data_50_50)
summary(cox)

cox2 <- coxph(Surv(Time.to.death.in.days, status) ~ Disease.at.diagnosis, 
             data = joined_data_50_50)
summary(cox2)

cox3 <- coxph(Surv(Time.to.death.in.days, status) ~ Histologic.response, 
              data = joined_data_50_50)
summary(cox3)
