#TCGAbiolinks(TCGA RNAseq data download)-----------------------------------

#BiocManager::install("TCGAbiolinks")
#browseVignettes("TCGAbiolinks")
#BiocManager::install("SummarizedExperiment")
#install.packages("DT")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data") 
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")


#Making TARGET-OS database transcriptome dataframe----------------------
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


#annotattion-----------------------------------------------
BiocManager::install("EnsDb.Hsapiens.v79")
library(EnsDb.Hsapiens.v79)
dataNorm <- 
  read.csv(file = "TARGET-OS_dataNorm.csv")
# 1. Convert from ensembl.gene to gene.symbol
names(dataNorm)[names(dataNorm) == 'X'] <- 'ensembl_GENEID'
ensembl.genes <- dataNorm$ensembl_GENEID
geneIDs1 <- as.data.frame((ensembldb::select(EnsDb.Hsapiens.v79, 
                                  keys= ensembl.genes, 
                                  keytype = "GENEID", 
                                  columns = c(colnames(),"GENEID"))))
dataNorm_annotated <- left_join(dataNorm, geneIDs1, 
                                by = c("ensembl_GENEID" = "GENEID"))
dataNorm_annotated <- dataNorm_annotated %>% 
  dplyr::select(ensembl_GENEID, SYMBOL, everything())
write.csv(dataNorm_annotated, 
          file = "dataNorm_annotated.csv")


#correlation analysis---------------------------------
dataNorm <- 
  read.csv(file = "TARGET-OS_dataNorm.csv")
library(dplyr)
library(ggplot2)
dataNorm_annotated_filtered <- 
  dataNorm_annotated %>% 
  dplyr::filter(is.na(SYMBOL)==FALSE) %>% 
  dplyr::select(-ensembl_GENEID)
dataNorm_t <- as.data.frame(t(dataNorm))

library(janitor)
dataNorm_t <- dataNorm_t %>% 
  row_to_names(row_number = 1)
dataNorm_t[,1:46264] <- lapply(dataNorm_t[,1:46264], as.numeric)

#YTHDF2, CDKN2B
ggplot(data = dataNorm_t,aes(x=ENSG00000147883,y=ENSG00000167508)) + 
  geom_point()+
  stat_smooth()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  )
cor_YTHDF2_CDKN2B <- 
  cor(dataNorm_t$ENSG00000198492,dataNorm_t$ENSG00000147883)
cor_CDKN2B_HMGCS1 <- 
  cor(dataNorm_t$ENSG00000147883,dataNorm_t$ENSG00000112972)
cor_CDKN2B_SQLE <- 
  cor(dataNorm_t$ENSG00000147883,dataNorm_t$ENSG00000104549)
cor_CDKN2B_MVK <- 
  cor(dataNorm_t$ENSG00000147883,dataNorm_t$ENSG00000110921)






dataMatrix <- 
  read.csv(file = "TARGET-OS_data_matrix.csv")
library(dplyr)
library(ggplot2)
dataNorm_annotated_filtered <- 
  dataNorm_annotated %>% 
  dplyr::filter(is.na(SYMBOL)==FALSE) %>% 
  dplyr::select(-ensembl_GENEID)
dataNorm_t <- as.data.frame(t(dataNorm))

library(janitor)
dataMatrix_t <- dataMatrix_t %>% 
  row_to_names(row_number = 1)
dataMatrix_t[,1:46264] <- lapply(dataMatrix_t[,1:46264], as.numeric)

#YTHDF2, CDKN2B
ggplot(data = dataMatrix_t,aes(x=ENSG00000147883,y=ENSG00000167508)) + 
  geom_point()+
  stat_smooth()+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
  )
cor_YTHDF2_CDKN2B <- 
  cor(dataMatrix_t$ENSG00000198492,dataMatrix_t$ENSG00000147883)
cor_CDKN2B_HMGCS1 <- 
  cor(dataMatrix_t$ENSG00000147883,dataMatrix_t$ENSG00000112972)
cor_CDKN2B_SQLE <- 
  cor(dataMatrix_t$ENSG00000147883,dataMatrix_t$ENSG00000104549)
cor_CDKN2B_MVK <- 
  cor(dataMatrix_t$ENSG00000147883,dataMatrix_t$ENSG00000110921)








#clinical data download and survival analysis----------------
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

#mutation data download------------------------------------------
query <- GDCquery(project = "TARGET-OS",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  experimental.strategy = "RNA-Seq",
                  legacy = FALSE)
GDCdownload(query, method = "api", files.per.chunk = 10)
RNAseq_data <- GDCprepare(query, save = TRUE, save.filename = "TARGET-OS.rda")
RNAseq_data_matrix= as.data.frame(assay(RNAseq_data))

query_snv <- GDCquery(
  project = "TARGET-OS",
  data.category = "Simple Nucleotide Variation",
  legacy = FALSE,
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query_snv)
snv_data <- GDCprepare(query_snv)

