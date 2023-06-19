setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OS_MeRIP_analysis")
getwd()

BiocManager::install("Rsamtools")
library(exomePeak)
library(Rsamtools)
library(rtracklayer)
GENE_ANNO_GTF <- import("gencode.v17.annotation.gtf")

GENE_ANNO_GTF <- 
IP_BAM <- scanBam("Galaxy53.bam")
INPUT_BAM <- scanBam("Galaxy52.bam")
result = exomepeak(GENE_ANNO_GTF=GENE_ANNO_GTF, 
                   IP_BAM=IP_BAM, 
                   INPUT_BAM=INPUT_BAM)