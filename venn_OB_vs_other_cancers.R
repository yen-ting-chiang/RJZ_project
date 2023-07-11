setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/venn_OB_vs_other_cancers")
getwd()
library(VennDiagram)
DF <- read.csv(file = "mutation_list_for_R_vennDiagram.csv")
vlist <- list(DF$OB_derived_tumor, 
              DF$Cell_2019_driver_genes, 
              DF$NC_breast_cancer_Trp53_R245W,
              DF$PNAS_breast_cancer_Trp53_R245W)
names(vlist) <- c("OB", 
                  "driver_genes", 
                  "list3", 
                  "PNAS_breast_cancer_Trp53_R245W")
venn.diagram(vlist[1:2], 
             filename="Venn_OB_vs_Driver_gene.png",
             imagetype="png",
             fill=c("green", "blue"),
             col=c("green", "blue"),
             cat.col=c("green", "blue"),
             resolution = 300)

