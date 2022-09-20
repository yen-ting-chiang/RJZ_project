setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/metascape_gene_list_pheatmap")
getwd()

library(pheatmap)
library(RColorBrewer)
library(dplyr)
bone_related_gene_list <- 
  read.csv(file = "bone_related_gene_list.csv", 
           header = T)

bone_related_gene_list_tmp = 
  bone_related_gene_list[,c(2:5)]
row.names(bone_related_gene_list_tmp) = 
  bone_related_gene_list[,1]


pheatmap(bone_related_gene_list_tmp,
         color = colorRampPalette(brewer.pal(n =3, 
                                             name = "Blues"))(100),
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color = "#A0A0A0",
         na_col = "#E0E0E0",
         show_rownames = FALSE,
         show_colnames = FALSE,
         fontsize_row = 8,
         legend = FALSE)
