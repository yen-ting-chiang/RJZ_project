setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/Y2_KD_genes_heatmap")
getwd()
library(data.table)
library(dplyr)
counts_file <- fread(file = "Normalized_counts_file.csv")
WNT_genes <- fread(file = "WNT.csv")
counts_file_WNT <- 
  counts_file %>% 
  filter(gene_name %in% WNT_genes$SYMBOL == TRUE)
counts_file_WNT[counts_file_WNT==0] <- NA
counts_file_WNT
counts_file_WNT<-counts_file_WNT[complete.cases(counts_file_WNT),]
write.csv(counts_file_WNT,"counts_file_WNT_remove_zero.csv")


Normalized_WNT <- read.csv(file = "normalized_WNT.csv",header = TRUE)
library(pheatmap)
library(RColorBrewer)
rownames(Normalized_WNT) <- Normalized_WNT[, 1]
Normalized_WNT <- Normalized_WNT[, -1]
plot_heatmap <- pheatmap(Normalized_WNT,
                         color = colorRampPalette(rev(brewer.pal(n =9, 
                                                                 name = "RdBu")))(100),
                         cluster_rows=FALSE, 
                         cluster_cols=FALSE,
                         border_color = "#A0A0A0",
                         na_col = "#F5F5F5",
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         fontsize_row = 7,
                         fontsize_col = 7)

png("Normalized_WNT.png",
    width = 5000,
    height = 3000,
    res = 600)
plot_heatmap
dev.off()
