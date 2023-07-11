setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/Y2_KD_genes_heatmap")
getwd()
library(data.table)
library(dplyr)
counts_file <- fread(file = "Normalized_counts_file.csv")
cell_cycle_genes <- fread(file = "KEGG_CELL_CYCLE.csv")
counts_file_cell_cycle <- 
  counts_file %>% 
  filter(gene_name %in% cell_cycle_genes$SYMBOL == TRUE) %>% 
  filter(ESNG != "ENSG00000269846")
library(pheatmap)
library(RColorBrewer)
counts_file_cell_cycle <- counts_file_cell_cycle[,2:8]
rownames(counts_file_cell_cycle) <- counts_file_cell_cycle[, 7]
counts_file_cell_cycle <- counts_file_cell_cycle[, -7]
counts_file_cell_cycle <- 
  tibble::column_to_rownames(counts_file_cell_cycle, counts_file_cell_cycle[,7])
write.csv(counts_file_cell_cycle, "counts_file_cell_cycle.csv")
counts_file_cell_cycle <- fread("counts_file_cell_cycle.csv")
counts_file_cell_cycle <- 
  tibble::column_to_rownames(counts_file_cell_cycle, counts_file_cell_cycle[,1])


cell_cycle_normalized <- read.csv(file = "cell_cycle_normalized.csv",
                                  header = TRUE)

rownames(cell_cycle_normalized) <- cell_cycle_normalized[, 1]
cell_cycle_normalized <- cell_cycle_normalized[, -1]

cell_cycle_normalized <- cell_cycle_normalized[, -1]

plot_heatmap <- pheatmap(cell_cycle_normalized,
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

png("data.png",
    width = 5000,
    height = 3000,
    res = 600)
plot_heatmap
dev.off()
