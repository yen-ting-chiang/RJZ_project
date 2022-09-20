setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/TP53_co-occurrence_table")
getwd()
library(dplyr)
cooccurrence_TP53 <- 
  read.csv(file = "cooccurrence_TP53.csv")
allele_frequency <- 
  read.csv(file = "T2_T4_T5_allele_frequency_sum.csv")

top_30 <- allele_frequency[-16,] %>% 
  head(30)
gene_order <- top_30$Gene 
cooccurrence_TP53 <- cooccurrence_TP53 %>% 
  arrange(factor(B, levels = gene_order))
cooccurrence_TP53_for_heatmap <- 
  cooccurrence_TP53 %>% 
  select(c(2,7))

library(pheatmap)
library(RColorBrewer)

cooccurrence_TP53_for_heatmap_tmp = 
  cooccurrence_TP53_for_heatmap[,c(2)]
cooccurrence_TP53_for_heatmap_tmp <- 
  as.data.frame(cooccurrence_TP53_for_heatmap_tmp)
row.names(cooccurrence_TP53_for_heatmap_tmp) = 
  cooccurrence_TP53_for_heatmap[,1]

cooccurrence_TP53_for_heatmap_tmp[27,1] <- NA
cooccurrence_TP53_for_heatmap_tmp[7,1] <- NA

cooccurrence_TP53_for_heatmap_tmp <- 
  log10(cooccurrence_TP53_for_heatmap_tmp)
pheatmap(cooccurrence_TP53_for_heatmap_tmp,
         color = colorRampPalette(brewer.pal(n =9, 
                                             name = "Greens"))(100),
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color = "#A0A0A0",
         na_col = "#E0E0E0",
         show_rownames = TRUE,
         show_colnames = FALSE)
