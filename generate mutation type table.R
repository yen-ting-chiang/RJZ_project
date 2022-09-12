setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/mutation_type_table")
getwd()
library(dplyr)

gene_list_combo_data <- read.csv("gene_list_combo_data.csv")
top_30 <- gene_list_combo_data %>% 
  filter(filename != "LFS-OBD7-MSC.anno_YTC.csv") %>% 
  group_by(Gene) %>% 
  summarise(sum(tumor))
colnames(top_30) <- c("Gene", "allele_frequency")
top_30 <- top_30 %>% 
  arrange(desc(allele_frequency)) %>% 
  head(30)

mutation_type <- gene_list_combo_data %>% 
  filter(filename != "LFS-OBD7-MSC.anno_YTC.csv") %>% 
  group_by(Gene, filename, ExonicFunction) %>% 
  summarise(n())

colnames(mutation_type) <- c("Gene", "filename", "ExonicFunction", "number")

mutation_type_top_30 <- mutation_type %>% 
  filter(Gene %in% top_30$Gene == TRUE)

gene_order <- top_30$Gene 
mutation_type_top_30 <- mutation_type_top_30 %>% 
  arrange(factor(Gene, levels = gene_order))
write.csv(mutation_type_top_30, file = "mutation_type_top_30.csv")
