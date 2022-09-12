library(dplyr)
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/enrichr_2")
getwd()

T2 <- read.csv("T2_allele_frequency_sum.csv")
T4 <- read.csv("T4_allele_frequency_sum.csv")
T5 <- read.csv("T5_allele_frequency_sum.csv")
T_all <- read.csv("T2_T4_T5_allele_frequency_sum.csv")

T2_0.1 <- T2 %>% 
  filter(allele_frequency_sum >= 0.1)
T4_0.1 <- T4 %>% 
  filter(allele_frequency_sum >= 0.1)
T5_0.1 <- T5 %>% 
  filter(allele_frequency_sum >= 0.1)
T_all_0.1 <- T_all %>% 
  filter(number >= 0.1)

library("enrichR")
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

# all_dbs_list <- dbs[,3]
# dbs_list <- all_dbs_list

# dbs_list <- c("GO_Biological_Process_2021", 
#               "GO_Molecular_Function_2021", 
#               "BioPlanet_2019",
#               "WikiPathway_2021_Human",
#               "KEGG_2021_Human",
#               "MSigDB_Hallmark_2020",
#               "Reactome_2016",
#               "ChEA_2016",
#               "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")


db <- c("GO_Biological_Process_2021")
T2_list <- T2[,2]
if (websiteLive) {
  T2_enrichr <- enrichr(T2_list, db)
}
T2_df <- T2_enrichr$GO_Biological_Process_2021

db <- c("GO_Biological_Process_2021")
T4_list <- T4[,2]
if (websiteLive) {
  T4_enrichr <- enrichr(T4_list, db)
}
T4_df <- T4_enrichr$GO_Biological_Process_2021

db <- c("GO_Biological_Process_2021")
T5_list <- T5[,2]
if (websiteLive) {
  T5_enrichr <- enrichr(T5_list, db)
}
T5_df <- T5_enrichr$GO_Biological_Process_2021

db <- c("GO_Biological_Process_2021")
T_all_list <- T_all[,2]
if (websiteLive) {
  T_all_enrichr <- enrichr(T_all_list, db)
}
T_all_df <- T_all_enrichr$GO_Biological_Process_2021



db <- c("GO_Biological_Process_2021")

T2_0.1_list <- T2_0.1[,2]
if (websiteLive) {
  T2_0.1_enrichr <- enrichr(T2_0.1_list, db)
}
T2_0.1_df <- T2_0.1_enrichr$GO_Biological_Process_2021

T4_0.1_list <- T4_0.1[,2]
if (websiteLive) {
  T4_0.1_enrichr <- enrichr(T4_0.1_list, db)
}
T4_0.1_df <- T4_0.1_enrichr$GO_Biological_Process_2021

T5_0.1_list <- T5_0.1[,2]
if (websiteLive) {
  T5_0.1_enrichr <- enrichr(T5_0.1_list, db)
}
T5_0.1_df <- T5_0.1_enrichr$GO_Biological_Process_2021

T_all_0.1_list <- T_all_0.1[,2]
if (websiteLive) {
  T_all_0.1_enrichr <- enrichr(T_all_0.1_list, db)
}
T_all_0.1_df <- T_all_0.1_enrichr$GO_Biological_Process_2021



T2_df_p_filtered <- T2_df %>% 
  filter(P.value < 0.05)
T4_df_p_filtered <- T4_df %>% 
  filter(P.value < 0.05)
T5_df_p_filtered <- T5_df %>% 
  filter(P.value < 0.05)
T_all_p_filtered <- T_all_df %>% 
  filter(P.value < 0.05)

T2_0.1_df_p_filtered <- T2_0.1_df %>% 
  filter(P.value < 0.05) %>% 
  filter(Adjusted.P.value < 0.25)
T4_0.1_df_p_filtered <- T4_0.1_df %>% 
  filter(P.value < 0.05) %>% 
  filter(Adjusted.P.value < 0.25)
T5_0.1_df_p_filtered <- T5_0.1_df %>% 
  filter(P.value < 0.05) %>% 
  filter(Adjusted.P.value < 0.25)
T_all_0.1_p_filtered <- T_all_0.1_df %>% 
  filter(P.value < 0.05) %>% 
  filter(Adjusted.P.value < 0.25)


DFlist <- list(T2_df_p_filtered,
               T4_df_p_filtered,
               T5_df_p_filtered,
               T_all_p_filtered)

DFlist_0.1 <- list(T2_0.1_df_p_filtered,
               T4_0.1_df_p_filtered,
               T5_0.1_df_p_filtered,
               T_all_0.1_p_filtered)



multi_full_0.1 <- Reduce(
  function(...) {
    full_join(..., 
              by = c("Term" = "Term"),
              keep = FALSE)
  },
  DFlist_0.1
)


#pheatmap----------------------------------------------------------

for_heatmap <- multi_full %>% 
  select(Term, starts_with("Combined")) %>% 
  arrange(desc(Combined.Score.y.y)) %>% 
  filter(Combined.Score.y.y>10)

for_heatmap_0.1 <- multi_full_0.1 %>% 
  select(Term, starts_with(c("Combined",
                           "adjust", 
                           "p.value",
                           "Genes",
                           "Overlap"))) %>% 
  arrange(desc(Combined.Score.y.y)) %>% 
  filter(Combined.Score.y.y>10)
write.csv(for_heatmap_0.1, file = "for_heatmap_0.1.csv")


#manually select terms for heat map generation

library(pheatmap)
library(RColorBrewer)

for_heatmap_selected <- 
  read.csv(file = "for_heatmap_selected.csv", 
           header = T)


for_heatmap_selected_tmp = 
  for_heatmap_selected[,c(2:5)]
row.names(for_heatmap_selected_tmp) = 
  for_heatmap_selected[,1]

for_heatmap_selected_tmp <- 
  log10(for_heatmap_selected_tmp)
pheatmap(for_heatmap_selected_tmp,
         color = colorRampPalette(brewer.pal(n =9, 
                                                 name = "Reds"))(100),
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         gaps_row=(c(9,13,17,22,27)),
         border_color = "#A0A0A0",
         na_col = "#E0E0E0",
         show_rownames = FALSE)
