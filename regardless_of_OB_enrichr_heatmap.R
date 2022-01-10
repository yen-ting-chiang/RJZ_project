#Files for enrichr heatmap generation (regardless of OB):-------------------

setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/regardless_of_OB_enrichr_heatmap")
getwd()

library(dplyr)
library("enrichR")

listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

# all_dbs_list <- dbs[,3]
# dbs_list <- all_dbs_list

dbs_list <- c("GO_Biological_Process_2021")
venny_result <- read.csv(file = "venny_result.csv")

if (websiteLive) {
  enriched_T2 <- enrichr(venny_result[,1], dbs_list)
}
save(enriched_T2,
     file = "enriched_T2.rdata")
if (websiteLive) {
  enriched_T4 <- enrichr(venny_result[,2], dbs_list)
}
save(enriched_T4,
        file = "enriched_T4.rdata")
if (websiteLive) {
  enriched_T5 <- enrichr(venny_result[,3], dbs_list)
}
saveRDS(enriched_T5,
        file = "enriched_T5.rdata")
if (websiteLive) {
  enriched_at_least_2 <- enrichr(venny_result[,4], dbs_list)
}
save(enriched_at_least_2,
        file = "enriched_at_least_2.rdata")

# load("enriched_T2.rdata")
# load("enriched_T4.rdata")
# load("enriched_T5.rdata")
# load("enriched_at_least_2.rdata")

T2 <- enriched_T2[["GO_Biological_Process_2021"]]%>% 
  mutate(log_P.value = -log10(P.value)) %>% 
  arrange(desc(log_P.value))
write.csv(T2, file = "T2_GO_Biological_Process_2021.csv")

T4 <- enriched_T4[["GO_Biological_Process_2021"]]%>% 
  mutate(log_P.value = -log10(P.value)) %>% 
  arrange(desc(log_P.value))
write.csv(T4, file = "T4_GO_Biological_Process_2021.csv")

T5 <- enriched_T5[["GO_Biological_Process_2021"]]%>% 
  mutate(log_P.value = -log10(P.value)) %>% 
  arrange(desc(log_P.value))
write.csv(T5, file = "T5_GO_Biological_Process_2021.csv")

at_least_2 <- enriched_at_least_2[["GO_Biological_Process_2021"]]%>% 
  mutate(log_P.value = -log10(P.value)) %>% 
  arrange(desc(log_P.value))
write.csv(at_least_2, file = "at_least_2_GO_Biological_Process_2021.csv")


T2_GO_BP <- read.csv("T2_GO_Biological_Process_2021.csv")
T4_GO_BP <- read.csv("T4_GO_Biological_Process_2021.csv")
T5_GO_BP <- read.csv("T5_GO_Biological_Process_2021.csv")
at_least_2_GO_BP <- read.csv("at_least_2_GO_Biological_Process_2021.csv")


join_1 = full_join(T2_GO_BP,
                   T4_GO_BP,
                   by = c("Term" = "Term"),
                   suffix = c(".T2", ".T4"),
                   keep = FALSE)
join_2 = full_join(T5_GO_BP,
                   at_least_2_GO_BP,
                   by = c("Term" = "Term"),
                   suffix = c(".T5", ".at_least_2"),
                   keep = FALSE)

join_3 = full_join(join_1, join_2, 
                  by = c("Term" = "Term"),
                  keep = FALSE)

write.csv(join_3, file = "for_enrichr_GO_BP_clustering_all.csv")

join_4 = join_3 %>% select(Term, starts_with("log_P.value")) %>% 
  arrange(desc(log_P.value.at_least_2))
write.csv(join_4, file = "for_enrichr_GO_BP_clustering_log_pvalue_at_least_2_arranged.csv")




#for creating BioPlanet_2019 figure (selected term)-------------------------

BioPlanet_2019_table_selected <- 
  read.csv("BioPlanet_2019_table_selected.csv")

library("enrichR")
listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

if (websiteLive) plotEnrich(BioPlanet_2019_table_selected, 
                            showTerms = 21, 
                            numChar = 60, 
                            y = "Count", 
                            orderBy = "Term")
