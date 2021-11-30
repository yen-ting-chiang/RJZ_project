library(dplyr)
setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/relaxed_gene_list")
getwd()


filenames  <- list.files(pattern = "*-MSC.anno_YTC.csv")

#combine files, dintinct gene------------------------------------
gene_list <- purrr::map_df(filenames,
                           ~read.csv(.x, stringsAsFactors = FALSE,
                                     colClasses = "character") %>%
                             filter(ExonicFunction %in% c("nonframeshift_substitution",
                                                          "frameshift_substitution",
                                                          "nonsynonymous_SNV",
                                                          "stopgain",
                                                          "stoploss")) %>% 
                             distinct(Gene, .keep_all = FALSE) %>% 
                             mutate(filename = .x))


#create anyone mutation gene list (containing OB)------------------
gene_anyone_single_name <- gene_list %>% 
  distinct(Gene, .keep_all = TRUE) %>% 
  dplyr::filter(!grepl(",",Gene))
write.csv(gene_anyone_single_name, 
          file = "gene_anyone_single_name.csv")

gene_anyone_multiple_name <- gene_list %>% 
  distinct(Gene, .keep_all = TRUE) %>% 
  dplyr::filter(grepl(",",Gene))
write.csv(gene_anyone_multiple_name, 
          file = "gene_anyone_multiple_name.csv")


#create anyone mutation gene list (not containing OB)------------------
gene_anyone_not_OB_single_name <- gene_list %>%
  filter(filename != "LFS-OBD7-MSC.anno_YTC.csv") %>% 
  distinct(Gene, .keep_all = TRUE) %>% 
  dplyr::filter(!grepl(",",Gene))
write.csv(gene_anyone_not_OB_single_name, 
          file = "gene_anyone_not_OB_single_name.csv")

gene_anyone_not_OB_multiple_name <- gene_list %>%
  filter(filename != "LFS-OBD7-MSC.anno_YTC.csv") %>% 
  distinct(Gene, .keep_all = TRUE) %>% 
  dplyr::filter(grepl(",",Gene))
write.csv(gene_anyone_not_OB_multiple_name, 
          file = "gene_anyone_not_OB_multiple_name.csv")



#create OB, T2, T4, and T5 mutation gene list -------------------
gene_list_OB = gene_list %>% 
  filter(filename == "LFS-OBD7-MSC.anno_YTC.csv")
gene_OB_single_name <- gene_list_OB %>% 
  dplyr::filter(!grepl(",",Gene))
write.csv(gene_OB_single_name, 
          file = "gene_OB_single_name.csv")
gene_OB_multiple_name <- gene_list_OB %>% 
  dplyr::filter(grepl(",",Gene))
write.csv(gene_OB_multiple_name, 
          file = "gene_OB_multiple_name.csv")


gene_list_T2 = gene_list %>% 
  filter(filename == "LFS-T2-MSC.anno_YTC.csv")
gene_T2_single_name <- gene_list_T2 %>% 
  dplyr::filter(!grepl(",",Gene))
write.csv(gene_T2_single_name, 
          file = "gene_T2_single_name.csv")
gene_T2_multiple_name <- gene_list_T2 %>% 
  dplyr::filter(grepl(",",Gene))
write.csv(gene_T2_multiple_name, 
          file = "gene_T2_multiple_name.csv")


gene_list_T4 = gene_list %>% 
  filter(filename == "LFS-T4-MSC.anno_YTC.csv")
gene_T4_single_name <- gene_list_T4 %>% 
  dplyr::filter(!grepl(",",Gene))
write.csv(gene_T4_single_name, 
          file = "gene_T4_single_name.csv")
gene_T4_multiple_name <- gene_list_T4 %>% 
  dplyr::filter(grepl(",",Gene))
write.csv(gene_T4_multiple_name, 
          file = "gene_T4_multiple_name.csv")



gene_list_T5 = gene_list %>% 
  filter(filename == "LFS-T5-MSC.anno_YTC.csv")
gene_T5_single_name <- gene_list_T5 %>% 
  dplyr::filter(!grepl(",",Gene))
write.csv(gene_T5_single_name, 
          file = "gene_T5_single_name.csv")
gene_T5_multiple_name <- gene_list_T5 %>% 
  dplyr::filter(grepl(",",Gene))
write.csv(gene_T5_multiple_name, 
          file = "gene_T5_multiple_name.csv")


# Summarize the gene occurrence number----------------------------
gene_number = gene_list %>% 
  group_by(Gene) %>%
  summarise(n())
write.csv(gene_number, file = "gene_number.csv")
gene_number <- read.csv(file = "gene_number.csv")

# read manual generated at_least_2 mutation gene list------------
gene_number_at_least_2 <- 
  read.csv(file = "gene_number_at_least_2.csv")
gene_number_at_least_2_list <- gene_number_at_least_2[,2]


# create anytwo mutation gene list------------------------------
gene_anytwo_single_name <- gene_number_at_least_2 %>% 
  dplyr::filter(!grepl(",",Gene))
write.csv(gene_anytwo_single_name, 
          file = "gene_anytwo_single_name.csv")

gene_anytwo_multiple_name <- gene_number_at_least_2 %>% 
  dplyr::filter(grepl(",",Gene))
write.csv(gene_anytwo_multiple_name, 
          file = "gene_anytwo_multiple_name.csv")

gene_anytwo_not_OB_single_name <- gene_number_at_least_2 %>% 
  dplyr::filter(!grepl(",",Gene)) %>% 
  filter(Gene %in% gene_list_OB[,1] == FALSE)
write.csv(gene_anytwo_not_OB_single_name, 
          file = "gene_anytwo_not_OB_single_name.csv")

gene_anytwo_not_OB_multiple_name <- gene_number_at_least_2 %>% 
  dplyr::filter(grepl(",",Gene)) %>% 
  filter(Gene %in% gene_list_OB[,1] == FALSE)
write.csv(gene_anytwo_not_OB_multiple_name, 
          file = "gene_anytwo_not_OB_multiple_name.csv")


#Files for enrichr heatmap generation:

# gene_anyone_single_name = read.csv(file = "gene_anyone_single_name.csv")
# gene_anyone_not_OB_single_name <- read.csv(file = "gene_anyone_not_OB_single_name.csv")
# gene_OB_single_name <- read.csv(file = "gene_OB_single_name.csv")
# gene_T2_single_name <- read.csv(file = "gene_T2_single_name.csv")
# gene_T4_single_name <- read.csv(file = "gene_T4_single_name.csv")
# gene_T5_single_name <- read.csv(file = "gene_T5_single_name.csv")
# gene_anytwo_single_name <- read.csv(file = "gene_anytwo_single_name.csv")
# gene_anytwo_not_OB_single_name <- read.csv(file = "gene_anytwo_not_OB_single_name.csv")





library("enrichR")

listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes

websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

# all_dbs_list <- dbs[,3]
# dbs_list <- all_dbs_list

dbs_list <- c("GO_Biological_Process_2021", 
              "GO_Molecular_Function_2021", 
              "BioPlanet_2019",
              "WikiPathway_2021_Human",
              "KEGG_2021_Human",
              "MSigDB_Hallmark_2020",
              "Reactome_2016",
              "ChEA_2016",
              "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")


input_gene_list <- gene_anytwo_single_name[,3]
if (websiteLive) {
  enriched <- enrichr(input_gene_list, dbs_list)
}
save(enriched,
     file = "enriched_gene_anytwo_single_name.rdata")


load("enriched_gene_OB_single_name.rdata")
enriched_gene_OB_single_name <- enriched
OB <- enriched_gene_OB_single_name[["KEGG_2021_Human"]]%>% 
  mutate(log_P.value = -log10(P.value)) %>% 
  arrange(desc(log_P.value))
write.csv(OB, file = "gene_OB_single_name_KEGG.csv")

load("enriched_gene_T2_single_name.rdata")
enriched_gene_T2_single_name <- enriched
T2 <- enriched_gene_T2_single_name[["KEGG_2021_Human"]]%>% 
  mutate(log_P.value = -log10(P.value)) %>% 
  arrange(desc(log_P.value))
write.csv(T2, file = "gene_T2_single_name_KEGG.csv")

load("enriched_gene_T4_single_name.rdata")
enriched_gene_T4_single_name <- enriched
T4 <- enriched_gene_T4_single_name[["KEGG_2021_Human"]]%>% 
  mutate(log_P.value = -log10(P.value)) %>% 
  arrange(desc(log_P.value))
write.csv(T4, file = "gene_T4_single_name_KEGG.csv")

load("enriched_gene_T5_single_name.rdata")
enriched_gene_T5_single_name <- enriched
T5 <- enriched_gene_T5_single_name[["KEGG_2021_Human"]]%>% 
  mutate(log_P.value = -log10(P.value)) %>% 
  arrange(desc(log_P.value))
write.csv(T5, file = "gene_T5_single_name_KEGG.csv")

load("enriched_gene_anytwo_not_OB_single_name.rdata")
enriched_gene_anytwo_not_OB_single_name <- enriched
anytwo_not_OB <- enriched_gene_anytwo_not_OB_single_name[["KEGG_2021_Human"]]%>% 
  mutate(log_P.value = -log10(P.value)) %>% 
  arrange(desc(log_P.value))
write.csv(anytwo_not_OB, file = "gene_anytwo_not_OB_single_name_KEGG.csv")

load("enriched_gene_anytwo_single_name.rdata")
enriched_gene_anytwo_single_name <- enriched
anytwo <- enriched_gene_anytwo_single_name[["KEGG_2021_Human"]]%>% 
  mutate(log_P.value = -log10(P.value)) %>% 
  arrange(desc(log_P.value))
write.csv(anytwo, file = "gene_anytwo_single_name_KEGG.csv")


# a <- enriched[["GO_Biological_Process_2021"]]
# b <- a %>% 
#   mutate(log_P.value = -log10(P.value)) %>% 
#   arrange(desc(log_P.value))
# write.csv(b, file = "gene_anytwo_single_name_GO_BP.csv")


gene_OB_single_name_KEGG <- read.csv("gene_OB_single_name_KEGG.csv")
gene_T2_single_name_KEGG <- read.csv("gene_T2_single_name_KEGG.csv")
gene_T4_single_name_KEGG <- read.csv("gene_T4_single_name_KEGG.csv")
gene_T5_single_name_KEGG <- read.csv("gene_T5_single_name_KEGG.csv")
gene_anytwo_not_OB_single_name_KEGG <- read.csv("gene_anytwo_not_OB_single_name_KEGG.csv")
gene_anytwo_single_name_KEGG <- read.csv("gene_anytwo_single_name_KEGG.csv")


join_1 = full_join(gene_OB_single_name_KEGG,
                   gene_T2_single_name_KEGG,
                            by = c("Term" = "Term"),
                            suffix = c(".OB", ".T2"),
                            keep = FALSE)
join_2 = full_join(gene_T4_single_name_KEGG,
                   gene_T5_single_name_KEGG,
                   by = c("Term" = "Term"),
                   suffix = c(".T4", ".T5"),
                   keep = FALSE)
join_3 = full_join(gene_anytwo_not_OB_single_name_KEGG,
                   gene_anytwo_single_name_KEGG,
                   by = c("Term" = "Term"),
                   suffix = c(".2-OB", ".two"),
                   keep = FALSE)
join4 = full_join(join_1, join_2, 
                  by = c("Term" = "Term"),
                  keep = FALSE)
join5 = full_join(join4, join_3, 
                  by = c("Term" = "Term"),
                  keep = FALSE)
write.csv(join5, file = "for_enrichr_KEGG_clustering_all.csv")

join6 = join5 %>% select(Term, starts_with("log_P.value")) %>% 
  arrange(desc(log_P.value.two))
write.csv(join6, file = "for_enrichr_KEGG_clustering_log_pvalue_anytwo.csv")

if (websiteLive) plotEnrich(gene_anytwo_single_name_KEGG, 
                            showTerms = 30, 
                            numChar = 100, 
                            y = "Count", 
                            orderBy = "P.value")

library(dplyr)
bind_all_result = enriched %>% 
  bind_rows() %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  arrange(desc(Combined.Score))

write.csv(bind_all_result, file = "LFS_OS_any2or3_enrichr_bind_all_result.csv")



