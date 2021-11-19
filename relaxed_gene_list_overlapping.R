library(dplyr)
setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/relaxed_gene_list")
getwd()


filenames  <- list.files(pattern = "*-MSC.anno_YTC.csv")

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

gene_number = gene_list %>% 
  group_by(Gene) %>%
  summarise(n())

write.csv(gene_number, file = "gene_number.csv")

combo_data <- purrr::map_df(filenames,
                            ~read.csv(.x, stringsAsFactors = FALSE,
                                      colClasses = "character") %>%
                              filter(ExonicFunction %in% c("nonframeshift_substitution",
                                                           "frameshift_substitution",
                                                           "nonsynonymous_SNV",
                                                           "stopgain",
                                                           "stoploss")) %>% 
                              mutate(filename = .x))

write.csv(combo_data, 
          file = "gene_list_combo_data.csv")

#need to be corrected
overlap_3 = combo_data %>% 
  filter(Gene %in% c("EPC2","FOXJ2","GABRR2","KCNMA1","LCORL",
                     "LRP1B","MATR3","MUC6","PCDH7","USH2A","ZFHX3"))

write.csv(overlap_3, 
          file = "overlap_3_combo_data.csv")

