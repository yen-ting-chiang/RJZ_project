library(dplyr)
setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/at_least_two_tumor")
getwd()


filenames  <- list.files(pattern = "*-MSC.anno_YTC.csv")

gene_list <- purrr::map_df(filenames,
                           ~read.csv(.x, stringsAsFactors = FALSE,
                                     colClasses = "character") %>%
                             filter(ExonicFunction %in% 
                                      c("nonframeshift_substitution",
                                        "frameshift_substitution",
                                        "nonsynonymous_SNV",
                                        "stopgain",
                                        "stoploss")) %>% 
                             distinct(Gene, .keep_all = FALSE) %>% 
                             mutate(filename = .x))

gene_list_T2_T4_T5 <- gene_list %>% 
  filter(filename != "LFS-OBD7-MSC.anno_YTC.csv")

gene_number_T2_T4_T5 <- gene_list_T2_T4_T5 %>% 
  group_by(Gene) %>%
  summarise(n())
write.csv(gene_number_T2_T4_T5, file = "gene_number_T2_T4_T5.csv")
gene_number <- read.csv(file = "gene_number.csv")


combo_data <- purrr::map_df(filenames,
                            ~read.csv(.x, stringsAsFactors = FALSE,
                                      colClasses = "character") %>%
                              filter(ExonicFunction %in% 
                                       c("nonframeshift_substitution",
                                         "frameshift_substitution",
                                         "nonsynonymous_SNV",
                                         "stopgain",
                                         "stoploss")) %>% 
                              mutate(filename = .x))

write.csv(combo_data, 
          file = "gene_list_combo_data.csv")
combo_data = read.csv(file = "gene_list_combo_data.csv")

at_least_two_tumor_gene_list <- read.csv("at_least_two_gene_list.txt")


at_least_two_tumor_combo_data = combo_data %>% 
  filter(Gene %in% at_least_two_tumor_gene_list[,1])

write.csv(at_least_two_tumor_combo_data, 
          file = "at_least_two_tumor_combo_data.csv")


