setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/combine_WES")
getwd()


filenames  <- list.files(pattern = "-MSC.anno_YTC.csv")
library(dplyr)
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

gene_number <- gene_list %>% 
  group_by(Gene) %>%
  summarise(n())
colnames(gene_number) <- c("Gene", "number")
gene_number <- gene_number %>% 
  arrange(desc(number))

write.csv(gene_number, file = "gene_number.csv")

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
