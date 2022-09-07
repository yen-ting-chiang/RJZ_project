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
# combo_data = read.csv(file = "gene_list_combo_data.csv")

expression_data <- 
  read.csv("OB FKPM more than 1 gene list.csv")

combo_data_combining_expression <- 
  left_join(combo_data, 
            expression_data, 
            by = c("Gene" = "NAME"))
write.csv(combo_data_combining_expression, 
          file = "combo_data_combining_expression.csv")

combo_data_combining_expression_D0_D24 <- 
  left_join(combo_data, 
            expression_data[,-(3:46)], 
            by = c("Gene" = "NAME"))
write.csv(combo_data_combining_expression_D0_D24, 
          file = "combo_data_combining_expression_D0_D24.csv")

combo_data = read.csv(file = "gene_list_combo_data.csv")

combo_data_filtered <- combo_data %>% 
  filter(filename != "LFS-OBD7-MSC.anno_YTC.csv") %>% 
  dplyr::filter(!grepl(",",Gene)) %>% 
  group_by(Gene) %>% 
  summarise(allele_frequency_sum = sum(tumor))
colnames(combo_data_filtered) <- c("Gene", "number")
combo_data_filtered <- combo_data_filtered %>% 
  arrange(desc(number))
write.csv(combo_data_filtered, file = "T2_T4_T5_allele_frequency_sum.csv")


