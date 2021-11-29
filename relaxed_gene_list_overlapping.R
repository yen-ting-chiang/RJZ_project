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



gene_number = gene_list %>% 
  group_by(Gene) %>%
  summarise(n())
write.csv(gene_number, file = "gene_number.csv")
gene_number <- read.csv(file = "gene_number.csv")

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


gene_number_at_least_2 <- 
  read.csv(file = "gene_number_at_least_2.csv")
gene_number_at_least_2_list <- gene_number_at_least_2[,2]


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
combo_data = read.csv(file = "gene_list_combo_data.csv")




overlap_3_and_4_gene_list <- c("FOXJ2", "KCNMA1", "ZFHX3", "ABCC1", "ANKRD17", "APOBEC3G", "ARID1B", "CHD3", "DAB2IP", "DHX37", "DLGAP3", "DLX6", "DUSP27", "EGR1", "EPC2", "GABRR2", "HERC2", "IGF1R", "KMT2D", "LCORL", "LILRB1", "LRP1B", "MATR3", "MUC16", "MUC6", "NFIX", "PCDH10", "PCDH19", "PCDH7", "PPM1E", "PRDM8", "RBM14", "REV3L", "RYR1", "SETDB1", "SFRP1", "TENM2", "TNRC6C", "TTN", "USH2A", "ZAK", "ZNF592")

overlap_3_and_4_combo_data = combo_data %>% 
  filter(Gene %in% overlap_3_and_4_gene_list)

write.csv(overlap_3_and_4_combo_data, 
          file = "overlap_3_and_4_combo_data.csv")

OB_contain_gene_data = overlap_3_and_4_combo_data %>% 
  filter(filename == "LFS-OBD7-MSC.anno_YTC.csv")
OB_contain_gene_list <- OB_contain_gene_data[,"Gene"]

overlap_3_and_4_combo_data_containing_OB = 
  overlap_3_and_4_combo_data %>% 
  filter(Gene %in% OB_contain_gene_list == TRUE)

write.csv(overlap_3_and_4_combo_data_containing_OB,
          file = "overlap_3_and_4_combo_data_containing_OB.csv")

overlap_3_and_4_combo_data_not_containing_OB = 
  overlap_3_and_4_combo_data %>% 
  filter(Gene %in% OB_contain_gene_list == FALSE)

write.csv(overlap_3_and_4_combo_data_not_containing_OB,
          file = "overlap_3_and_4_combo_data_not_containing_OB.csv")
