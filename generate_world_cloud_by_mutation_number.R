setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/word_cloud")
getwd()
library(dplyr)

gene_list_combo_data <- read.csv("gene_list_combo_data.csv")

T2 <- gene_list_combo_data %>% 
  filter(filename == "LFS-T2-MSC.anno_YTC.csv") %>% 
  group_by(Gene) %>% 
  summarise(n())
colnames(T2) <- c("Gene", "number")
T2 <- T2 %>% 
  arrange(desc(number))
T2_top30 <- 
  T2 %>% head(30)
#add YTHDF2 to T2----------------------------
YTHDF2_T2_num <- T2 %>% 
  filter(Gene == "YTHDF2")
T2_top30_add_Y2 <- 
  rbind(T2_top30, YTHDF2_T2_num)
wordcloud2(data = T2_top30_add_Y2, 
           size = 0.3,
           fontFamily = 'Arial',
           shuffle = 'FALSE',
           color = 'random-dark',
           shape = 'circle')

T4 <- gene_list_combo_data %>% 
  filter(filename == "LFS-T4-MSC.anno_YTC.csv") %>% 
  group_by(Gene) %>% 
  summarise(n())
colnames(T4) <- c("Gene", "number")
T4 <- T4 %>% 
  arrange(desc(number))
T4_top30 <- 
  T4 %>% head(30)
wordcloud2(data = T4_top30, 
           size = 0.3,
           fontFamily = 'Arial',
           shuffle = 'FALSE',
           color = 'random-dark',
           shape = 'circle')

T5 <- gene_list_combo_data %>% 
  filter(filename == "LFS-T5-MSC.anno_YTC.csv") %>% 
  group_by(Gene) %>% 
  summarise(n())
colnames(T5) <- c("Gene", "number")
T5 <- T5 %>% 
  arrange(desc(number))
T5_top30 <- 
  T5 %>% head(30)
#add YTHDF2 to T5----------------------------
YTHDF2_T5_num <- T5 %>% 
  filter(Gene == "YTHDF2")
T5_top30_add_Y2 <- 
  rbind(T5_top30, YTHDF2_T5_num)
wordcloud2(data = T5_top30_add_Y2, 
           size = 0.3,
           fontFamily = 'Arial',
           shuffle = 'FALSE',
           color = 'random-dark',
           shape = 'circle')


T_all <- gene_list_combo_data %>% 
  group_by(Gene) %>% 
  summarise(n())
colnames(T_all) <- c("Gene", "number")
T_all <- T_all %>% 
  arrange(desc(number))
T_all_top30 <- 
  T_all %>% head(30)
#add YTHDF2 to T_all----------------------------
YTHDF2_T_all_num <- T_all %>% 
  filter(Gene == "YTHDF2")
T_all_top30_add_Y2 <- 
  rbind(T_all_top30, YTHDF2_T_all_num)
wordcloud2(data = T_all_top30_add_Y2, 
           size = 0.3,
           fontFamily = 'Arial',
           shuffle = 'FALSE',
           color = 'random-dark',
           shape = 'circle')

