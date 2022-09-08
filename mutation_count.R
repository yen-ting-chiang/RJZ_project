setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/YTHDF2_lollipop_plot")
getwd()
library(dplyr)
df <- read.csv(file = "at_least_two_tumor_combo_data_combining_expression_average.csv")
df2 <- df %>% 
  group_by(Gene) %>% 
  summarise(n())
colnames(df2) <- c("Gene", "number")
df2 <- df2 %>% 
  arrange(desc(number))
write.csv(df2, file = "mutation_count.csv")
