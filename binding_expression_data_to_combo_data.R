library(dplyr)
setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/at_least_two_tumor")
getwd()


at_least_two_tumor_combo_data <- 
  read.csv("at_least_two_tumor_combo_data.csv")

expression_data <- 
  read.csv("OB FKPM more than 1 gene list.csv")


at_least_two_tumor_combo_data_combining_expression <- 
  left_join(at_least_two_tumor_combo_data, 
          expression_data, 
          by = c("Gene" = "NAME"))
write.csv(at_least_two_tumor_combo_data_combining_expression, 
          file = "at_least_two_tumor_combo_data_combining_expression_all.csv")

