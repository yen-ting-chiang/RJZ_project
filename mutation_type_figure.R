library(dplyr)
setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/at_least_two_tumor")
getwd()

filenames  <- list.files(pattern = "*-MSC.anno_YTC.csv")

combo_data_no_filter <- purrr::map_df(filenames,
                            ~read.csv(.x, stringsAsFactors = FALSE,
                                      colClasses = "character") %>% 
                              mutate(filename = .x))

setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/RJZ_ggplot2")
getwd()
write.csv(combo_data_no_filter, 
          file = "gene_list_combo_data_no_filter.csv")
combo_data_no_filter = read.csv(file = "gene_list_combo_data_no_filter.csv")

library(ggplot2)
ggplot(data = combo_data_no_filter) + 
  geom_bar(aes(x = filename, 
               fill = RefGene),
           position = "stack")
