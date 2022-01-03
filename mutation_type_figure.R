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
           position = "stack",
           alpha = 0.8) + 
  scale_fill_manual(values = c("#a6cee3","#1f78b4","#b2df8a",
                               "#33a02c","#fb9a99","#e31a1c",
                               "#fdbf6f","#ff7f00","#cab2d6",
                               "#6a3d9a","#ffff99","#b15928",
                               "#543005")) + 
  scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a",
                                "#33a02c","#fb9a99","#e31a1c",
                                "#fdbf6f","#ff7f00","#cab2d6",
                                "#6a3d9a","#ffff99","#b15928",
                                "#543005")) +
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))

ggplot(data = combo_data_no_filter) + 
  geom_bar(aes(x = filename, 
               fill = RefGene),
           position = "fill",
           alpha = 0.8) + 
  scale_fill_manual(values = c("#a6cee3","#1f78b4","#b2df8a",
                               "#33a02c","#fb9a99","#e31a1c",
                               "#fdbf6f","#ff7f00","#cab2d6",
                               "#6a3d9a","#ffff99","#b15928",
                               "#543005")) + 
  scale_color_manual(values = c("#a6cee3","#1f78b4","#b2df8a",
                                "#33a02c","#fb9a99","#e31a1c",
                                "#fdbf6f","#ff7f00","#cab2d6",
                                "#6a3d9a","#ffff99","#b15928",
                                "#543005")) +
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))
