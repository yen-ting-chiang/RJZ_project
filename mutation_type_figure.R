library(dplyr)
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_mutation_type_figure")
getwd()

filenames  <- list.files(pattern = "-MSC.anno_YTC.csv")

combo_data_no_filter <- purrr::map_df(filenames,
                            ~read.csv(.x, stringsAsFactors = FALSE,
                                      colClasses = "character") %>% 
                              mutate(filename = .x))

write.csv(combo_data_no_filter, 
          file = "gene_list_combo_data_no_filter.csv")
combo_data_no_filter = read.csv(file = "gene_list_combo_data_no_filter.csv")

library(ggplot2)

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}


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
