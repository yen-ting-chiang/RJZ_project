setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_allele_frequency_figure")
getwd()

combo_data = read.csv(file = "gene_list_combo_data.csv")

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

combo_data$filename <- 
  factor(combo_data$filename, 
         levels = c("LFS-OBD7-MSC.anno_YTC.csv", 
                    "LFS-T2-MSC.anno_YTC.csv", 
                    "LFS-T4-MSC.anno_YTC.csv", 
                    "LFS-T5-MSC.anno_YTC.csv"))
levels(combo_data$filename)

ggplot(data = combo_data) + 
  geom_density(aes(x = tumor, 
                   fill = filename),
               alpha = 0.8) +
  scale_colour_Publication() + 
  scale_fill_Publication() + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))

ggplot(data = combo_data) + 
  geom_histogram(aes(x = tumor, 
                     fill = filename),
                 alpha = 0.8) + 
  facet_grid(cols = vars(filename)) + 
  scale_colour_Publication() + 
  scale_fill_Publication() + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))

ggplot(data = combo_data) + 
  geom_violin(aes(x = filename, 
                  y = tumor,
                  fill = filename),
              alpha = 0.8) + 
  scale_colour_Publication() + 
  scale_fill_Publication() + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))

ggplot(combo_data, aes(x= filename, y=tumor,
                       fill = filename)) + 
  geom_boxplot(alpha = 0.8) + 
  xlab("filename")+ 
  scale_colour_Publication() + 
  scale_fill_Publication() + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))


library(dplyr)
combo_data %>% 
  group_by(filename) %>% 
  summarise(median(tumor))
