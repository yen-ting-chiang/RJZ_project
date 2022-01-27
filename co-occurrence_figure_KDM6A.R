library(dplyr)
setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/co-occurrence_figure")
getwd()
TCGA_mutation_KDM6A_genomic_co_occurrence = 
  read.csv(file = "TCGA_mutation_KDM6A_genomic_co_occurrence.csv")
TCGA_mutation_KDM6A_genomic_co_occurrence_filtered <- 
  TCGA_mutation_KDM6A_genomic_co_occurrence %>% 
  filter(Gene %in% 
           c("TP53", "BRCA2","POT1", "MEN1","VHL", "MSH2", "APC",
             "RECQL4","FAH") == TRUE) %>% 
  arrange(Log.Ratio) %>% 
  mutate(Gene = factor(Gene, levels = unique(Gene))) %>% 
  mutate(Log.Ratio = as.numeric(Log.Ratio))
library(ggplot2)
library(ggthemes)
theme_Publication <- function(base_size=14, base_family="Arial") {
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
            legend.key.size= unit(0.9, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

ggplot(data = TCGA_mutation_KDM6A_genomic_co_occurrence_filtered, 
       aes(x = Gene, y = Log.Ratio, fill = q.Value)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = terrain.colors(2),
                       trans = "log10") + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))
