library(dplyr)
setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/co-occurrence_figure/co-occurrence")
getwd()
TCGA_TP53_OB_top_15_genes = 
  read.csv(file = "TCGA_TP53_OB_top_15_genes.csv")
TCGA_TP53_OB_top_15_genes[,7] =
  gsub('>', '', TCGA_TP53_OB_top_15_genes[,7])
TCGA_TP53_OB_top_15_genes[,8] =
  gsub('<', '', TCGA_TP53_OB_top_15_genes[,8])
TCGA_TP53_OB_top_15_genes[,9] =
  gsub('<', '', TCGA_TP53_OB_top_15_genes[,9])
TCGA_TP53_OB_top_15_genes_filtered <- 
  TCGA_TP53_OB_top_15_genes %>% 
  filter(A == "TP53") %>% 
  mutate(Log2.Odds.Ratio = as.numeric(Log2.Odds.Ratio)) %>% 
  mutate(p.Value = as.numeric(p.Value)) %>% 
  mutate(q.Value = as.numeric(q.Value)) %>% 
  arrange(Log2.Odds.Ratio) %>% 
  mutate(B = factor(B, levels = unique(B)))
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

ggplot(data = TCGA_TP53_OB_top_15_genes_filtered, 
       aes(x = B, y = Log2.Odds.Ratio, fill = q.Value)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = terrain.colors(10),
                       trans = "log10") + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))
