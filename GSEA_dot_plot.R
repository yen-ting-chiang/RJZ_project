setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/GSEA_dot_plot")
getwd()
library(dplyr)
library(data.table)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(tidyverse)

neg_df <- fread(file = "KEGG_C6_neg.csv")
neg_KEGG <- neg_df %>% 
  filter(GeneSets == "KEGG")
neg_C6 <- neg_df %>% 
  filter(GeneSets == "Oncogenic Signature")

theme_Publication <- function(base_size=7, base_family="helvetica") {
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

plot_KEGG <- ggplot(data = neg_KEGG, 
       aes(x = fct_rev(factor(NAME, level=unique(NAME))),
           y = -NES)) + 
  geom_col(fill = "#FF9933") + 
  coord_flip() + 
  theme_Publication()

png("Y2_KD_RNA-seq_KEGG_barplot_orange.png",
    width = 4500,
    height = 3000,
    res = 600)
plot_KEGG
dev.off()

plot_C6 <- ggplot(data = neg_C6, 
                    aes(x = fct_rev(factor(NAME, level=unique(NAME))),
                        y = -NES)) + 
  geom_col(fill = "#FF00FF") + 
  coord_flip() + 
  theme_Publication()

png("Y2_KD_RNA-seq_C6_barplot_neg.png",
    width = 4500,
    height = 3000,
    res = 600)
plot_C6
dev.off()



#pvalue-------------------------
neg_KEGG_pvalue <- fread(file = "gsea_report_for_shCtrl_1688281337834(KEGG).tsv")
neg_KEGG_pvalue <- neg_KEGG_pvalue %>% 
  arrange(`NOM p-val`)
plot_KEGG_pvalue <- ggplot(data = neg_KEGG_pvalue, 
                    aes(x = fct_rev(factor(NAME, level=unique(NAME))),
                        y = -log10(`NOM p-val`), 
                        fill = -NES)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = 
                         colorRampPalette(
                           brewer.pal(n = 10, 
                                      name = "Greens"))(100)) + 
  theme_Publication()

png("Y2_KD_RNA-seq_KEGG_pvalue_neg.png",
    width = 4500,
    height = 3000,
    res = 600)
plot_KEGG_pvalue
dev.off()


neg_C6_pvalue <- fread(file = "gsea_report_for_shCtrl_1688282183749(C6 oncogenic).tsv")
neg_C6_pvalue <- neg_C6_pvalue %>% 
  arrange(`NOM p-val`)
plot_C6_pvalue <- ggplot(data = neg_C6_pvalue, 
                           aes(x = fct_rev(factor(NAME, level=unique(NAME))),
                               y = -log10(`NOM p-val`), 
                               fill = -NES)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = 
                         colorRampPalette(
                           brewer.pal(n = 10, 
                                      name = "Blues"))(100)) + 
  theme_Publication()

png("Y2_KD_RNA-seq_C6_pvalue_neg.png",
    width = 4500,
    height = 3000,
    res = 600)
plot_C6_pvalue
dev.off()
