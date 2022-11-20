library(dplyr)
library(data.table)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

#GSEA_shARID1B_vs_con_C6---------------------------------------
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_RNAseq_figure/GSEA_shARID1B_vs_con_C6")
getwd()
  ## GSEA_shARID1B_vs_con_C6_v2022_pos----------------------------
GSEA_shARID1B_vs_con_C6_v2022_pos <- 
  fread(file = "gsea_report_for_shARID1B_1668714521799.tsv")
GSEA_shARID1B_vs_con_C6_v2022_pos_filtered <- 
  GSEA_shARID1B_vs_con_C6_v2022_pos %>% 
  filter(`NOM p-val` < 0.05) %>% 
  filter(`FDR q-val` <0.25) %>% 
  arrange(NES) %>% 
  mutate(NAME = factor(NAME, levels = unique(NAME)))

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

ggplot(data = GSEA_shARID1B_vs_con_C6_v2022_pos_filtered, 
       aes(x = NAME, y = NES, fill = `NOM p-val`)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = 
                         colorRampPalette(
                           brewer.pal(n = 5, 
                                      name = "RdBu"))(100)) + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))

  ## GSEA_shARID1B_vs_con_C6_v2022_neg----------------------------
GSEA_shARID1B_vs_con_C6_v2022_neg <- 
  fread(file = "gsea_report_for_control_1668714521799.tsv")
GSEA_shARID1B_vs_con_C6_v2022_neg_filtered <- 
  GSEA_shARID1B_vs_con_C6_v2022_neg %>% 
  filter(`NOM p-val` < 0.05) %>% 
  filter(`FDR q-val` <0.25) %>% 
  arrange(NES) %>% 
  mutate(NAME = factor(NAME, levels = unique(NAME)))

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

ggplot(data = GSEA_shARID1B_vs_con_C6_v2022_neg_filtered, 
       aes(x = NAME, y = NES, fill = `NOM p-val`)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = 
                         colorRampPalette(
                           brewer.pal(n = 5, 
                                      name = "RdBu"))(100)) + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))

  ##combine_pos_neg--------------------------------
GSEA_shARID1B_vs_con_C6_v2022 <- 
  rbind(GSEA_shARID1B_vs_con_C6_v2022_pos,
        GSEA_shARID1B_vs_con_C6_v2022_neg)
GSEA_shARID1B_vs_con_C6_v2022_filtered <- 
  GSEA_shARID1B_vs_con_C6_v2022 %>% 
  filter(`NOM p-val` < 0.05) %>% 
  filter(`FDR q-val` <0.25) %>% 
  arrange(NES) %>% 
  mutate(NAME = factor(NAME, levels = unique(NAME)))

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

plotting <- ggplot(data = GSEA_shARID1B_vs_con_C6_v2022_filtered, 
       aes(x = NAME, y = NES, fill = `NOM p-val`)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = 
                         colorRampPalette(
                           brewer.pal(n = 5, 
                                      name = "RdBu"))(100)) + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))


png("GSEA_shARID1B_vs_con_C6_v2022.png",
    width = 3000,
    height = 3000,
    res = 600)
plotting
dev.off()


#GSEA_shARID1B_vs_con_H---------------------------------------
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_RNAseq_figure/GSEA_shARID1B_vs_con_H")
getwd()
##combine_pos_neg--------------------------------
GSEA_shARID1B_vs_con_H_v2022_pos <- 
  fread(file = "gsea_report_for_shARID1B_1668714383317.tsv")
GSEA_shARID1B_vs_con_H_v2022_neg <- 
  fread(file = "gsea_report_for_control_1668714383317.tsv")
GSEA_shARID1B_vs_con_H_v2022 <- 
  rbind(GSEA_shARID1B_vs_con_H_v2022_pos,
        GSEA_shARID1B_vs_con_H_v2022_neg)
GSEA_shARID1B_vs_con_H_v2022_filtered <- 
  GSEA_shARID1B_vs_con_H_v2022 %>% 
  filter(`NOM p-val` < 0.05) %>% 
  filter(`FDR q-val` <0.25) %>% 
  arrange(NES) %>% 
  mutate(NAME = factor(NAME, levels = unique(NAME)))

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

plotting_H_ <- ggplot(data = GSEA_shARID1B_vs_con_H_v2022_filtered, 
                   aes(x = NAME, y = NES, fill = `NOM p-val`)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = 
                         colorRampPalette(
                           brewer.pal(n = 5, 
                                      name = "RdBu"))(100)) + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))


png("GSEA_shARID1B_vs_con_H_v2022.png",
    width = 3000,
    height = 3000,
    res = 600)
plotting_H_
dev.off()


#GSEA_shZFHX3_vs_con_H---------------------------------------
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_RNAseq_figure/GSEA_shZFHX3_vs_con_H")
getwd()
##combine_pos_neg--------------------------------
GSEA_shZFHX3_vs_con_H_v2022_pos <- 
  fread(file = "gsea_report_for_shZFHX3_1668716121393.tsv")
GSEA_shZFHX3_vs_con_H_v2022_neg <- 
  fread(file = "gsea_report_for_control_1668716121393.tsv")
GSEA_shZFHX3_vs_con_H_v2022 <- 
  GSEA_shZFHX3_vs_con_H_v2022_pos
GSEA_shZFHX3_vs_con_H_v2022_filtered <- 
  GSEA_shZFHX3_vs_con_H_v2022 %>% 
  filter(`NOM p-val` < 0.05) %>% 
  filter(`FDR q-val` <0.25) %>% 
  arrange(NES) %>% 
  mutate(NAME = factor(NAME, levels = unique(NAME)))

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

plotting_H_ <- ggplot(data = GSEA_shZFHX3_vs_con_H_v2022_filtered, 
                      aes(x = NAME, y = NES, fill = `NOM p-val`)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = 
                         colorRampPalette(
                           brewer.pal(n = 5, 
                                      name = "RdBu"))(100)) + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))


png("GSEA_shZFHX3_vs_con_H_v2022.png",
    width = 3000,
    height = 3000,
    res = 600)
plotting_H_
dev.off()


#GSEA_shZFHX3_vs_con_C6---------------------------------------
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_RNAseq_figure/GSEA_shZFHX3_vs_con_C6")
getwd()
##combine_pos_neg--------------------------------
GSEA_shZFHX3_vs_con_C6_v2022_pos <- 
  fread(file = "gsea_report_for_shZFHX3_1668716797937.tsv")
GSEA_shZFHX3_vs_con_C6_v2022_neg <- 
  fread(file = "gsea_report_for_control_1668716797937.tsv")
GSEA_shZFHX3_vs_con_C6_v2022 <- 
  GSEA_shZFHX3_vs_con_C6_v2022_pos
GSEA_shZFHX3_vs_con_C6_v2022_filtered <- 
  GSEA_shZFHX3_vs_con_C6_v2022 %>% 
  filter(`NOM p-val` < 0.05) %>% 
  filter(`FDR q-val` <0.25) %>% 
  arrange(NES) %>% 
  mutate(NAME = factor(NAME, levels = unique(NAME)))

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

plotting_C6_ <- ggplot(data = GSEA_shZFHX3_vs_con_C6_v2022_filtered, 
                      aes(x = NAME, y = NES, fill = `NOM p-val`)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = 
                         colorRampPalette(
                           brewer.pal(n = 5, 
                                      name = "RdBu"))(100)) + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))


png("GSEA_shZFHX3_vs_con_C6_v2022.png",
    width = 3000,
    height = 3000,
    res = 600)
plotting_C6_
dev.off()

