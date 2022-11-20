library(dplyr)
library(data.table)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
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

setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_ZFHX3_RNAseq_DEG_figure")
getwd()
#DEG_ARID1B------------------------------------------
DEG_ARID1B <- 
  fread(file = "DEG_ARID1B.tabular")
DEG_ARID1B_arranged <- 
  DEG_ARID1B %>% 
  arrange(`log2(FC)`) %>% 
  mutate(`Gene name` = 
           factor(`Gene name`, 
                  levels = unique(`Gene name`)))

plotting_DEG_ARID1B <- 
  ggplot(data = DEG_ARID1B_arranged, 
         aes(x = `Gene name`, 
             y = `log2(FC)`, 
             fill = `P-adj`)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = 
                         colorRampPalette(
                           rev(brewer.pal(n = 9, 
                                      name = "YlGnBu")))(100)) + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))

png("DEG_ARID1B.png",
    width = 5000,
    height = 3000,
    res = 900)
plotting_DEG_ARID1B
dev.off()


#DEG_ZFHX3------------------------------------------
DEG_ZFHX3 <- 
  fread(file = "DEG_ZFHX3.tabular")
DEG_ZFHX3_arranged <- 
  DEG_ZFHX3 %>% 
  arrange(`log2(FC)`) %>% 
  mutate(`Gene name` = 
           factor(`Gene name`, 
                  levels = unique(`Gene name`)))

plotting_DEG_ZFHX3 <- 
  ggplot(data = DEG_ZFHX3_arranged, 
         aes(x = `Gene name`, 
             y = `log2(FC)`, 
             fill = `P-adj`)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = 
                         colorRampPalette(
                           rev(brewer.pal(n = 9, 
                                          name = "YlGnBu")))(100)) + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))

png("DEG_ZFHX3.png",
    width = 5000,
    height = 3000,
    res = 300)
plotting_DEG_ZFHX3
dev.off()


DEG_ZFHX3_top30 <- 
  DEG_ZFHX3 %>% 
  mutate(`log2(FC)_abs` = abs(`log2(FC)`)) %>% 
  arrange(`log2(FC)_abs`) %>% 
  tail(30) %>% 
  arrange(`log2(FC)`) %>% 
  mutate(`Gene name` = 
           factor(`Gene name`, 
                  levels = unique(`Gene name`)))

plotting_DEG_ZFHX3_top30 <- 
  ggplot(data = DEG_ZFHX3_top30, 
         aes(x = `Gene name`, 
             y = `log2(FC)`, 
             fill = `P-adj`)) + 
  geom_col() + 
  coord_flip() + 
  scale_fill_gradientn(colours = 
                         colorRampPalette(
                           rev(brewer.pal(n = 9, 
                                          name = "YlGnBu")))(100)) + 
  theme_Publication() + 
  theme(panel.border = 
          element_rect(colour = "black", 
                       fill = NA, 
                       size = 1))

png("DEG_ZFHX3_top30.png",
    width = 5000,
    height = 3000,
    res = 600)
plotting_DEG_ZFHX3_top30
dev.off()
