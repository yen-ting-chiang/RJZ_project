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

#enrichr_shARID1B_vs_con_GOBP2021---------------------------------------
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_ZFHX3_RNAseq_enrichr_figure/ARID1B enrichr")
getwd()
enrichr_shARID1B_vs_con_GOBP2021 <- 
  fread(file = "GO_Biological_Process_2021_table.txt")
enrichr_shARID1B_vs_con_GOBP2021_filtered <- 
  enrichr_shARID1B_vs_con_GOBP2021 %>% 
  filter(`Adjusted P-value` < 0.05) %>% 
  arrange(`Combined Score`) %>% 
  mutate(Term = factor(Term, levels = unique(Term)))

plotting_enrichr_shARID1B_vs_con_GOBP2021 <- 
  ggplot(data = enrichr_shARID1B_vs_con_GOBP2021_filtered, 
                       aes(x = Term, 
                           y = `Combined Score`, 
                           fill = `Adjusted P-value`)) + 
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

png("enrichr_shARID1B_vs_con_GOBP2021.png",
    width = 5000,
    height = 3000,
    res = 900)
plotting_enrichr_shARID1B_vs_con_GOBP2021
dev.off()



#enrichr_shARID1B_vs_con_Human_Phenotype_Ontology_table---------------------------------------
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_ZFHX3_RNAseq_enrichr_figure/ARID1B enrichr")
getwd()
enrichr_shARID1B_vs_con_Human_Phenotype_Ontology_table <- 
  fread(file = "Human_Phenotype_Ontology_table.txt")
enrichr_shARID1B_vs_con_Human_Phenotype_Ontology_table_filtered <- 
  enrichr_shARID1B_vs_con_Human_Phenotype_Ontology_table %>% 
  filter(`Adjusted P-value` < 0.05) %>% 
  arrange(`Combined Score`) %>% 
  mutate(Term = factor(Term, levels = unique(Term)))

plotting_enrichr_shARID1B_vs_con_Human_Phenotype_Ontology_table <- 
  ggplot(data = enrichr_shARID1B_vs_con_Human_Phenotype_Ontology_table_filtered, 
         aes(x = Term, 
             y = `Combined Score`, 
             fill = `Adjusted P-value`)) + 
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

png("enrichr_shARID1B_vs_con_Human_Phenotype_Ontology_table.png",
    width = 5000,
    height = 3000,
    res = 600)
plotting_enrichr_shARID1B_vs_con_Human_Phenotype_Ontology_table
dev.off()


#enrichr_shARID1B_vs_con_MGI_Mammalian_Phenotype_Level_4_2021_table---------------------------------------
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_ZFHX3_RNAseq_enrichr_figure/ARID1B enrichr")
getwd()
enrichr_shARID1B_vs_con_MGI_Mammalian_Phenotype_Level_4_2021_table <- 
  fread(file = "MGI_Mammalian_Phenotype_Level_4_2021_table.txt")
enrichr_shARID1B_vs_con_MGI_Mammalian_Phenotype_Level_4_2021_table_filtered <- 
  enrichr_shARID1B_vs_con_MGI_Mammalian_Phenotype_Level_4_2021_table %>% 
  filter(`Adjusted P-value` < 0.05) %>% 
  arrange(`Combined Score`) %>% 
  mutate(Term = factor(Term, levels = unique(Term)))

plotting_enrichr_shARID1B_vs_con_MGI_Mammalian_Phenotype_Level_4_2021_table <- 
  ggplot(data = enrichr_shARID1B_vs_con_MGI_Mammalian_Phenotype_Level_4_2021_table_filtered, 
         aes(x = Term, 
             y = `Combined Score`, 
             fill = `Adjusted P-value`)) + 
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

png("enrichr_shARID1B_vs_con_MGI_Mammalian_Phenotype_Level_4_2021_table.png",
    width = 5000,
    height = 3000,
    res = 600)
plotting_enrichr_shARID1B_vs_con_MGI_Mammalian_Phenotype_Level_4_2021_table
dev.off()

#enrichr_shARID1B_vs_con_GOBP2021---------------------------------------
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_ZFHX3_RNAseq_enrichr_figure/ARID1B enrichr")
getwd()
enrichr_shARID1B_vs_con_MSigDB_Hallmark_2020_table <- 
  fread(file = "MSigDB_Hallmark_2020_table.txt")
enrichr_shARID1B_vs_con_MSigDB_Hallmark_2020_table_filtered <- 
  enrichr_shARID1B_vs_con_MSigDB_Hallmark_2020_table %>% 
  filter(`Adjusted P-value` < 0.05) %>% 
  arrange(`Combined Score`) %>% 
  mutate(Term = factor(Term, levels = unique(Term)))

plotting_enrichr_shARID1B_vs_con_MSigDB_Hallmark_2020_table <- 
  ggplot(data = enrichr_shARID1B_vs_con_MSigDB_Hallmark_2020_table_filtered, 
         aes(x = Term, 
             y = `Combined Score`, 
             fill = `Adjusted P-value`)) + 
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

png("enrichr_shARID1B_vs_con_MSigDB_Hallmark_2020_table.png",
    width = 5000,
    height = 3000,
    res = 600)
plotting_enrichr_shARID1B_vs_con_MSigDB_Hallmark_2020_table
dev.off()




#enrichr_shZFHX3_vs_con_GOBP2021---------------------------------------
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_ZFHX3_RNAseq_enrichr_figure/ZFHX3 enrichr")
getwd()
enrichr_shZFHX3_vs_con_GOBP2021 <- 
  fread(file = "GO_Biological_Process_2021_table.txt")
enrichr_shZFHX3_vs_con_GOBP2021_filtered <- 
  enrichr_shZFHX3_vs_con_GOBP2021 %>% 
  filter(`Adjusted P-value` < 0.05) %>% 
  arrange(`Combined Score`) %>% 
  mutate(Term = factor(Term, levels = unique(Term)))

plotting_enrichr_shZFHX3_vs_con_GOBP2021 <- 
  ggplot(data = enrichr_shZFHX3_vs_con_GOBP2021_filtered, 
         aes(x = Term, 
             y = `Combined Score`, 
             fill = `Adjusted P-value`)) + 
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

png("enrichr_shZFHX3_vs_con_GOBP2021.png",
    width = 5000,
    height = 3000,
    res = 300)
plotting_enrichr_shZFHX3_vs_con_GOBP2021
dev.off()


enrichr_shZFHX3_vs_con_GOBP2021_filtered_top30 <- 
  enrichr_shZFHX3_vs_con_GOBP2021_filtered %>% 
  tail(30)
plotting_enrichr_shZFHX3_vs_con_GOBP2021_top30 <- 
  ggplot(data = enrichr_shZFHX3_vs_con_GOBP2021_filtered_top30, 
         aes(x = Term, 
             y = `Combined Score`, 
             fill = `Adjusted P-value`)) + 
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

png("enrichr_shZFHX3_vs_con_GOBP2021_top30.png",
    width = 5000,
    height = 3000,
    res = 600)
plotting_enrichr_shZFHX3_vs_con_GOBP2021_top30
dev.off()

#enrichr_shZFHX3_vs_con_BioPlanet_2019_table---------------------------------------
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_ZFHX3_RNAseq_enrichr_figure/ZFHX3 enrichr")
getwd()
enrichr_shZFHX3_vs_con_BioPlanet_2019_table <- 
  fread(file = "BioPlanet_2019_table.txt")
enrichr_shZFHX3_vs_con_BioPlanet_2019_table_filtered <- 
  enrichr_shZFHX3_vs_con_BioPlanet_2019_table %>% 
  filter(`Adjusted P-value` < 0.05) %>% 
  arrange(`Combined Score`) %>% 
  mutate(Term = factor(Term, levels = unique(Term)))

plotting_enrichr_shZFHX3_vs_con_BioPlanet_2019_table <- 
  ggplot(data = enrichr_shZFHX3_vs_con_BioPlanet_2019_table_filtered, 
         aes(x = Term, 
             y = `Combined Score`, 
             fill = `Adjusted P-value`)) + 
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

png("enrichr_shZFHX3_vs_con_BioPlanet_2019_table.png",
    width = 5000,
    height = 3000,
    res = 600)
plotting_enrichr_shZFHX3_vs_con_BioPlanet_2019_table
dev.off()


#enrichr_shZFHX3_vs_con_MSigDB_Hallmark_2020_table---------------------------------------
setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_KD_ARID1B_ZFHX3_RNAseq_enrichr_figure/ZFHX3 enrichr")
getwd()
enrichr_shZFHX3_vs_con_MSigDB_Hallmark_2020_table <- 
  fread(file = "MSigDB_Hallmark_2020_table.txt")
enrichr_shZFHX3_vs_con_MSigDB_Hallmark_2020_table_filtered <- 
  enrichr_shZFHX3_vs_con_MSigDB_Hallmark_2020_table %>% 
  filter(`Adjusted P-value` < 0.05) %>% 
  arrange(`Combined Score`) %>% 
  mutate(Term = factor(Term, levels = unique(Term)))

plotting_enrichr_shZFHX3_vs_con_MSigDB_Hallmark_2020_table <- 
  ggplot(data = enrichr_shZFHX3_vs_con_MSigDB_Hallmark_2020_table_filtered, 
         aes(x = Term, 
             y = `Combined Score`, 
             fill = `Adjusted P-value`)) + 
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

png("enrichr_shZFHX3_vs_con_MSigDB_Hallmark_2020_table.png",
    width = 5000,
    height = 3000,
    res = 1200)
plotting_enrichr_shZFHX3_vs_con_MSigDB_Hallmark_2020_table
dev.off()

