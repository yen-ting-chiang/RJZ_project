setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/Y2_KD_RNA-seq_m6A_gene_expression_compare")
getwd()
library(data.table)
library(dplyr)
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
fpkm <- fread(file = "fpkm.csv")
m6A_genes <- fread(file = "m6A_genes.csv")

df_shCtrl_3 <- fpkm[,1:2]
df_shCtrl_3 <- df_shCtrl_3 %>% 
  filter(`shCtrl-1`>10)
df_shCtrl_3_m6A <- df_shCtrl_3 %>% 
  filter(V1 %in% m6A_genes$m6A_genes == TRUE) %>% 
  mutate(m6A_status = "m6A_genes")
df_shCtrl_3_non_m6A <- df_shCtrl_3 %>% 
  filter(V1 %in% m6A_genes$m6A_genes == FALSE) %>% 
  mutate(m6A_status = "non_m6A_genes")
df_shCtrl_3 <- rbind(df_shCtrl_3_m6A,df_shCtrl_3_non_m6A) %>% 
  na.omit()

library(ggplot2)
plot_df_shCtrl_3 <- ggplot(data = df_shCtrl_3,
                           aes(x = m6A_status,
                               y = log10(`shCtrl-1`),
                               fill = m6A_status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+ 
  theme_Publication()+
  scale_fill_Publication() + 
  scale_fill_brewer(palette = "Set1")

png("violin_chart_m6A_Y2_strict.png",
    width = 4500,
    height = 3000,
    res = 600)
plot_violin_chart_m6A_Y2_strict
dev.off()

m6A_Y2_strict_binding_df %>% 
  group_by(m6A_status) %>% 
  summarise(mean(`shCtrl-1`))

t.test(`shCtrl-1` ~ `m6A_status`, data = df_shCtrl_2)
wilcox.test(`shCtrl-1` ~ `m6A_status`, data = df_shCtrl_2)

#shCtrl-2------------------------------
df_shCtrl_2 <- fpkm[,c(1,3)]
df_shCtrl_2 <- df_shCtrl_2 %>% 
  filter(`shCtrl-2`>10)
df_shCtrl_2_m6A <- df_shCtrl_2 %>% 
  filter(V1 %in% m6A_genes$m6A_genes == TRUE) %>% 
  mutate(m6A_status = "m6A_genes")
df_shCtrl_2_non_m6A <- df_shCtrl_2 %>% 
  filter(V1 %in% m6A_genes$m6A_genes == FALSE) %>% 
  mutate(m6A_status = "non_m6A_genes")
df_shCtrl_2 <- rbind(df_shCtrl_2_m6A,df_shCtrl_2_non_m6A) %>% 
  na.omit()

plot_df_shCtrl_2 <- ggplot(data = df_shCtrl_2,
                           aes(x = m6A_status,
                               y = log10(`shCtrl-2`),
                               fill = m6A_status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+ 
  theme_Publication()+
  scale_fill_Publication() + 
  scale_fill_brewer(palette = "Set1")

png("violin_chart_m6A_Y2_strict.png",
    width = 4500,
    height = 3000,
    res = 600)
plot_violin_chart_m6A_Y2_strict
dev.off()

m6A_Y2_strict_binding_df %>% 
  group_by(m6A_status) %>% 
  summarise(mean(`shCtrl-2`))

t.test(`shCtrl-2` ~ `m6A_status`, data = df_shCtrl_2)
wilcox.test(`shCtrl-2` ~ `m6A_status`, data = df_shCtrl_2)




#shCtrl-3--------------------------------------
df_shCtrl_3 <- fpkm[,c(1,4)]
df_shCtrl_3 <- df_shCtrl_3 %>% 
  filter(`shCtrl-3`>10)
df_shCtrl_3_m6A <- df_shCtrl_3 %>% 
  filter(V1 %in% m6A_genes$m6A_genes == TRUE) %>% 
  mutate(m6A_status = "m6A_genes")
df_shCtrl_3_non_m6A <- df_shCtrl_3 %>% 
  filter(V1 %in% m6A_genes$m6A_genes == FALSE) %>% 
  mutate(m6A_status = "non_m6A_genes")
df_shCtrl_3 <- rbind(df_shCtrl_3_m6A,df_shCtrl_3_non_m6A) %>% 
  na.omit()
# plot_df_shCtrl_3 <- 
  ggplot(data = df_shCtrl_3,
                           aes(x = m6A_status,
                               y = log10(`shCtrl-3`),
                               fill = m6A_status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+ 
  theme_Publication()+
  scale_fill_Publication() + 
  scale_fill_brewer(palette = "Set1")
# png("violin_chart_m6A_Y2_strict.png",
#     width = 4500,
#     height = 3000,
#     res = 600)
# plot_violin_chart_m6A_Y2_strict
# dev.off()
t.test(`shCtrl-3` ~ `m6A_status`, data = df_shCtrl_3)
wilcox.test(`shCtrl-3` ~ `m6A_status`, data = df_shCtrl_3)

#shYTHDF2-3---------------------------------------
df_shYTHDF2_3 <- fpkm[,c(1,7)]
df_shYTHDF2_3 <- df_shYTHDF2_3 %>% 
  filter(`shYTHDF2-3`>10)
df_shYTHDF2_3_m6A <- df_shYTHDF2_3 %>% 
  filter(V1 %in% m6A_genes$m6A_genes == TRUE) %>% 
  mutate(m6A_status = "m6A_genes")
df_shYTHDF2_3_non_m6A <- df_shYTHDF2_3 %>% 
  filter(V1 %in% m6A_genes$m6A_genes == FALSE) %>% 
  mutate(m6A_status = "non_m6A_genes")
df_shYTHDF2_3 <- rbind(df_shYTHDF2_3_m6A,df_shYTHDF2_3_non_m6A) %>% 
  na.omit()
# plot_df_shYTHDF2_3 <- 
ggplot(data = df_shYTHDF2_3,
       aes(x = m6A_status,
           y = log10(`shYTHDF2-3`),
           fill = m6A_status)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+ 
  theme_Publication()+
  scale_fill_Publication() + 
  scale_fill_brewer(palette = "Set1")
# png("violin_chart_m6A_Y2_strict.png",
#     width = 4500,
#     height = 3000,
#     res = 600)
# plot_violin_chart_m6A_Y2_strict
# dev.off()
t.test(`shYTHDF2-3` ~ `m6A_status`, data = df_shYTHDF2_3)
wilcox.test(`shYTHDF2-3` ~ `m6A_status`, data = df_shYTHDF2_3)
