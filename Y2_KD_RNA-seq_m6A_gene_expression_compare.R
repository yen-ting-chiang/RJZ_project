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

df <- fread(file = "Annotated_DESeq2_results.csv")
m6A_genes <- fread(file = "m6A_genes.csv")
df_m6A <- df %>% 
  filter(`Gene name` %in% m6A_genes$m6A_genes == TRUE) %>% 
  mutate(m6A_status = "m6A_genes")
df_non_m6A <- df %>% 
  filter(`Gene name` %in% m6A_genes$m6A_genes == FALSE) %>% 
  mutate(m6A_status = "non_m6A_genes")
new_df <- rbind(df_m6A,df_non_m6A) %>% 
  select(`log2(FC)`, m6A_status) %>% 
  na.omit()
library(ggplot2)
plot_violin_chart <- ggplot(data = new_df,
                            aes(x = m6A_status,
                                y = `log2(FC)`)) +
  geom_violin()
new_df %>% 
  group_by(m6A_status) %>% 
  summarise(mean(`log2(FC)`))

#compare strict Y2 targets and non-targets in m6A genes---------------------
YTHDF2_targets <- fread(file = "YTHDF2_targets.csv")
df_m6A_Y2_strict <- df_m6A %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$eCLIP_and_Hela_cells_and_m6A == TRUE) %>% 
  mutate(Y2_status = "YTHDF2_targets_genes")
df_m6A_non_Y2_strict <- df_m6A %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$eCLIP_and_Hela_cells_and_m6A == FALSE) %>% 
  mutate(Y2_status = "non_YTHDF2_targets_genes")
m6A_Y2_strict_binding_df <- 
  rbind(df_m6A_Y2_strict,df_m6A_non_Y2_strict) %>% 
  select(`log2(FC)`, Y2_status) %>% 
  na.omit()
library(ggplot2)
plot_violin_chart_m6A_Y2_strict <- ggplot(data = m6A_Y2_strict_binding_df,
                            aes(x = Y2_status,
                                y = `log2(FC)`,
                                fill = Y2_status)) +
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
  group_by(Y2_status) %>% 
  summarise(mean(`log2(FC)`))

t.test(`log2(FC)` ~ `Y2_status`, data = m6A_Y2_strict_binding_df)
wilcox.test(`log2(FC)` ~ `Y2_status`, data = m6A_Y2_strict_binding_df)



#compare loose Y2 targets and non-targets in m6A genes---------------------
YTHDF2_targets <- fread(file = "YTHDF2_targets.csv")
df_m6A_Y2_loose <- df_m6A %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$eCLIP_or_Hela_cells_and_m6A == TRUE) %>% 
  mutate(Y2_status = "YTHDF2_targets_genes")
df_m6A_non_Y2_loose <- df_m6A %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$eCLIP_or_Hela_cells_and_m6A == FALSE) %>% 
  mutate(Y2_status = "non_YTHDF2_targets_genes")
m6A_Y2_loose_binding_df <- 
  rbind(df_m6A_Y2_loose,df_m6A_non_Y2_loose) %>% 
  select(`log2(FC)`, Y2_status) %>% 
  na.omit()

plot_violin_chart_m6A_Y2_loose <- ggplot(data = m6A_Y2_loose_binding_df,
                                          aes(x = Y2_status,
                                              y = `log2(FC)`)) +
  geom_violin()

m6A_Y2_loose_binding_df %>% 
  group_by(Y2_status) %>% 
  summarise(mean(`log2(FC)`))
t.test(`log2(FC)` ~ `Y2_status`, data = m6A_Y2_loose_binding_df)


#compare Hela Y2 targets and non-targets in m6A genes---------------------
YTHDF2_targets <- fread(file = "YTHDF2_targets.csv")
df_m6A_Y2_Hela <- df_m6A %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$Hela_cells == TRUE) %>% 
  mutate(Y2_status = "YTHDF2_targets_genes")
df_m6A_non_Y2_Hela <- df_m6A %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$Hela_cells == FALSE) %>% 
  mutate(Y2_status = "non_YTHDF2_targets_genes")
m6A_Y2_Hela_binding_df <- 
  rbind(df_m6A_Y2_Hela,df_m6A_non_Y2_Hela) %>% 
  select(`log2(FC)`, Y2_status) %>% 
  na.omit()

plot_violin_chart_m6A_Y2_Hela <- ggplot(data = m6A_Y2_Hela_binding_df,
                                         aes(x = Y2_status,
                                             y = `log2(FC)`)) +
  geom_violin()

m6A_Y2_Hela_binding_df %>% 
  group_by(Y2_status) %>% 
  summarise(mean(`log2(FC)`))
t.test(`log2(FC)` ~ `Y2_status`, data = m6A_Y2_Hela_binding_df)
wilcox.test(`log2(FC)` ~ `Y2_status`, data = m6A_Y2_Hela_binding_df)


#compare astrocyte_Y2 Y2 targets and non-targets in m6A genes---------------------
YTHDF2_targets <- fread(file = "YTHDF2_targets.csv")
df_m6A_Y2_astrocyte_Y2 <- df_m6A %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$astrocyte_Y2_targets == TRUE) %>% 
  mutate(Y2_status = "YTHDF2_targets_genes")
df_m6A_non_Y2_astrocyte_Y2 <- df_m6A %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$astrocyte_Y2_targets == FALSE) %>% 
  mutate(Y2_status = "non_YTHDF2_targets_genes")
m6A_Y2_astrocyte_Y2_binding_df <- 
  rbind(df_m6A_Y2_astrocyte_Y2,df_m6A_non_Y2_astrocyte_Y2) %>% 
  select(`log2(FC)`, Y2_status) %>% 
  na.omit()

plot_violin_chart_m6A_Y2_astrocyte_Y2 <- ggplot(data = m6A_Y2_astrocyte_Y2_binding_df,
                                        aes(x = Y2_status,
                                            y = `log2(FC)`,
                                            fill = Y2_status),
                                        ) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+ 
  theme_Publication()+
  scale_fill_Publication()+ 
  scale_fill_brewer(palette = "Set3")

png("violin_chart_m6A_Y2_astrocyte_Y2.png",
    width = 4500,
    height = 3000,
    res = 600)
plot_violin_chart_m6A_Y2_astrocyte_Y2
dev.off()


m6A_Y2_astrocyte_Y2_binding_df %>% 
  group_by(Y2_status) %>% 
  summarise(mean(`log2(FC)`))
t.test(`log2(FC)` ~ `Y2_status`, data = m6A_Y2_astrocyte_Y2_binding_df)
wilcox.test(`log2(FC)` ~ `Y2_status`, data = m6A_Y2_astrocyte_Y2_binding_df)














#Y2 targets expression comparison---------------------
YTHDF2_targets <- fread(file = "YTHDF2_targets.csv")
df_Y2_strict <- df %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$eCLIP_and_Hela_cells_and_m6A == TRUE) %>% 
  mutate(Y2_status = "YTHDF2_targets_genes")
df_non_Y2_strict <- df %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$eCLIP_and_Hela_cells_and_m6A == FALSE) %>% 
  mutate(Y2_status = "non_YTHDF2_targets_genes")
new_df_strict <- rbind(df_Y2_strict,df_non_Y2_strict) %>% 
  select(`log2(FC)`, Y2_status) %>% 
  na.omit()
library(ggplot2)
plot_violin_chart <- ggplot(data = new_df_strict,
                            aes(x = Y2_status,
                                y = `log2(FC)`)) +
  geom_violin()
new_df_strict %>% 
  group_by(Y2_status) %>% 
  summarise(mean(`log2(FC)`))

#loose Y2 targets expression comparison---------------------
df_Y2_loose <- df %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$eCLIP_or_Hela_cells_and_m6A == TRUE) %>% 
  mutate(Y2_status = "YTHDF2_targets_genes")
df_non_Y2_loose <- df %>% 
  filter(`Gene name` %in% 
           YTHDF2_targets$eCLIP_or_Hela_cells_and_m6A == FALSE) %>% 
  mutate(Y2_status = "non_YTHDF2_targets_genes")
new_df_loose <- rbind(df_Y2_loose,df_non_Y2_loose) %>% 
  select(`log2(FC)`, Y2_status) %>% 
  na.omit()
plot_violin_chart <- ggplot(data = new_df_strict,
                            aes(x = Y2_status,
                                y = `log2(FC)`)) +
  geom_violin()
new_df_loose %>% 
  group_by(Y2_status) %>% 
  summarise(mean(`log2(FC)`))
