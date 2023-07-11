setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/soft_agar_colony_size_boxplot")
getwd()
library(data.table)
df <- fread(file = "for_colony_size_boxplot_NA_by_type.csv")
library(ggplot2)

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = "black"),
            axis.title = element_text(size = rel(1)),
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
            strip.background=element_rect(colour="black",fill="#f0f0f0"),
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

ggplot(df, aes(x=Stage, 
               y=Size, 
               fill=Stage)) + 
  geom_boxplot(notch = FALSE,
               outlier.shape = NA) +
  facet_wrap(~Line) + 
  scale_colour_Publication() + 
  scale_fill_Publication() + 
  theme_Publication() + 
  theme(axis.text = element_text(angle = 90)) + 
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) + 
  labs(y= "Size of colony (μm)")
ggsave("filename.png", dpi=300, height=4, width=5, units="in")

ggplot(df, aes(x=Stage, 
               y=Size, 
               fill=Stage)) + 
  geom_violin(draw_quantiles = TRUE) +
  facet_wrap(~Line) + 
  scale_colour_Publication() + 
  scale_fill_Publication() + 
  theme_Publication() + 
  theme(axis.text = element_text(angle = 90)) + 
  theme(axis.title.x = element_blank(),
        legend.title = element_blank()) + 
  labs(y= "Size of colony (μm)")
ggsave("violin.png", dpi=300, height=4, width=5, units="in")



install.packages("palmerpenguins")
install.packages("ggstatsplot")
library(palmerpenguins)
library(ggstatsplot)
plt <- ggbetweenstats(
  data = df,
  x = Stage,
  y = Size,
  pairwise.comparisons = FALSE
)
plt

ggsave(
  filename = "soft-agar-violinplot-with-ggstatsplot.png",
  plot = plt,
  dpi=300, height=5, width=8, units="in"
)
