setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/Y2_target_Venn")
getwd()
library(VennDiagram)

seq_venn <- 
  read.csv(file = "Y2_seq_overlapping_venn.csv",
           header = TRUE)

library(RColorBrewer)
library(ggplot2)
grid.newpage() 
venn.diagram(
  x = list(seq_venn$RNA.seq, 
           seq_venn$meRIP.seq, 
           seq_venn$eCLIP,
           seq_venn$confident_targets),
  category.names = c("RNA-seq" , 
                     "meRIP.seq" , 
                     "eCLIP",
                     "confident_targets"),
  filename = 'seq_venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 2500 , 
  width = 2500 , 
  resolution = 600,
  compression = "lzw",
  lwd = 0.7,
  col= "black",
  fill = c(alpha("#0080FF",1), 
           alpha('#FF8000',1), 
           alpha('#CC00CC',1),
           alpha('#00FF00',1)),
  cex = 0.7,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0,
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.col = c("#0080FF", '#FF8000', '#CC00CC','#00FF00')
  #,print.mode=c("raw","percent")
)

