setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/OB_tumor_Venn")
getwd()
library(VennDiagram)
T2_T4_T5_venn <- 
  read.csv(file = "T2_T4_T5_venn.csv",
           header = TRUE)
grid.newpage() 
library(RColorBrewer)
myCol <- brewer.pal(3, "Dark2")
library(ggplot2)

# Chart
venn.diagram(
  x = list(T2_T4_T5_venn$T2, 
           T2_T4_T5_venn$T4, 
           T2_T4_T5_venn$T5),
  category.names = c("Tumor #2" , 
                     "Tumor #4" , 
                     "Tumor #5"),
  filename = 'T2_T4_T5_venn_diagramm_2.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 0.8,
  col="black",
  fill = c(alpha("#0080FF",1), 
           alpha('#FF8000',1), 
           alpha('#CC00CC',1)),
  cex = 0,
  fontfamily = "sans",
  cat.cex = 0,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#0080FF", '#FF8000', '#CC00CC'),
  rotation = 1
)

# Chart
venn.diagram(
  x = list(T2_T4_T5_venn$OB, 
           T2_T4_T5_venn$T2),
  category.names = c("OB D7" , 
                     "Tumor #2"),
  filename = 'OB_T2_venn_diagramm_2.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 0.8,
  col="black",
  fill = c(alpha("#66CC00",1), 
           alpha("#0080FF",1)),
  cex = 0,
  fontfamily = "sans",
  cat.cex = 0,
  cat.default.pos = "outer",
  cat.pos = c(27, -27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("#66CC00", '#0080FF'),
  inverted = TRUE,
  scaled = FALSE
)

# Chart
venn.diagram(
  x = list(T2_T4_T5_venn$OB, 
           T2_T4_T5_venn$T4),
  category.names = c("OB D7" , 
                     "Tumor #4"),
  filename = 'OB_T4_venn_diagramm_2.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 0.8,
  col="black",
  fill = c(alpha("#66CC00",1), 
           alpha("#FF8000",1)),
  cex = 0,
  fontfamily = "sans",
  cat.cex = 0,
  cat.default.pos = "outer",
  cat.pos = c(27, -27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("#66CC00", '#FF8000'),
  inverted = TRUE,
  scaled = FALSE
)

# Chart
venn.diagram(
  x = list(T2_T4_T5_venn$OB, 
           T2_T4_T5_venn$T5),
  category.names = c("Tumor #5" , 
                     "OB D7"),
  filename = 'OB_T5_venn_diagramm_2.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 0.8,
  col="black",
  fill = c(alpha("#66CC00",1), 
           alpha("#CC00CC",1)),
  cex = 0,
  fontfamily = "sans",
  cat.cex = 0,
  cat.default.pos = "outer",
  cat.pos = c(27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = c("#66CC00", '#CC00CC'),
  inverted = TRUE,
  scaled = FALSE
)


OB_mutation_related_paper <- 
  read.csv(file = "OB_mutation_related_paper.csv",
           header = TRUE)
venn.diagram(
  x = list(OB_mutation_related_paper$LFS.OBs.derived.tumor, 
           OB_mutation_related_paper$Sleeping.Beauty.genetic.screening..Moriarity.et.al..2015., 
           OB_mutation_related_paper$Validated.point.mutations.in.osteosarcoma.samples..Chen.et.al.2014.),
  category.names = c("" , 
                     "" , 
                     ""),
  filename = 'OB_mutation_related_paper_venn_diagramm_2.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 0.8,
  col="black",
  fill = c(alpha("#FF0000",1), 
           alpha('#FF8000',1), 
           alpha('#006666',1)),
  cex = 0,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#CC0000", '#006666', '#660066'),
  rotation = 1
)


