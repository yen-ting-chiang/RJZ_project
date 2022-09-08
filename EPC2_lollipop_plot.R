setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/YTHDF2_lollipop_plot/EPC2_lollipop_plot")
getwd()
library(trackViewer)
library(dplyr)

EPC2_cosmic_mutation <- read.csv(file = "EPC2_NM_015630.3.csv")
EPC2_cosmic_mutation <- EPC2_cosmic_mutation %>% 
  filter(Type %in% c("Unknown","Substitution - coding silent") == FALSE)
SNP <- as.numeric(EPC2_cosmic_mutation[,1])
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, 
                                     names=EPC2_cosmic_mutation[,3]))

EPL1 <- "Enhancer of polycomb-like"
E_Pc_C <- "Enhancer of Polycomb C-terminus"

features <- GRanges("chr1", 
                    IRanges(c(6,577), 
                            width = c(142,230),
                            names = c(EPL1,E_Pc_C)))

sample.gr$SNPsideID <- "top"
sample.gr$SNPsideID[50:56] <- "bottom"

sample.gr$legend <- "COSMIC database"
sample.gr$legend[50:56] <- "LFS-OBs-derived mouse tumors"

features$fill <- c("#51C6E6","#FF8833")
features$height <- c(0.05,0.05)
sample.gr$color <- "#4C9900"
sample.gr$color[50:56] <- "#0080FF"
sample.gr$border <- "gray30"
sample.gr$alpha <- 0.6
sample.gr$score <- EPC2_cosmic_mutation[,5]
xaxis <- c(1, 200, 400, 600, 807)
sample.gr$dashline.col <- "gray30"

lolliplot(sample.gr, 
          features, 
          ranges = GRanges("chr1", IRanges(1, 807)),
          xaxis = xaxis,
          yaxis = FALSE,
          type="circle",
          cex = 1,
          jitter="label",
          label_on_feature = FALSE,
          legend="legend")
png("EPC2_1.png", width=1433, height = 765)
dev.off()


