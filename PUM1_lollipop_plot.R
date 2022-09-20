setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/YTHDF2_lollipop_plot/PUM1_lollipop_plot")
getwd()
library(trackViewer)
library(dplyr)

PUM1_cosmic_mutation <- read.csv(file = "PUM1_NM_001020658.1.csv")
PUM1_cosmic_mutation <- PUM1_cosmic_mutation %>% 
  filter(Type %in% c("Unknown","Substitution - coding silent") == FALSE)
SNP <- as.numeric(PUM1_cosmic_mutation[,1])
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, 
                                     names=PUM1_cosmic_mutation[,3]))
PUF <- "PUF: Pumilio-family RNA binding repeat"
features <- GRanges("chr1", 
                    IRanges(c(849,885,921,957,993,1029,1064,1112), 
                            width = c(33,29,30,27,33,32,29,27),
                            names = c(PUF,PUF,PUF,PUF,PUF,PUF,PUF,PUF)))

sample.gr$SNPsideID <- "top"
sample.gr$SNPsideID[44:49] <- "bottom"

sample.gr$legend <- "COSMIC database"
sample.gr$legend[44:49] <- "LFS-OBs-derived mouse tumors"

features$fill <- c("#FF8833","#FF8833","#FF8833","#FF8833",
                   "#FF8833","#FF8833","#FF8833","#FF8833")
features$height <- c(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05)
sample.gr$color <- "#4C9900"
sample.gr$color[44:49] <- "#0080FF"
sample.gr$border <- "gray30"
sample.gr$alpha <- 0.6
sample.gr$score <- PUM1_cosmic_mutation[,5]
xaxis <- c(1, 300, 600, 900, 1186)
sample.gr$dashline.col <- "gray30"

lolliplot(sample.gr, 
          features, 
          ranges = GRanges("chr1", IRanges(1, 1186)),
          xaxis = xaxis,
          yaxis = FALSE,
          type="circle",
          cex = 1,
          jitter="label",
          label_on_feature = FALSE,
          legend="legend")
png("PUM1_1.png", width=1433, height = 765)
dev.off()


