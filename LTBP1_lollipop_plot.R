setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/YTHDF2_lollipop_plot/LTBP1_lollipop_plot")
getwd()
library(trackViewer)
library(dplyr)

LTBP1_cosmic_mutation <- read.csv(file = "LTBP1_NM_206943.2.csv")
LTBP1_cosmic_mutation <- LTBP1_cosmic_mutation %>% 
  filter(Type %in% c("Unknown","Substitution - coding silent") == FALSE)
SNP <- as.numeric(LTBP1_cosmic_mutation[,1])
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, 
                                     names=LTBP1_cosmic_mutation[,3]))

TB <- "TB domain"
EGF_CA <- "Calcium-binding EGF domain"
cEGF <- "Complement Clr-like EGF-like"
features <- GRanges("chr1", 
                    IRanges(c(566,626,687,915,957,998,1059,1120,1161,
                              1202,1244,1286,1357,1424,1467,1534,1662), 
                            width = c(44,39,42,40,39,38,33,39,39,40,40,
                                      41,44,41,39,43,43),
                            names = c(TB,EGF_CA,TB,EGF_CA,EGF_CA,EGF_CA,cEGF,EGF_CA,
                                      EGF_CA,EGF_CA,EGF_CA,EGF_CA,TB,EGF_CA,
                                      EGF_CA,TB,EGF_CA)))

sample.gr$SNPsideID <- "top"
sample.gr$SNPsideID[95:102] <- "bottom"

sample.gr$legend <- "COSMIC database"
sample.gr$legend[95:102] <- "LFS-OBs-derived mouse tumors"

features$fill <- c("#51C6E6","#FF8833","#51C6E6","#FF8833",
                   "#FF8833","#FF8833","#0000FF","#FF8833",
                   "#FF8833","#FF8833","#FF8833","#FF8833",
                   "#51C6E6","#FF8833","#FF8833","#51C6E6",
                   "#FF8833")
features$height <- c(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
                     0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,
                     0.05)
sample.gr$color <- "#4C9900"
sample.gr$color[95:102] <- "#0080FF"
sample.gr$border <- "gray30"
sample.gr$alpha <- 0.6
sample.gr$score <- LTBP1_cosmic_mutation[,5]
xaxis <- c(1, 500, 1000, 1500, 1721)
sample.gr$dashline.col <- "gray30"

lolliplot(sample.gr, 
          features, 
          ranges = GRanges("chr1", IRanges(1, 1721)),
          xaxis = xaxis,
          yaxis = FALSE,
          type="circle",
          cex = 0.8,
          jitter="label",
          label_on_feature = FALSE,
          legend="legend")
png("LTBP1_0.8.png", width=1433, height = 765)
dev.off()


