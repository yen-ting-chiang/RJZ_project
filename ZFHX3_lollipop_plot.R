setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/YTHDF2_lollipop_plot/ZFHX3_lollipop_plot")
getwd()
library(trackViewer)
library(dplyr)

ZFHX3_cosmic_mutation <- read.csv(file = "ZFHX3_NM_006885.3.csv")
ZFHX3_cosmic_mutation <- ZFHX3_cosmic_mutation %>% 
  filter(Type %in% c("Unknown","Substitution - coding silent") == FALSE)
SNP <- as.numeric(ZFHX3_cosmic_mutation[,1])
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, 
                                     names=ZFHX3_cosmic_mutation[,3]))
features <- GRanges("chr1", 
                    IRanges(c(672, 726, 1596, 2146, 2243, 2642, 2947), 
                            width = c(22, 24, 24, 56, 56, 56, 56),
                            names = c("zf-C2H2", 
                                      "zf-met", 
                                      "zf-C2H2",
                                      "Homeobox",
                                      "Homeobox",
                                      "Homeobox",
                                      "Homeobox")))

sample.gr$SNPsideID <- "top"
sample.gr$SNPsideID[35:41] <- "bottom"

sample.gr$legend <- "COSMIC database"
sample.gr$legend[35:41] <- "LFS-OBs-derived mouse tumors"

features$fill <- c("#FF8833", "#51C6E6")
sample.gr$color <- "#4C9900"
sample.gr$color[35:41] <- "#0080FF"
sample.gr$border <- "gray30"
sample.gr$alpha <- 0.6
features$height <- c(0.07, 0.05)
sample.gr$score <- ZFHX3_cosmic_mutation[,5]
xaxis <- c(1, 100, 200, 300, 400, 500, 579)
yaxis <- c(0, 1, 2, 3)
sample.gr$dashline.col <- "gray30"

lolliplot(sample.gr, 
          features, 
          ranges = GRanges("chr1", IRanges(1, 3703)),
          yaxis = FALSE,
          type="circle",
          cex = 0.3,
          jitter="label",
          ylab="Count",
          yaxis.gp = gpar(fontsize=5, lwd=0),
          xaxis.gp = gpar(fontsize=10, lwd=0),
          label_on_feature = FALSE,
          legend="legend")
png("testplottt2.png", width=1433, height = 765)
dev.off()


