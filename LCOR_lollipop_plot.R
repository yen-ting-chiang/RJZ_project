setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/YTHDF2_lollipop_plot/LCOR_lollipop_plot")
getwd()
library(trackViewer)
library(dplyr)

LCOR_cosmic_mutation <- read.csv(file = "LCOR_NM_001170766.1.csv")
LCOR_cosmic_mutation <- LCOR_cosmic_mutation %>% 
  filter(Type %in% c("Unknown","Substitution - coding silent") == FALSE)
SNP <- as.numeric(LCOR_cosmic_mutation[,1])
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, 
                                     names=LCOR_cosmic_mutation[,3]))

HTH_psq <- "helix-turn-helix, Psq domain (351 - 394)"
features <- GRanges("chr1", 
                    IRanges(c(351), 
                            width = c(43),
                            names = c(HTH_psq)))

sample.gr$SNPsideID <- "top"
sample.gr$SNPsideID[19:25] <- "bottom"

sample.gr$legend <- "COSMIC database"
sample.gr$legend[19:25] <- "LFS-OBs-derived mouse tumors"

features$fill <- c("#FF8833")
features$height <- c(0.05)
sample.gr$color <- "#4C9900"
sample.gr$color[19:25] <- "#0080FF"
sample.gr$border <- "gray30"
sample.gr$alpha <- 0.6
sample.gr$score <- LCOR_cosmic_mutation[,5]
xaxis <- c(1, 100, 200, 300, 433)
sample.gr$dashline.col <- "gray30"

lolliplot(sample.gr, 
          features, 
          ranges = GRanges("chr1", IRanges(1, 433)),
          xaxis = xaxis,
          yaxis = FALSE,
          type="circle",
          cex = 1,
          jitter="label",
          label_on_feature = FALSE,
          legend="legend")
png("LCOR_1.png", width=1433, height = 765)
dev.off()


