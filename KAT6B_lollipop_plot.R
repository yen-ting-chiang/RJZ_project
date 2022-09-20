setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/YTHDF2_lollipop_plot/KAT6B_lollipop_plot")
getwd()
library(trackViewer)
library(dplyr)

KAT6B_cosmic_mutation <- read.csv(file = "KAT6B_NM_001256468.1.csv")
KAT6B_cosmic_mutation <- KAT6B_cosmic_mutation %>% 
  filter(Type %in% c("Unknown","Substitution - coding silent") == FALSE)
SNP <- as.numeric(KAT6B_cosmic_mutation[,1])
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, 
                                     names=KAT6B_cosmic_mutation[,3]))

Linker_histone <- "linker histone H1 and H5 family "
PHD <- "PHD-finger"
MOZ_SAS <- "MOZ/SAS family"
features <- GRanges("chr1", 
                    IRanges(c(110,272,773), 
                            width = c(57,48,184),
                            names = c(Linker_histone,PHD,MOZ_SAS)))

sample.gr$SNPsideID <- "top"
sample.gr$SNPsideID[88:94] <- "bottom"

sample.gr$legend <- "COSMIC database"
sample.gr$legend[88:94] <- "LFS-OBs-derived mouse tumors"

features$fill <- c("#51C6E6","#FF8833","#0000FF")
features$height <- c(0.05,0.05,0.05)
sample.gr$color <- "#4C9900"
sample.gr$color[88:94] <- "#0080FF"
sample.gr$border <- "gray30"
sample.gr$alpha <- 0.6
sample.gr$score <- KAT6B_cosmic_mutation[,5]
xaxis <- c(1, 500, 1000, 1500, 2073)
sample.gr$dashline.col <- "gray30"

lolliplot(sample.gr, 
          features, 
          ranges = GRanges("chr1", IRanges(1, 2073)),
          xaxis = xaxis,
          yaxis = FALSE,
          type="circle",
          cex = 0.8,
          jitter="label",
          label_on_feature = FALSE,
          legend="legend")
png("KAT6B_0.8.png", width=1433, height = 765)
dev.off()


