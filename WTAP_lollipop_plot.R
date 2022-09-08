setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/YTHDF2_lollipop_plot/WTAP_lollipop_plot")
getwd()
library(trackViewer)
library(dplyr)

WTAP_cosmic_mutation <- read.csv(file = "WTAP_NM_001270531.1.csv")
WTAP_cosmic_mutation <- WTAP_cosmic_mutation %>% 
  filter(Type %in% c("Unknown","Substitution - coding silent") == FALSE)
SNP <- as.numeric(WTAP_cosmic_mutation[,1])
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, 
                                     names=WTAP_cosmic_mutation[,3]))
features <- GRanges("chr1", 
                    IRanges(c(1), 
                            width = c(245),
                            names = c("coiled-coil")))

sample.gr$SNPsideID <- "top"
sample.gr$SNPsideID[23:29] <- "bottom"

sample.gr$legend <- "COSMIC database"
sample.gr$legend[23:29] <- "LFS-OBs-derived mouse tumors"

features$fill <- c("#FF8833")
features$height <- c(0.05)
sample.gr$color <- "#4C9900"
sample.gr$color[23:29] <- "#0080FF"
sample.gr$border <- "gray30"
sample.gr$alpha <- 0.6
sample.gr$score <- WTAP_cosmic_mutation[,5]
xaxis <- c(1, 100, 200, 300, 396)
sample.gr$dashline.col <- "gray30"

lolliplot(sample.gr, 
          features, 
          ranges = GRanges("chr1", IRanges(1, 396)),
          xaxis = xaxis,
          yaxis = FALSE,
          type="circle",
          cex = 1.2,
          jitter="label",
          label_on_feature = FALSE,
          legend="legend")
png("1.2.png", width=1433, height = 765)
dev.off()


