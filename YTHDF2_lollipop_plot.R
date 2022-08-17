setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/YTHDF2_lollipop_plot")
getwd()
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("trackViewer")
library(trackViewer)
YTHDF2_cosmic_mutation <- read.csv(file = "YTHDF2 mutaiton in COSMIC.csv")
library(dplyr)
YTHDF2_cosmic_mutation <- YTHDF2_cosmic_mutation %>% 
  filter(Mutation.Type %in% c("Unknown","Substitution - coding silent") == FALSE)
SNP <- as.numeric(YTHDF2_cosmic_mutation[,1])
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, 
                                     names=YTHDF2_cosmic_mutation[,3]))
features <- GRanges("chr1", 
                    IRanges(c(101, 411), 
                                    width = c(100, 135),
                                    names = c("C1BD", "YTH")))


features$fill <- c("#FF8833", "#51C6E6")
sample.gr$color <- "#4C9900"
sample.gr$border <- "gray30"
sample.gr$alpha <- 0.6
features$height <- c(0.05, 0.05)
sample.gr$score <- YTHDF2_cosmic_mutation[,5]
xaxis <- c(1, 100, 200, 300, 400, 500, 579)
yaxis <- c(0, 5, 10, 15)
names(yaxis) <- yaxis # define the labels
names(yaxis)[3] <- "y-axis" 
sample.gr$dashline.col <- "gray30"

lolliplot(sample.gr, 
          features, 
          ranges = GRanges("chr1", IRanges(1, 579)),
          xaxis = xaxis,
          yaxis = FALSE,
          type="circle",
          cex = 0.4,
          jitter="label",
          ylab="Count",
          yaxis.gp = gpar(fontsize=5, lwd=0),
          xaxis.gp = gpar(fontsize=7, lwd=0))
