setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/YTHDF2_lollipop_plot")
getwd()
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("trackViewer")
library(trackViewer)
library(dplyr)

df <- read.csv(file = "YTHDF2_mutaiton_in_COSMIC_patients_and_T2_T5.csv")
df <- df %>% 
  filter(Type %in% c("Unknown","Substitution - coding silent") == FALSE)
SNP <- as.numeric(df[,1])
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, 
                                     names=df[,3]))
features <- GRanges("chr1", 
                    IRanges(c(101, 411), 
                            width = c(100, 135),
                            names = c("C1BD", "YTH")))

sample.gr$SNPsideID <- "top"
sample.gr$SNPsideID[163:169] <- "bottom"

sample.gr$legend <- "COSMIC database"
sample.gr$legend[163:169] <- "LFS-OBs-derived mouse tumors"

features$fill <- c("#FF8833", "#51C6E6")
sample.gr$color <- "#4C9900"
sample.gr$color[163:169] <- "#0080FF"
sample.gr$border <- "gray30"
sample.gr$alpha <- 0.6
features$height <- c(0.03, 0.03)
sample.gr$score <- df[,5]
xaxis <- c(1, 100, 200, 300, 400, 500, 579)
yaxis <- c(0, 3, 6, 9)
sample.gr$dashline.col <- FALSE


plot_lolipop <- lolliplot(sample.gr, 
          features, 
          ranges = GRanges("chr1", IRanges(1, 579)),
          xaxis = xaxis,
          yaxis = yaxis,
          type="circle",
          cex = 0.4,
          jitter="label",
          ylab="Count",
          yaxis.gp = gpar(fontsize=3, lwd=1),
          xaxis.gp = gpar(fontsize=3, lwd=1),
          label_on_feature = TRUE,
          legend = NULL)
png("testplottt22.png", width=1433, height = 765)
dev.off()


png("YTHDF2_lolipop_patient.png",
    width = 4500,
    height = 3000,
    res = 600)
plot_lolipop
dev.off()


