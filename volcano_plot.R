setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/volcano_plot")
getwd()
df <- read.csv(file = "Annotated_DESeq2_results.csv")

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('EnhancedVolcano')
library(ggplot2)
library(EnhancedVolcano)

keyvals <- ifelse(
  df$log2.FC. < -1 & df$P.adj<0.05, 'royalblue',
  ifelse(df$log2.FC. > 1 & df$P.adj<0.05, 'gold',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'gold'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'

plot_volcano <- EnhancedVolcano(df,
                lab = df$Gene.name,
                x = 'log2.FC.',
                y = 'P.adj',
                xlim = c(min(df[['log2.FC.']], na.rm = TRUE) - 1.5, max(df[['log2.FC.']], na.rm = TRUE) +
                           1.5),
                ylim = c(-5, max(-log10(df[['P.adj']]), na.rm = TRUE) + 5),
                selectLab = c('CHRNA1','IGSF5', 'SCN2A',
                              'DNER','SIDT1','IFITM10',
                              'GPR115','C5AR2','CDH6','ENPP5',
                              'EHD1','CDC42EP2',
                              'MSMO1',
                              'CENPE',	'NEK2',	'KIF15',
                              'FAM111B',	'KIF20A',	
                              'WDR76',	'PTGS2',	'C1orf110',
                              'CCDC15',	'MKI67',	'ESCO2',
                              'ARC',	'UGT8',	
                              'SAPCD2',	'CCNA2'),
                title = 'shYTHDF2 versus shCtrl',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 3.8,
                colCustom = keyvals,
                legendPosition = 'none',
                drawConnectors = TRUE,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = FALSE,
                max.overlaps = Inf,
                border = 'full')
png("Y2_KD_RNA-seq_volcano_plot.png",
    width = 4500,
    height = 3000,
    res = 600)
plot_volcano
dev.off()


