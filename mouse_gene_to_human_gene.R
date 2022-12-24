setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/mouse_gene_to_human_gene")
getwd()
gene_df <- read.csv(file = "gigi_PNAS_mouse_gene_list.csv")
musGenes <- gene_df$Gene

library(dplyr)
mouse_human_genes = read.csv("HOM_MouseHumanSequence.rpt",sep="\t")
convert_mouse_to_human <- function(gene_list){
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  
  return (output)
}
humanGenes <- convert_mouse_to_human(musGenes)
write.csv(humanGenes, file = "gigi_PNAS_humanGenes.csv")
