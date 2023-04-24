setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/m6A_mutation_rate_across_cancer_type/tsv_files")
getwd()
library(dplyr)
filenames  <- list.files()

library(data.table)
all_files <- lapply(filenames,function(x) {
  fread(x,data.table=FALSE)
})

names(all_files) <- filenames
filelist <- list(filenames)
combined_data <- Reduce(
  function(...) {
    full_join(..., 
              by = c("Project" = "Project"),
              keep = FALSE)
  },
  all_files
)

filtered_data <- combined_data %>% 
  filter(grepl('TCGA', Project)) %>% 
  select(Project, starts_with("# SSM"))

colnames(filtered_data) <- c("Project", filenames)

filtered_data_2 <- gsub("^.[:punct:]","",filtered_data)

clean_df <- purrr::map_df(filtered_data, 
                          ~ gsub("^.*\\(", "", .x))
clean_df <- purrr::map_df(clean_df, 
                          ~ gsub("%\\)", "", .x))

clean_df <- type.convert(clean_df, as.is = TRUE)

colnames(clean_df) <- 
  gsub("\\.tsv","",colnames(clean_df))


setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/m6A_mutation_rate_across_cancer_type")
getwd()

m6A_gene_order <- 
  read.csv(file = "m6A_gene_order.csv")
gene_order <- m6A_gene_order$Gene 



clean_df_t <- as.data.frame(t(clean_df))


clean_df_t_tmp = 
  clean_df_t[-1,]
colnames(clean_df_t_tmp) <- clean_df_t[1,]
clean_df_t_tmp <- 
  tibble::rownames_to_column(clean_df_t_tmp, "Project")
clean_df_t_tmp <- clean_df_t_tmp %>% 
  arrange(factor(Project, levels = gene_order))

library(pheatmap)
library(RColorBrewer)


# add TARGET OS data manually
# clean_df_t_tmp <-
#   read.csv(file=
#              "clean_df_t_tmp_tmp_manually_add_TARGETG_OS.csv",
#            check.names=FALSE)


clean_df_t_tmp_tmp = 
  clean_df_t_tmp[,c(2:34)]
row.names(clean_df_t_tmp_tmp) = 
  clean_df_t_tmp[,1]

clean_df_t_tmp_tmp <- 
  type.convert(clean_df_t_tmp_tmp, as.is = TRUE)

plot_heatmap <- pheatmap(clean_df_t_tmp_tmp,
         color = colorRampPalette(brewer.pal(n =5, 
                                             name = "RdPu"))(20),
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color = "#A0A0A0",
         na_col = "#E0E0E0",
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 7,
         fontsize_col = 7,
         gaps_row = c(6,8))

png("m6A_mutation_rate_across_cancer_type_heatmap_add_target.png",
    width = 5000,
    height = 3000,
    res = 600)
plot_heatmap
dev.off()

write.csv(clean_df_t_tmp_tmp, file = "clean_df_t_tmp_tmp.csv")
