setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/mutation_rate_across_cancer_type/tsv_files")
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


setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/mutation_rate_across_cancer_type")
getwd()

allele_frequency <- 
  read.csv(file = "T2_T4_T5_allele_frequency_sum.csv")

top_30 <- allele_frequency[-16,] %>% 
  head(30)
gene_order <- top_30$Gene 
gene_order_2 <- c("YTHDF2", gene_order)

clean_df_t <- as.data.frame(t(clean_df))
clean_df_t_tmp = 
  clean_df_t[-1,]
colnames(clean_df_t_tmp) <- clean_df_t[1,]
clean_df_t_tmp <- 
  tibble::rownames_to_column(clean_df_t_tmp, "Project")
clean_df_t_tmp <- clean_df_t_tmp %>% 
  arrange(factor(Project, levels = gene_order_2))

library(pheatmap)
library(RColorBrewer)

clean_df_t_tmp_tmp = 
  clean_df_t_tmp[,c(2:34)]
row.names(clean_df_t_tmp_tmp) = 
  clean_df_t_tmp[,1]

clean_df_t_tmp_tmp <- 
  type.convert(clean_df_t_tmp_tmp, as.is = TRUE)

pheatmap(clean_df_t_tmp_tmp,
         color = colorRampPalette(brewer.pal(n =9, 
                                             name = "Oranges"))(100),
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         border_color = "#A0A0A0",
         na_col = "#E0E0E0",
         show_rownames = FALSE,
         show_colnames = FALSE,
         fontsize_row = 5)

write.csv(clean_df_t_tmp_tmp, file = "clean_df_t_tmp_tmp.csv")
