setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/word_cloud")
getwd()
library(dplyr)

#generate dataframe containing word and frequency column------------------
T2 <- read.csv(file = "LFS-T2-MSC.anno_YTC.csv")
T2_filtered <- T2 %>% 
  filter(ExonicFunction %in% c("nonframeshift_substitution",
                               "frameshift_substitution",
                               "nonsynonymous_SNV",
                               "stopgain",
                               "stoploss")) %>% 
  dplyr::filter(!grepl(",",Gene))
T2_allele_frequency_sum <- 
  T2_filtered %>% 
  group_by(Gene) %>% 
  summarise(allele_frequency_sum = sum(tumor)) %>% 
  arrange(desc(allele_frequency_sum))
write.csv(T2_allele_frequency_sum, file = "T2_allele_frequency_sum.csv")


T4 <- read.csv(file = "LFS-T4-MSC.anno_YTC.csv")
T4_filtered <- T4 %>% 
  filter(ExonicFunction %in% c("nonframeshift_substitution",
                               "frameshift_substitution",
                               "nonsynonymous_SNV",
                               "stopgain",
                               "stoploss")) %>% 
  dplyr::filter(!grepl(",",Gene))

T4_allele_frequency_sum <- 
  T4_filtered %>% 
  group_by(Gene) %>% 
  summarise(allele_frequency_sum = sum(tumor)) %>% 
  arrange(desc(allele_frequency_sum))
write.csv(T4_allele_frequency_sum, file = "T4_allele_frequency_sum.csv")


T5 <- read.csv(file = "LFS-T5-MSC.anno_YTC.csv")
T5_filtered <- T5 %>% 
  filter(ExonicFunction %in% c("nonframeshift_substitution",
                               "frameshift_substitution",
                               "nonsynonymous_SNV",
                               "stopgain",
                               "stoploss")) %>% 
  dplyr::filter(!grepl(",",Gene))

T5_allele_frequency_sum <- 
  T5_filtered %>% 
  group_by(Gene) %>% 
  summarise(allele_frequency_sum = sum(tumor)) %>% 
  arrange(desc(allele_frequency_sum))
write.csv(T5_allele_frequency_sum, file = "T5_allele_frequency_sum.csv")


#generate wordcloud-------------------------------------------------

#install.packages("wordcloud")
library(wordcloud)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("wordcloud2")
library(wordcloud2)

# T2_allele_frequency_sum <- read.csv(file = "T2_allele_frequency_sum.csv")
# T4_allele_frequency_sum <- read.csv(file = "T4_allele_frequency_sum.csv")
# T5_allele_frequency_sum <- read.csv(file = "T5_allele_frequency_sum.csv")

T2_allele_frequency_sum_top <- T2_allele_frequency_sum %>% 
  head(30)
wordcloud2(data = T2_allele_frequency_sum_top, 
                      size = 0.3,
                      fontFamily = 'Arial',
                      shuffle = 'FALSE',
                      color = 'random-dark',
                      shape = 'circle')

T4_allele_frequency_sum_top <- T4_allele_frequency_sum %>% 
  head(30)
wordcloud2(data = T4_allele_frequency_sum_top, 
           size = 0.2,
           fontFamily = 'Arial',
           shuffle = 'FALSE',
           color = 'random-dark',
           shape = 'circle')


T5_allele_frequency_sum_top <- T5_allele_frequency_sum %>% 
  head(30)
wordcloud2(data = T5_allele_frequency_sum_top, 
           size = 0.3,
           fontFamily = 'Arial',
           shuffle = 'FALSE',
           color = 'random-dark',
           shape = 'circle')






setwd("C:/Users/danny/OneDrive - 中國醫藥大學/文件/R_project/RJZ_project/word_cloud")
getwd()
library(dplyr)
library(wordcloud2)

T5_allele_frequency_sum_top <- T5_allele_frequency_sum %>% 
  head(30)
wordcloud2(data = T5_allele_frequency_sum_top, 
           size = 0.3,
           fontFamily = 'Arial',
           shuffle = 'FALSE',
           color = 'random-dark',
           shape = 'circle')


#wordcloud package-----------------------------------------------------------
# set.seed(1234) # for reproducibility
# wordcloud(words = T2_allele_frequency_sum$Gene,
#           freq = T2_allele_frequency_sum$allele_frequency_sum,
#           min.freq = 0,
#           max.words = 30,
#           scale = c(1.5,.15),
#           random.order = FALSE,
#           rot.per = 0,
#           colors = brewer.pal(8, "Dark2"))
