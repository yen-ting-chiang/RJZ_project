setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/world_cloud_by_allele_frequency")
getwd()
library(dplyr)

T2 <- read.csv(file = "LFS-T2-MSC.anno_YTC.csv")
T2_filtered <- T2 %>% 
  filter(ExonicFunction %in% c("nonframeshift_substitution",
                               "frameshift_substitution",
                               "nonsynonymous_SNV",
                               "stopgain",
                               "stoploss"))
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
                               "stoploss"))
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
                               "stoploss"))
T5_allele_frequency_sum <- 
  T5_filtered %>% 
  group_by(Gene) %>% 
  summarise(allele_frequency_sum = sum(tumor)) %>% 
  arrange(desc(allele_frequency_sum))
write.csv(T5_allele_frequency_sum, file = "T5_allele_frequency_sum.csv")



install.packages("wordcloud")
library(wordcloud)
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("wordcloud2")
library(wordcloud2)


set.seed(1234) # for reproducibility
wordcloud(words = df$word, freq = df$freq, min.freq = 1,
          max.words=200, random.order=FALSE, rot.per=0.35,
          colors=brewer.pal(8, "Dark2"))

