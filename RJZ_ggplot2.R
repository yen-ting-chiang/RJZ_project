setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/RJZ_ggplot2")
getwd()
library(dplyr)
library(ggplot2)

allele_frequency_exonic_nonsynonymous_2 <- 
  read.csv("allele_frequency_exonic_nonsynonymous_2.csv")
ggplot(data = allele_frequency_exonic_nonsynonymous_2) + 
  geom_density(aes(x = allele_frequency, 
                     fill = sample_name),
                 alpha = 0.8)

ggplot(data = allele_frequency_exonic_nonsynonymous_2) + 
  geom_histogram(aes(x = allele_frequency, 
                   fill = sample_name)) + 
  facet_grid(cols = vars(sample_name))

ggplot(data = allele_frequency_exonic_nonsynonymous_2) + 
  geom_violin(aes(x = sample_name, 
                  y = allele_frequency,
                  fill = sample_name))

