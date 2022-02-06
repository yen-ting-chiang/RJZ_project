setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/canine_database_figure")
getwd()
install.packages("ranger")
install.packages("ggfortify")
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
survival_input <- read.csv("survival_data_cluster.csv")
# Kaplan Meier Survival Curve
km <- with(survival_input, Surv(time, status))
head(km,80)

km_fit <- survfit(Surv(time, status) ~ 1, 
                  data=survival_input)
summary(km_fit, 
        times = c(1,10,20,30,40,50,60,70*(1:10)))
autoplot(km_fit)

km_group_fit <- survfit(Surv(time, status) ~ group, 
                        data=survival_input)
autoplot(km_group_fit)

cox <- coxph(Surv(time, status) ~ group, 
             data = survival_input)
summary(cox)
