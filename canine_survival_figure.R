setwd("C:/Users/dannyj/Documents/Rproject/RJZ_project/canine_database_figure")
getwd()
install.packages("ranger")
install.packages("ggfortify")
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
ATP2A2_survival_input <- read.csv("ATP2A2_survival_input.csv")
# Kaplan Meier Survival Curve
km <- with(ATP2A2_survival_input, Surv(time, status))
head(km,80)

km_fit <- survfit(Surv(time, status) ~ 1, 
                  data=ATP2A2_survival_input)
summary(km_fit, 
        times = c(1,10,20,30,40,50,60,70*(1:10)))
autoplot(km_fit)

km_group_fit <- survfit(Surv(time, status) ~ group, 
                        data=ATP2A2_survival_input)
autoplot(km_group_fit)

cox <- coxph(Surv(time, status) ~ group, 
             data = ATP2A2_survival_input)
summary(cox)
