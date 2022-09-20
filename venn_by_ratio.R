
install.packages("eulerr")
library(eulerr)

OB_T2 <- euler(c("A" = 442, "B" = 570, "A&B" = 39))
plot(OB_T2, quantities = FALSE, alpha=0.5,labels = FALSE,
     fill=c("#E0E0E0", "#66B2FF"), col = c("#202020", "#003366"))

OB_T4 <- euler(c("A" = 401, "B" = 488, "A&B" = 41))
plot(OB_T4, quantities = FALSE, alpha=0.5,labels = FALSE,
     fill=c("#E0E0E0", "#CCFF99"), col = c("#003366", "#336600"))

OB_T5 <- euler(c("A" = 359, "B" = 1219, "A&B" = 83))
plot(OB_T5, quantities = FALSE, alpha=0.5,labels = FALSE,
     fill=c("#E0E0E0", "#FFCC99"), col = c("#003366", "#663300"))


T2_T4_T5 <- euler(c("A" = 395, "B" = 419, "C" = 1078, 
                    "A&B" = 25, "B&C" = 74, 
                    "A&C" = 139, "A&B&C" = 11), shape = "circle")


plot(T2_T4_T5, quantities = FALSE, alpha=0.5,labels = FALSE,
     fill=c("#66B2FF", "#CCFF99","#FFCC99"), 
     col = c("#003366", "#336600","#663300"))



article_venn <- euler(c("A" = 503, "B" = 193, "C" = 664, 
                    "A&B" = 23, "B&C" = 10, 
                    "A&C" = 39, "A&B&C" = 6), shape = "ellipse")

plot(article_venn, quantities = FALSE, alpha=0.5,labels = FALSE,
     fill=c("#66CC00", "#FFFF00","#7F00FF"), 
     col = c("#000000", "#000000","#000000"),
     lwd = c(2,2,2))


# VennDiag <- euler(c("A" = 481, "B" = 180, "C" = 648, "D" = 240 ,
#                     "A&B" = 18, "A&C" = 38, "A&D" = 22,
#                     "B&C" = 8, "B&D" = 13, "C&D" = 16,
#                     "A&B&C" = 6, "A&B&D" = 5,"A&C&D" = 1,"B&C&D" = 2,
#                     "A&B&C&D" = 6), shape = "circle")
# plot(VennDiag, quantities = TRUE, font=0.1, cex=0.1, alpha=0.5, 
#      labels = TRUE)

