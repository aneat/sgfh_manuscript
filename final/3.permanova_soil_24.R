##read in packages
library(tidyverse)
library(factoextra)
library(AICcPermanova)
library(vegan)
library(ecole)
library(Hmisc)
library(corrplot)
library(ggpubr)

###FOR SOILS###
##read in matrices
new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_soil_24.rds")
env1 <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_soil_24.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_soil_24.rds")
dist <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\dist_soil_24.rds")
#per.table <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.soil_24.rds")

##subset data frames for soil, foliar, and climate PCA. DO NOT INCLUDE SULFUR. IT LOOKS LIKE A DATA MISSENTRY.
soil.characteristics <- select(env1, c(27:44))
rownames(soil.characteristics) <- env1$row.names
foliar.nutrients <- select(env1, c(13:26))
rownames(foliar.nutrients) <- env1$row.names
climate <- select(env1, c(47:52))
rownames(climate) <- env1$row.names

##soil characteristics pca
data.pca <- princomp(soil.characteristics)
summary(data.pca)
soil.pca.scores <- data.pca$scores
soil.pca.scores <- as.data.frame(soil.pca.scores)
data.pca$loadings[, 1:2]
soil1 <- (fviz_pca_var(data.pca, col.var = "#993404", labelsize = 5, title = "", arrowsize = 1.5, repel = TRUE, ggtheme = theme_classic(base_size = 20)))
soil2 <- (fviz_cos2(data.pca, choice = "var", axes = 1, ggtheme = theme_classic(base_size=30), labelsize = 8, fill = "#993404", color = "black", title = FALSE))
soil1
##ggarrange(soil1, soil2, labels = "AUTO", font.label = list(size = 30))

##foliar nutrients pca
foliar.pca <- princomp(foliar.nutrients)
summary(foliar.pca)
foliar.pca.scores <- foliar.pca$scores
foliar.pca.scores <- as.data.frame(foliar.pca.scores)
foliar.pca$loadings[, 1:2]
ned1 <- (fviz_pca_var(foliar.pca, col.var = "#006D2C", ggtheme = theme_classic(base_size = 20), arrowsize = 1.5, labelsize = 5, title = "", repel = TRUE))
ned2 <- (fviz_cos2(foliar.pca, choice = "var", axes = 1, ggtheme = theme_classic(base_size = 30), labelsize = 8, fill = "#006D2C", color = "black", title = FALSE))
ned1
##ned2
##ggarrange(ned1, ned2, labels = "AUTO", font.label = list(size = 30))

##climate pca
climate.pca <- princomp(climate)
summary(climate.pca)
climate.pca.scores <- climate.pca$scores
climate.pca.scores <- as.data.frame(climate.pca.scores)
climate.pca$loadings[, 1:2]
c1 <- (fviz_pca_var(climate.pca, col.var = "royalblue", ggtheme = theme_classic(base_size = 20), labelsize = 5, title = "", arrowsize = 1.5, repel = TRUE))
c2 <- (fviz_cos2(climate.pca, choice = "var", axes = 1, ggtheme = theme_classic(base_size = 30), title = FALSE, labelsize = 8))
#ggarrange(c1, c2, labels = "AUTO", font.label = list(size = 30))

ggarrange(c1,soil1, ned1, cl.ele, labels = c("A", "B", "C", "D"))
ggarrange(cl.ele,cor.plot, labels = c("A", "B"))

##create permanova table
soil.pca.scores <- select(soil.pca.scores, c(1))
climate.pca.scores <- select(climate.pca.scores, c(1))
foliar.pca.scores <- select(foliar.pca.scores, c(1))
per.table <- select(env1, c(5,6,10,34,4))
per.table <- merge(per.table, soil.pca.scores, 
                   by = 'row.names', all = TRUE)
per.table <- column_to_rownames(per.table, "Row.names")
per.table <- merge(per.table, foliar.pca.scores, 
                   by = 'row.names', all = TRUE) 
per.table <- column_to_rownames(per.table, "Row.names")
per.table <- merge(per.table, climate.pca.scores, 
                   by = 'row.names', all = TRUE)
per.table <- column_to_rownames(per.table, "Row.names")
per.table$soil <- per.table$Comp.1.x
per.table$foliar <- per.table$Comp.1.y
per.table$climate <- per.table$Comp.1

per.table <- select(per.table, -c(6:8))
saveRDS(per.table, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.soil_24.rds")

corr_data <- subset(per.table, select = c(3,4,6:8))
p_values <- rcorr(as.matrix(corr_data))
print(p_values)
c <- cor(corr_data)
cor.plot <- corrplot(c, method = "number", type = "lower", tl.col = "black")
cl.ele <- ggplot(per.table, aes(x = climate, y = elevation)) + stat_smooth(method = "lm", size = 3, color = "royalblue") + theme_classic(base_size = 20) + xlab("climate PCA axis 1") + ylab("elevation (meters)") + annotate("text", x=-2, y=500, label= "r = 0.96", size = 5)
cl.ele

ph.ele <- ggplot(per.table, aes(x = pH, y = elevation)) + stat_smooth(method = "lm", size = 3, color = "brown") + theme_classic(base_size = 20) + xlab("soil pH") + ylab("elevation (meters)") + annotate("text", x=4.3, y=200, label= "r = 0.33", size = 5)
ph.ele

soil.ele <- ggplot(per.table, aes(x = soil, y = elevation)) + stat_smooth(method = "lm", size = 3, color = "brown") + theme_classic(base_size = 20) + xlab("soil PCA axis 1") + ylab("elevation (meters)") + annotate("text", x=-60, y=50, label= "r = 0.29", size = 5)
soil.ele
##permanvoa
a <- adonis2(dist ~ species + climate + soil + pH, data=per.table, by='margin')
b <- adonis2(dist ~ species + climate + pH, data=per.table, by='margin')
c <- adonis2(dist ~ species + climate, data=per.table, by='margin')
d <- adonis2(dist ~ climate, data = per.table, by= "margin")
saveRDS(b, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\permanova_needle_24.rds")

##determine bestfit model. The lower the AIC value, the better the model fit.

AICc_permanova2(a)
AICc_permanova2(b)
AICc_permanova2(c)
AICc_permanova2(d)
soil_perm_tab <- d
write.csv(soil_perm_tab, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.perm_24.csv")
c
##THE END##

