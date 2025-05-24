##load packages##
library(ggpubr)
library(gridExtra)
library(broom)
library(tidyverse)
library(AICcmodavg)
library(MuMIn)

##load data
resid_soil_hc <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\resid_soil_hc_24.rds")
resid_foliar_hc <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\resid_foliar_hc_24.rds")

##SOIL GRAPHS##
##SOIL DIVERSITY##
p_climate <- lm(p_resid_d ~ climate + species + climate:species, data = resid_soil_hc)
summary(p_climate)
spd <- ggplot(resid_soil_hc, aes(x = climate, y = p_resid_d, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species), linetype = "twodash") + ylab("shannon's diversity") + theme_classic(base_size=15) + scale_color_brewer(palette = "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.11", size = 5, vjust = -0.5, hjust = -0.2) + ggtitle("soil pathogens") + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) 
spd

e_climate <- lm(e_resid_d ~ climate + species + climate:species, data = resid_soil_hc)
summary(e_climate)
sed <- ggplot(resid_soil_hc, aes(x = climate, y = e_resid_d, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species)) + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p < 0.01", size = 5, vjust = -0.5, hjust = -0.2) + ggtitle("emf") + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.4))
sed

a_climate <- lm(a_resid_d ~ climate + species + climate:species, data = resid_soil_hc)
summary(a_climate)
sad <- ggplot(resid_soil_hc, aes(x = climate, y = a_resid_d, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species)) + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p < 0.01", size = 5, vjust = -0.5, hjust = -0.2) + ggtitle("am") + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.4))
sad

s_climate <- lm(s_resid_d ~ climate + species + climate:species, data = resid_soil_hc)
summary(s_climate)
ssd <- ggplot(resid_soil_hc, aes(x = climate, y = s_resid_d, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species), linetype = "twodash") + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.77", size = 5, vjust = -0.5, hjust = -0.2) + ggtitle("saprobe") + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.4))
ssd

##SOIL CWM##
p_climate <- lm(p_resid_cwm ~ climate + species + climate:species, data = resid_soil_hc)
summary(p_climate)
spc <- ggplot(resid_soil_hc, aes(x = climate, y = p_resid_cwm, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species), linetype = "twodash") + ylab("community weighted mean") + theme_classic(base_size=15) + scale_color_brewer(palette = "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.29", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) 
spc

e_climate <- lm(e_resid_cwm ~ climate + species + climate:species, data = resid_soil_hc)
summary(e_climate)
sec <- ggplot(resid_soil_hc, aes(x = climate, y = e_resid_cwm, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species), linetype = "twodash") + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.13", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.4))
sec

a_climate <- lm(a_resid_cwm ~ climate + species + climate:species, data = resid_soil_hc)
summary(a_climate)
sac <- ggplot(resid_soil_hc, aes(x = climate, y = a_resid_cwm, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species)) + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.04", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("")
sac

s_climate <- lm(s_resid_cwm ~ climate + species + climate:species, data = resid_soil_hc)
summary(s_climate)
ssc <- ggplot(resid_soil_hc, aes(x = climate, y = s_resid_cwm, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species), linetype = "twodash") + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.51", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.4))
ssc

##SOIL SPECIFICITY##
p_climate <- lm(p_resid_spe ~ climate + species + climate:species, data = resid_soil_hc)
summary(p_climate)
sps <- ggplot(resid_soil_hc, aes(x = climate, y = p_resid_spe, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species), linetype = "twodash") + ylab("shannon's diversity") + theme_classic(base_size=15) + scale_color_brewer(palette = "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.15", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) 
sps

e_climate <- lm(e_resid_spe ~ climate + species + climate:species, data = resid_soil_hc)
summary(e_climate)
ses <- ggplot(resid_soil_hc, aes(x = climate, y = e_resid_spe, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species)) + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.05", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("")
ses

a_climate <- lm(a_resid_spe ~ climate + species + climate:species, data = resid_soil_hc)
summary(a_climate)
sas <- ggplot(resid_soil_hc, aes(x = climate, y = a_resid_spe, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species), linetype = "twodash") + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.16", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.4))
sas

s_climate <- lm(s_resid_spe ~ climate + species + climate:species, data = resid_soil_hc)
summary(s_climate)
sss <- ggplot(resid_soil_hc, aes(x = climate, y = s_resid_spe, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species)) + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.04", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("")
sss

#FOLIAR GRAPHS##
##foliar DIVERSITY##
p_climate <- lm(p_resid_d ~ climate + species + climate:species, data = resid_foliar_hc)
summary(p_climate)
fpd <- ggplot(resid_foliar_hc, aes(x = climate, y = p_resid_d, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species), linetype = "twodash") + ylab(" ") + theme_classic(base_size=15) + scale_color_brewer(palette = "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.21", size = 5, vjust = -0.5, hjust = -0.2) + ggtitle("foliar pathogens") + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + xlab(" ")
fpd

s_climate <- lm(s_resid_d ~ climate + species + climate:species, data = resid_foliar_hc)
summary(s_climate)
fsd <- ggplot(resid_foliar_hc, aes(x = climate, y = s_resid_d, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species)) + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.03", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("")
fsd

##FOLIAR CWM##
p_climate <- lm(p_resid_cwm ~ climate + species + climate:species, data = resid_foliar_hc)
summary(p_climate)
fpc <- ggplot(resid_foliar_hc, aes(x = climate, y = p_resid_cwm, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species), linetype = "twodash") + theme_classic(base_size=15) + scale_color_brewer(palette = "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.15", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + xlab("") + ylab("")
fpc

s_climate <- lm(s_resid_cwm ~ climate + species + climate:species, data = resid_foliar_hc)
summary(s_climate)
ssc <- ggplot(resid_foliar_hc, aes(x = climate, y = s_resid_cwm, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species)) + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p < 0.01", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + xlab(" ")
ssc

##FOLIAR SPECIFICITY##
p_climate <- lm(p_resid_spe ~ climate + species + climate:species, data = resid_foliar_hc)
summary(p_climate)
fps <- ggplot(resid_foliar_hc, aes(x = climate, y = p_resid_spe, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species), linetype = "twodash") + ylab("") + theme_classic(base_size=15) + scale_color_brewer(palette = "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.15", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + xlab(" ")
fps

s_climate <- lm(s_resid_spe ~ climate + species + climate:species, data = resid_foliar_hc)
summary(s_climate)
fss <- ggplot(resid_foliar_hc, aes(x = climate, y = s_resid_spe, color = species)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = species)) + theme_classic(base_size=15) + scale_color_brewer(palette= "Dark2") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p < 0.01", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + xlab(" ")
sss

plots <- ggarrange(spd, sed, sad, fpd, spc, sec, sac, fpc, sps, ses, sas, fps,
                   ncol = 4, nrow = 3, align = c("hv"), common.legend = TRUE, legend = "right")
plots




