##read in packages
library(tidyverse)
library(vegan)
library(ecole)
library(RColorBrewer)
library(ComplexHeatmap)
library(MuMIn)
library(AICcmodavg)

##read in new_otu table, env1, and trait matrix
new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_soil_24.rds")
env1 <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_soil_24.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_soil_24.rds")
per.table <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.soil_24.rds")
cwm.per <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\cwm.per.soil_24.rds")

##subset guilds
path <- subset(tra, pathogen == 1)
sap <- subset(tra, saprobe == 1)
ecm <- subset(tra, ecm == 1)
am <- subset(tra, am == 1)

path_names <- row.names(path)
sap_names <- row.names(sap)
ecm_names <- row.names(ecm)
am_names <- row.names(am)

path_otu <- new_otu[path_names]
sap_otu <- new_otu[sap_names]
ecm_otu <- new_otu[ecm_names]
am_otu <- new_otu[am_names]

x_path <- specnumber(path_otu)
y_path <- diversity(path_otu)

x_sap <- specnumber(sap_otu)
y_sap <- diversity(sap_otu)

x_ecm <- specnumber(ecm_otu)
y_ecm <- diversity(ecm_otu)

x_am <- specnumber(am_otu)
y_am <- diversity(am_otu)

x <- specnumber(new_otu)
y <- diversity(new_otu)

##match order of matrices
cwm.per <- cwm.per[match(row.names(new_otu), row.names(cwm.per)),]

cwm.per$path.diversity <- y_path
cwm.per$path.richness <- x_path

cwm.per$sap.diversity <- y_sap
cwm.per$sap.richness <- x_sap

cwm.per$ecm.diversity <- y_ecm
cwm.per$ecm.richness <- x_ecm

cwm.per$am.diversity <- y_am
cwm.per$am.richness <- x_am

cwm.per$total.diversity <- y
cwm.per$total.richness <- x

saveRDS(cwm.per, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soil.diversity_24.rds")

##investigate##
##GENERAL##

##guild linear regression analyses##
##PATHOGENS##
##diversity##
p <- lm(path.diversity ~ species + pH + soil + climate, data = cwm.per)
p1 <- lm(path.diversity ~ species + climate + soil, data = cwm.per)
p2 <- lm(path.diversity ~ climate + soil, data = cwm.per)
p3 <- lm(path.diversity ~ climate, data = cwm.per)
summary(p2)
resid <- resid(p)
p.resid <- lm(resid ~ climate, data = cwm.per)
summary(p3)
models_p <- list(p, p1, p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
model_avg <- summary(model.avg(models_p, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.path.diversity_24.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.path.diversity_24.aic.csv")

##richness##
p <- lm(path.richness ~ species + climate + pH + soil, data = cwm.per)
p1 <- lm(path.richness ~ species + climate + soil, data = cwm.per)
p2 <- lm(path.richness ~ climate + soil, data = cwm.per)
p3 <- lm(path.richness ~ climate, data = cwm.per)
summary(p)
models_p <- list(p, p1, p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
model_avg <- summary(model.avg(models_p, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.path.richness_24.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.path.richness.aic_24.csv")

##SAPROBES##
##diversity##
s <- lm(sap.diversity ~ species + climate + pH + soil, data = cwm.per)
s1 <- lm(sap.diversity ~ species + climate + soil, data = cwm.per)
s2 <- lm(sap.diversity ~ climate + soil, data = cwm.per)
s3 <- lm(sap.diversity ~ climate, data = cwm.per)
summary(s3)
models_s <- list(s, s1, s2, s3)
aic_s <- aictab(cand.set = models_s)
aic_s
model_avg <- summary(model.avg(models_s, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.sap.diversity_24.csv")
write.csv(aic_s, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.sap.diversity.aic_24.csv")

##richness##
s <- lm(sap.richness ~ species + climate + pH + soil, data = cwm.per)
s1 <- lm(sap.richness ~ species + climate + soil, data = cwm.per)
s2 <- lm(sap.richness ~ climate + soil, data = cwm.per)
s3 <- lm(sap.richness ~ climate, data = cwm.per)
summary(s)
models_s <- list(s, s1, s2, s3)
aic_s <- aictab(cand.set = models_s)
model_avg <- summary(model.avg(models_s, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.sap.richness_24.csv")
write.csv(aic_s, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.sap.richness.aic_24.csv")

##ECM##
##diversity##
ecm <- lm(ecm.diversity ~ species + climate + pH + soil, data = cwm.per)
ecm1 <- lm(ecm.diversity ~ species + climate + soil, data = cwm.per)
ecm2 <- lm(ecm.diversity ~ climate + soil, data = cwm.per)
ecm3 <- lm(ecm.diversity ~ climate, data = cwm.per)
summary(ecm)
models_ecm <- list(ecm, ecm1, ecm2, ecm3)
aic_ecm <- aictab(cand.set = models_ecm)
aic_ecm
model_avg <- summary(model.avg(models_ecm, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.ecm.diversity_24.csv")
write.csv(aic_ecm, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.ecm.diversity.aic_24.csv")

##richness##
ecm <- lm(ecm.richness ~ species + climate + pH + soil, data = cwm.per)
ecm1 <- lm(ecm.richness ~ species + climate + soil, data = cwm.per)
ecm2 <- lm(ecm.richness ~ climate + soil, data = cwm.per)
ecm3 <- lm(ecm.richness ~ climate, data = cwm.per)
summary(ecm)
models_ecm <- list(ecm, ecm1, ecm2, ecm3)
aic_ecm <- aictab(cand.set = models_ecm)
aic_ecm
model_avg <- summary(model.avg(models_ecm, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.ecm.richness_24.csv")
write.csv(aic_ecm, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.ecm.richness.aic_24.csv")

##AM##
##diversity##
am <- lm(am.diversity ~ species + climate + pH + soil, data = cwm.per)
am1 <- lm(am.diversity ~ species + climate + soil, data = cwm.per)
am2 <- lm(am.diversity ~ climate + soil, data = cwm.per)
am3 <- lm(am.diversity ~ climate, data = cwm.per)
summary(am)
models_am <- list(am, am1, am2, am3)
aic_am <- aictab(cand.set = models_am)
aic_am
model_avg <- summary(model.avg(models_am, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.am.diversity_24.csv")
write.csv(aic_am, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.am.diversity.aic_24.csv")

##richness##
am <- lm(am.richness ~ species + climate + pH + soil, data = cwm.per)
am1 <- lm(am.richness ~ species + climate + soil, data = cwm.per)
am2 <- lm(am.richness ~ climate + soil, data = cwm.per)
am3 <- lm(am.richness ~ climate, data = cwm.per)
summary(am)
models_am <- list(am, am1, am2, am3)
aic_am <- aictab(cand.set = models_am)
model_avg <- summary(model.avg(models_am, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.am.richness_24.csv")
write.csv(aic_am, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.am.richness.aic_24.csv")

lm <- lm(am.richness ~ species, data = cwm.per)
summary(lm)

am <- ggplot(cwm.per, aes(x = species, y = am.diversity)) + geom_boxplot()
am

sap <- ggplot(cwm.per, aes(x = elevation, y = am.diversity)) + geom_point()
sap

