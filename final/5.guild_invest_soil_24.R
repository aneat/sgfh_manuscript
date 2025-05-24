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

##perform cwm##
cwm <- makecwm(new_otu, tra, na.rm = TRUE)
cwm <- subset(cwm)
head(cwm)

##combine per.table and cwm
cwm.per <- merge(cwm, per.table, by = 0)
cwm.per <- column_to_rownames(cwm.per, "Row.names")
env1 <- subset(env1, select = -c(4,6,10,34))
cwm.per <- merge(cwm.per, env1, by = 0)
cwm.per <- column_to_rownames(cwm.per, "Row.names")
saveRDS(cwm.per, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\cwm.per.soil_24.rds")

##find counts associated with each guild
guild_counts_soil <- colSums(tra, na.rm = FALSE, dims = 1)

guild_counts_soil <- as.data.frame(guild_counts_soil)
write.csv(guild_counts_soil, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild_counts_soil_24.csv")

##guild linear regression analyses##
##PATHOGENS##
p <- lm(pathogen ~ species + pH + soil + climate, data = cwm.per)
p1 <- lm(pathogen ~ species + climate + soil, data = cwm.per)
p2 <- lm(pathogen ~ climate + soil, data = cwm.per)
p3 <- lm(pathogen ~ climate, data = cwm.per)
summary(p3)
resid <- resid(p)
models_p <- list(p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
model_avg <- summary(model.avg(models_p, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\soil\\soil.pathogens_24.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\soil\\soil.pathogens.aic_24.csv")

##SAPROBES##
s <- lm(saprobe ~ species + climate + pH + soil, data = cwm.per)
s1 <- lm(saprobe ~ species + pH + soil, data = cwm.per)
s2 <- lm(saprobe ~ pH + soil, data = cwm.per)
s3 <- lm(saprobe ~ pH, data = cwm.per)
summary(s3)
models_s <- list(s, s1, s2, s3)
aic_s <- aictab(cand.set = models_s)
model_avg <- summary(model.avg(models_s, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\soil\\soil.saprobes_24.csv")
write.csv(aic_s, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\soil\\soil.saprobes.aic_24.csv")

##ECM##
ecm <- lm(ecm ~ species + climate + pH + soil, data = cwm.per)
ecm1 <- lm(ecm ~ species + climate + soil, data = cwm.per)
ecm2 <- lm(ecm ~ species + climate, data = cwm.per)
ecm3 <- lm(ecm ~ species, data = cwm.per)
models_ecm <- list(ecm, ecm1, ecm2, ecm3)
aic_ecm <- aictab(cand.set = models_ecm)
aic_ecm
summary(ecm)
model_avg <- summary(model.avg(models_ecm, subset = delta <= 2))
model_avg
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\soil\\soil.ecm_24.csv")
write.csv(aic_ecm, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\soil\\aic.ecm_24.csv")

##AM##
am <- lm(am ~ species + climate + pH + soil, data = cwm.per)
am1 <- lm(am ~ species + pH + soil, data = cwm.per)
am2 <- lm(am ~ species + pH, data = cwm.per)
am3 <- lm(am ~ pH, data = cwm.per)
models_am <- list(am, am1, am2, am3)
aic_am <- aictab(cand.set = models_am)
aic_am
summary(am1)
model_avg <- summary(model.avg(models_am, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
model_avg
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\soil\\soil.am_24.csv")
write.csv(aic_am, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\soil\\soil.am.aic_24.csv")

