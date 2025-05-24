##read in packages
library(tidyverse)
library(vegan)
library(ecole)
library(phyloseq)
library(RColorBrewer)
library(AICcmodavg)
library(MuMIn)

##read in new_otu table, env1, and trait matrix
new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_needle_24.rds")
spe_needle <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\spe_needle_24.rds")
env1 <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_needle_24.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_needle_24.rds")
taxa <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\taxa_needle_24.rds")
per.table <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.needle_24.rds")
cwm.per <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\cwm.per.needle_24.rds")

##subset guilds
path <- subset(tra, pathogen == 1)
sap <- subset(tra, saprobe == 1)

path_names <- row.names(path)
sap_names <- row.names(sap)

path_otu <- new_otu[path_names]
sap_otu <- new_otu[sap_names]

x_path <- specnumber(path_otu)
y_path <- diversity(path_otu)

x_sap <- specnumber(sap_otu)
y_sap <- diversity(sap_otu)

x <- specnumber(new_otu)
y <- diversity(new_otu)

##match order of matrices
cwm.per <- cwm.per[match(row.names(new_otu), row.names(cwm.per)),]

cwm.per$path.diversity <- y_path
cwm.per$path.richness <- x_path

cwm.per$sap.diversity <- y_sap
cwm.per$sap.richness <- x_sap

cwm.per$total.diversity <- y
cwm.per$total.richness <- x

saveRDS(cwm.per,"\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needle.div_24.rds")

##investigate##
##GENERAL##

##diversity##
d <- lm(total.diversity ~ species + climate + pH + soil, data = cwm.per)
d1 <- lm(total.diversity ~ species + climate + soil, data = cwm.per)
d2 <- lm(total.diversity ~ climate + soil, data = cwm.per)
d3 <- lm(total.diversity ~ climate, data = cwm.per)
AICc(d3)
summary(d3)
models_d <- list(d, d1, d2, d3)
aic_d <- aictab(cand.set = models_d)
model_avg <- summary(model.avg(models_d, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.all.div_24.csv")
write.csv(aic_d, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.all.div.aic_24.csv")
all <- ggplot(cwm.per, aes(x = elevation, y = total.diversity)) + geom_point()
all

##richness##
r <- lm(total.richness ~ species + climate + pH + soil, data = cwm.per)
r1 <- lm(total.richness ~ species + climate + soil, data = cwm.per)
r2 <- lm(total.richness ~ climate + soil, data = cwm.per)
r3 <- lm(total.richness ~ climate, data = cwm.per)
AICc(r3)
summary(r3)
models_r <- list(r, r1, r2, r3)
aic_r <- aictab(cand.set = models_r)
model_avg <- summary(model.avg(models_r, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.all.rich_24.csv")
write.csv(aic_r, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.all.div.aic_24.csv")
all <- ggplot(cwm.per, aes(x = elevation, y = total.diversity)) + geom_point()
all

##PATHOGENS##
##diversity##
p <- lm(path.diversity ~ species + climate + pH + soil, data = cwm.per)
p1 <- lm(path.diversity ~ species + climate + soil, data = cwm.per)
p2 <- lm(path.diversity ~ climate + soil, data = cwm.per)
p3 <- lm(path.diversity ~ climate, data = cwm.per)
AICc(p)
summary(p)
models_p <- list(p, p1, p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
model_avg <- summary(model.avg(models_p, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.pathogens.div_24.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.pathogens.div.aic_24.csv")

##richness##
p <- lm(path.richness ~ species + climate + pH + soil, data = cwm.per)
p1 <- lm(path.richness ~ species + climate + soil, data = cwm.per)
p2 <- lm(path.richness ~ climate + soil, data = cwm.per)
p3 <- lm(path.richness ~ climate, data = cwm.per)
AICc(p3)
summary(p3)
models_p <- list(p, p1, p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
model_avg <- summary(model.avg(models_p, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.pathogens.rich_24.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.pathogens.rich.aic_24.csv")

##SAPROBES##
##diversity##
s <- lm(sap.diversity ~ species + climate + pH + soil, data = cwm.per)
s1 <- lm(sap.diversity ~ species + climate + soil, data = cwm.per)
s2 <- lm(sap.diversity ~ climate + soil, data = cwm.per)
s3 <- lm(sap.diversity ~ climate, data = cwm.per)
AICc(s3)
summary(s)
models_s <- list(s, s1, s2, s3)
aic_s <- aictab(cand.set = models_s)
model_avg <- summary(model.avg(models_s, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.sap.div_24.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.sap.div.aic_24.csv")

##richness##
s <- lm(sap.richness ~ species + climate + pH + soil, data = cwm.per)
s1 <- lm(sap.richness ~ species + climate + soil, data = cwm.per)
s2 <- lm(sap.richness ~ climate + soil, data = cwm.per)
s3 <- lm(sap.richness ~ climate, data = cwm.per)
AICc(s3)
summary(s)
models_s <- list(s, s1, s2, s3)
aic_s <- aictab(cand.set = models_s)
model_avg <- summary(model.avg(models_s, subset = delta <= 2))
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.sap.rich_24.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar\\foliar.sap.rich.aic_24.csv")


sap <- ggplot(cwm.per, aes(x = elevation, y = sap.richness)) + geom_point()
sap

