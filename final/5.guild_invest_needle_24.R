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

##perform cwm##
cwm <- makecwm(new_otu, tra, na.rm = TRUE)
cwm <- subset(cwm, select = -c(8:9, 11:12,18:19))
head(cwm)

##combine per.table and cwm
cwm.per <- merge(cwm, per.table, by = 0)
cwm.per <- column_to_rownames(cwm.per, "Row.names")
env1 <- subset(env1, select = -c(4,6,10,34))
cwm.per <- merge(cwm.per, env1, by = 0)
cwm.per <- column_to_rownames(cwm.per, "Row.names")
saveRDS(cwm.per, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\cwm.per.needle_24.rds")

##find counts associated with each guild
guild_counts_needle <- colSums(tra, na.rm = FALSE, dims = 1)

guild_counts_needle <- as.data.frame(guild_counts_needle)
write.csv(guild_counts_needle, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild_counts_needle_24.csv")

##guild linear regression analyses##
##PATHOGENS##
p <- lm(pathogen ~ species + climate + pH + foliar, data = cwm.per)
p1 <- lm(pathogen ~ species + pH + foliar, data = cwm.per)
p2 <- lm(pathogen ~ pH + climate, data = cwm.per)
p3 <- lm(pathogen ~ climate, data = cwm.per)
AICc(p3)
summary(p3)
models_p <- list(p, p1, p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
model_avg <- summary(model.avg(models_p, subset = delta <= 2))
model_avg
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\needles\\foliar.pathogens_24.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\needles\\foliar.pathogens.aic_24.csv")
path <- ggplot(cwm.per, aes(x = elevation, y = pathogen, color = species)) + geom_point()
path

##SAPROBES##
s <- lm(saprobe ~ species + climate + pH + foliar, data = cwm.per)
s1 <- lm(saprobe ~ species + pH + foliar, data = cwm.per)
s2 <- lm(saprobe ~ pH + foliar, data = cwm.per)
s3 <- lm(saprobe ~ foliar, data = cwm.per)
models_s <- list(s, s1, s2, s3)
aictab(cand.set = models_s)
model_avg <- summary(model.avg(models_s, subset = delta <= 2))
model_avg
mod_avg_full <- model_avg$coefmat.full
mod_avg_full
write.csv(mod_avg_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\needles\\foliar.saprobes_24.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild models\\needles\\foliar.saprobes.aic_24.csv")

