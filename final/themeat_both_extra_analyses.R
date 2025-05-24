##load packages##
library(ggpubr)
library(gridExtra)
library(broom)
library(tidyverse)
library(AICcmodavg)
library(MuMIn)

##load soil data 
soil.diversity <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soil.diversity.rds")
soil.spec <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soil.spec.rds")

##load foliar data 
foliar.diversity <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needle.div.rds")
foliar.spec <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\foliar.spec.rds")

soil.diversity <- subset(soil.diversity, select = c(1,2,3,8,33,35:39,89:98))

merge_soil <- merge(soil.diversity, soil.spec, by = 0)
merge_soil <- column_to_rownames(merge_soil, "Row.names")
resid_soil <- subset(merge_soil, select = c(5,6,7,8,10))

foliar.diversity <- subset(foliar.diversity, select = c(2,3,27:33,83:88, 43:50))

merge_foliar <- merge(foliar.diversity, foliar.spec, by = 0)
merge_foliar <- column_to_rownames(merge_foliar, "Row.names")
resid_foliar <- subset(merge_foliar, select = c(3,4,6,8,9))

##ECM DIVERSITY##
e2 <- lm(ecm.diversity ~ species + pH, data = merge_soil)
e_ph <- lm(ecm.diversity ~ climate + species + soil, data = merge_soil)
e_sp <- lm(ecm.diversity ~ climate + pH + soil, data = merge_soil)
summary(e_ph)
summary(e_sp)
e_resid_ph <- resid(e_ph)
e_resid_sp <- resid(e_sp)
e_ph_lm <- lm(e_resid_ph ~ pH, data = merge_soil)
e_sp_lm <- lm(e_resid_sp ~ species, data = merge_soil)
summary(e_ph_lm)
summary(e_sp_lm)
e_resid_ph <- as.data.frame(e_resid_ph)
e_resid_sp <- as.data.frame(e_resid_sp)
resid_soil_extra <- cbind(resid_soil, e_resid_sp, e_resid_ph)

##AM DIVERSITY##
am_sp <- lm(am.diversity ~ climate + pH + soil, data = merge_soil)
summary(am_sp)
a_resid_sp <- resid(am_sp)
a_sp_lm <- lm(a_resid_sp ~ species, data = merge_soil)
summary(a_sp_lm)
a_resid_sp <- as.data.frame(a_resid_sp)
resid_soil_extra <- cbind(resid_soil_extra, a_resid_sp)

##AM CWM##
am_sp_cwm <- lm(am ~ climate + pH + soil, data = merge_soil)
summary(am_sp_cwm)
a_resid_cwm <- resid(am_sp_cwm)
a_cwm_lm <- lm(a_resid_cwm ~ species, data = merge_soil)
summary(a_cwm_lm)
a_resid_cwm <- as.data.frame(a_resid_cwm)
resid_soil_extra <- cbind(resid_soil_extra, a_resid_cwm)

##SAPROBE CWM##
sap_cwm <- lm(saprobe ~ climate + species + soil, data = merge_soil)
summary(sap_cwm)
sap_resid_cwm <- resid(sap_cwm)
sap_cwm_lm <- lm(sap_resid_cwm ~ pH, data = merge_soil)
summary(sap_cwm_lm)
sap_resid_cwm <- as.data.frame(sap_resid_cwm)
resid_soil_extra <- cbind(resid_soil_extra, sap_resid_cwm)

##ECM HOST SPECIFICITY##
ecm_hs <- lm(ecm.host.spe ~ climate + species + soil, data = merge_soil)
summary(ecm_hs)
ecm_resid_hs <- resid(ecm_hs)
ecm_hs_lm <- lm(ecm_resid_hs ~ pH, data = merge_soil)
summary(ecm_hs_lm)
ecm_resid_hs <- as.data.frame(ecm_resid_hs)
resid_soil_extra <- cbind(resid_soil_extra, ecm_resid_hs)


##FOLIAR PATHOGEN DIVERSITY##
fpath_ph <- lm(path.diversity ~ climate + species + foliar, data = merge_foliar)
summary(fpath_ph)
fpath_div <- resid(fpath_ph)
fpath_lm <- lm(fpath_div ~ pH, data = merge_foliar)
summary(fpath_lm)
fpath_div <- as.data.frame(fpath_div)
resid_foliar_extra <- cbind(resid_foliar, fpath_div)

##FOLIAR SAPROBE DIVERSITY##
fsap_sp <- lm(sap.diversity ~ climate + pH + foliar, data = merge_foliar)
summary(fsap_sp)
fsap_div <- resid(fsap_sp)
fsap_lm <- lm(fsap_div ~ species, data = merge_foliar)
summary(fsap_lm)
fpath_div <- as.data.frame(fsap_div)
resid_foliar_extra <- cbind(resid_foliar_extra, fsap_div)

##FOLIAR PATHOGEN CWM##
fpath_cwm <- lm(pathogen ~ climate + species + foliar, data = merge_foliar)
summary(fpath_cwm)
fpath_cwm <- resid(fpath_cwm)
fpath_lm <- lm(fpath_cwm ~ pH, data = merge_foliar)
summary(fpath_lm)
fpath_cwm <- as.data.frame(fpath_cwm)
resid_foliar_extra <- cbind(resid_foliar_extra, fpath_cwm)

##FOLIAR SAPROBE CWM##
fsap_cwm <- lm(saprobe ~ climate + pH + foliar, data = merge_foliar)
summary(fsap_cwm)
fsap_cwm <- resid(fsap_cwm)
fsap_lm <- lm(fsap_cwm ~ species, data = merge_foliar)
summary(fsap_lm)
fsap_cwm_sp <- as.data.frame(fsap_cwm)
resid_foliar_extra <- cbind(resid_foliar_extra, fsap_cwm_sp)

##FOLIAR SAPROBE CWM##
fsap_cwm <- lm(saprobe ~ climate + pH + species, data = merge_foliar)
summary(fsap_cwm)
fsap_cwm_nut <- resid(fsap_cwm)
fsap_lm <- aov(fsap_cwm_nut ~ foliar, data = merge_foliar)
summary(fsap_lm)
fsap_cwm_nut <- as.data.frame(fsap_cwm_nut)
resid_foliar_extra <- cbind(resid_foliar_extra, fsap_cwm_nut)

##FOLIAR PATHOGEN SPECIFICITY##
fpath_spec <- lm(path.host.spe ~ climate + pH + foliar, data = merge_foliar)
summary(fpath_spec)
fpath_spec <- resid(fpath_spec)
fpath_lm <- lm(fpath_spec ~ species, data = merge_foliar)
summary(fpath_lm)
fpath_spec <- as.data.frame(fpath_spec)
resid_foliar_extra <- cbind(resid_foliar_extra, fpath_spec)

##FOLIAR SAPROBE SPECIFICITY##
fsap_spec <- lm(sap.host.spe ~ climate + pH + foliar, data = merge_foliar)
summary(fsap_spec)
fsap_spec <- resid(fsap_spec)
fsap_lm <- lm(fsap_spec ~ species, data = merge_foliar)
summary(fsap_lm)
fsap_spec <- as.data.frame(fsap_spec)
resid_foliar_extra <- cbind(resid_foliar_extra, fsap_spec)

##SAVE FILES FOR FIGURES##
saveRDS(resid_foliar_extra, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\resid_foliar_extra.rds")
saveRDS(resid_soil_extra, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\resid_soil_extra.rds")


