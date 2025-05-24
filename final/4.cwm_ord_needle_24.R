##read in packages
library(tidyverse)
library(vegan)
library(ecole)
library(phyloseq)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)

##read in new_otu table, env1, and trait matrix
new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_needle_24.rds")
env1 <-  readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_needle_24.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_needle_24.rds")
dist <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\dist_needle_24.rds")
taxa <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\taxa_needle_24.rds")
per.table <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.needle_24.rds")

##ordinate
nmds <- metaMDS(dist)
nmds
m0 <- MDSrotate(nmds, as.numeric(per.table$elevation))
ecole:::fitstats_nms(m0)
m0
plot(m0, "sites")
env_sub <- subset(per.table, select =  c(2:4))
e <- envfit(m0, env_sub, na.rm = T)
e
plot(e)

##Create a df in ggplot with NMS scores and environmental variables
m_ggplot <- as.data.frame(vegan::scores(m0, display = c("sites")))

##GGplot matrix with environmental scores
env.scores <- as.data.frame(scores(e, display = "vectors"))
head(env.scores)
env.scores <- cbind(env.scores, env.variables = rownames(env.scores))
env.scores <- cbind(env.scores, pval = e$vectors$pvals)
env.scores

##isolate only the significant environmental scores
env.scores <- subset(env.scores, pval <= .1)
head(env.scores)

##GGplot base plot
nms_ggplot <- ggplot(m_ggplot, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(NMDS1, NMDS2, fill = per.table$species, size = per.table$elevation, stroke =2), color = "black", pch = 21) + scale_size_continuous(range = c(2, 15)) + scale_alpha(guide="none") + stat_ellipse(geom = "polygon", aes(fill = per.table$species), alpha = 0.25) + coord_equal()
nms_ggplot
new_plot <- (nms_ggplot + scale_fill_brewer(palette = "Dark2") + theme_classic(base_size=20) + labs( fill = "tree host species", size = "site elevation (m)"))
new_plot
new_plot_1 <- (new_plot + geom_segment(data = env.scores, aes(x = 0, xend=NMDS1*1.5, y=0, yend=NMDS2*1.5), arrow = arrow(length = unit(0.6, "cm")), colour = "grey10", lwd=1.5) + ggrepel::geom_label_repel(data = env.scores, aes(x=NMDS1, y=NMDS2, label = env.variables), size = 6, cex = 4, segment.size = 0.5, nudge_x = 1))
new_plot_1
needle_plot <- (new_plot_1 + guides(fill = guide_legend(override.aes = list(size = 5))) + theme(legend.text=element_text(size=20), legend.title=element_text(size=20)))
needle_plot
soil_plot
ggarrange(soil_plot, needle_plot, ncol=2, labels = c("C", "D"), align = "hv", font.label = list(size = 30), legend = "right", common.legend = TRUE)

##perform cwm##
cwm <- makecwm(new_otu, tra, na.rm = TRUE)
cwm <- subset(cwm, select = -c(8:9, 11:12,18:19))
head(cwm)

##combine per.table and cwm
cwm.per <- merge(cwm, per.table, by = 0)
cwm.per <- column_to_rownames(cwm.per, "Row.names")

##Heat Map
##heatmap based on elevation
heat <- subset(cwm.per, select = c(2:7,9:13,17:25,29))
heat_1 <- (heat %>% 
             group_by(elevation) %>%
             summarize_all(~ mean(.x, na.rm = TRUE)))
rownames(heat_1) <- heat_1$elevation
heat_needle <- heat_1
heat_needle <- as.data.frame(heat_needle)
heat_needle <- (heat_needle %>% map_df(rev))
rownames(heat_needle) <- heat_needle$elevation
heat_needle <- as.matrix(heat_needle)
heat_needle <- subset(heat_needle, select = -c(1))
heat_needle_final <- ComplexHeatmap::pheatmap(heat_needle, scale = "column", treeheight_row = 0, treeheight_col = 0, col = brewer.pal(8, "Greens"), cluster_cols = FALSE, cluster_rows = FALSE, border_color = NA, legend = FALSE, labels_col = c("foliar pathogens", "foliar saprobes"), fontsize = 25)
heat_needle_final <- ComplexHeatmap::pheatmap(heat_needle, scale = "column", treeheight_row = 0, treeheight_col = 0, col = brewer.pal(8, "Greens"), cluster_cols = FALSE, cluster_rows = FALSE, border_color = NA, legend = FALSE, fontsize = 25)
heat_needle_final
heat_soil_final + heat_needle_final 

##Look at dataset richness and diversity
env1 <- env1[match(row.names(new_otu), row.names(env1)),]
per.table <- env1[match(row.names(new_otu), row.names(per.table)),]

x <- specnumber(new_otu)
y <- diversity(new_otu)

per.table$total.fungal.diversity <- y
per.table$richness <- x
env1$total.fungal.diversity <- y
env1$richness <- x
needle_d <- subset(per.table, select = c(1,4,5,9,10))
needle_d
saveRDS(needle_d, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needle_richness_24.rds")

##look at diversity and tree host
anova <- aov(total.fungal.diversity) ~ elevation, data = per.table)
summary(anova)

lm <- lm(pathogen ~ elevation, data = cwm.per)
summary(lm)
plot1 <- ggplot(cwm.per, aes(x = elevation, y = pathogen)) + geom_point()
plot1
plot2 <- ggplot(env1, aes(x = elevation, y = total.fungal.diversity)) + geom_point()
grid.arrange(plot1,plot2)

