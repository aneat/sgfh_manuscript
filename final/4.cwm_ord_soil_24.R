  ##read in packages
  library(tidyverse)
  library(vegan)
  library(ecole)
  library(RColorBrewer)
  library(ComplexHeatmap)
  
  ##read in new_otu table, env1, and trait matrix
  new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_soil_24.rds")
  env1 <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_soil_24.rds")
  tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_soil_24.rds")
  dist <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\dist_soil_24.rds")
  per.table.soil <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.soil_24.rds")
  
  ##ordinate
  nmds <- metaMDS(dist)
  m0_soil <- MDSrotate(nmds, as.numeric(per.table.soil$elevation))
  plot(m0_soil, "sites")
  m0_soil
  plot(m0_soil, "sites")
  
  ##cumulative variance explained is metric fit, R2m##
  # metric fit = squared correlation between original dissimilarities and ordination distances. Foreign to NMS.  Null: no linear relationship between dissimilarities and ordination distances.
  ecole:::fitstats_nms(m0_soil)
  
  env_sub <- subset(per.table.soil, select =  c(2,3,4))
  es <- envfit(m0_soil, env_sub, na.rm = T)
  plot(es)
  es
  ##Create a df in ggplot with NMS scores and environmental variables
  m_ggplot_soil <- as.data.frame(vegan::scores(m0_soil, display = "sites"))
  
  ##GGplot matrix with environmental scores
  env.scores.soil <- as.data.frame(vegan::scores(es, display = "vectors"))
  head(env.scores.soil)
  env.scores.soil <- cbind(env.scores.soil, env.variables = rownames(env.scores.soil))
  env.scores.soil <- cbind(env.scores.soil, pval = es$vectors$pvals)
  env.scores.soil
  
  ##isolate only the significant environmental scores
  env.scores.soil <- subset(env.scores.soil, env.variables == "elevation")
  env.scores.soil
  
  ##GGplot base plot
  nms_ggplot_soil <- ggplot(m_ggplot_soil, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(NMDS1, NMDS2, fill = per.table.soil$species, size = per.table.soil$elevation, stroke =2), color = "black", pch = 21) +  scale_size_continuous(range = c(2, 15)) + scale_alpha(guide="none") + stat_ellipse(geom = "polygon", aes(fill = per.table.soil$species), alpha = 0.25) + coord_equal()
                                                                                                                                                                                                                                                                                                                    
  nms_ggplot_soil
  new_plot_soil <- (nms_ggplot_soil + scale_fill_brewer(palette = "Dark2") + theme_classic(base_size=20) + labs( fill = "tree host species", size = "site elevation (m)"))
  new_plot_soil
  new_plot_1_soil <- (new_plot_soil + geom_segment(data = env.scores.soil, aes(x = 0, xend=NMDS1*1.5, y=0, yend=NMDS2*1.5), arrow = arrow(length = unit(0.6, "cm")), colour = "grey10", lwd=1.5) + ggrepel::geom_label_repel(data = env.scores.soil, aes(x=NMDS1*.4, y=NMDS2*.4, label = env.variables), size = 6, cex = 4, segment.size = 0.5))
  soil_plot <- (new_plot_1_soil + guides(fill = guide_legend(override.aes = list(size = 5))) + theme(legend.text=element_text(size=20), legend.title=element_text(size=20)))
  soil_plot
  ##soil_plot
  saveRDS(soil_plot,"\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soil_plot_24.rds")
  soil_plot <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soil_plot_24.rds")
  ##soil_plot
  ##needle_plot
  
  ##perform cwm
  cwm <- makecwm(new_otu, tra, na.rm = TRUE)
  cwm <- subset(cwm, select = -c(9, 11:12,18:19))
  
  ##combine per.table and cwm
  cwm.per <- merge(cwm, per.table.soil, by = 0)
  cwm.per <- column_to_rownames(cwm.per, "Row.names")
  saveRDS(cwm.per,"\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\cwm.per.soil_24.rds")
  
  
  ##heatmap based on elevation
  heat <- subset(cwm.per, select = c(1:14,17:26, 29))
  heat_1 <- (heat %>% 
                  group_by(elevation) %>%
                  summarize_all(~ mean(.x, na.rm = TRUE)))
  heat_soil <- heat_1
  heat_soil <- as.data.frame(heat_soil)
  heat_soil <- (heat_soil %>% map_df(rev))
  rownames(heat_soil) <- heat_soil$elevation
  heat_soil <- as.matrix(heat_soil)
  heat_soil <- subset(heat_soil, select = -c(1))
  heat_soil_final <- ComplexHeatmap::pheatmap(heat_soil, scale = "column", treeheight_row = 0, treeheight_col = 0, col = brewer.pal(8, "YlOrBr"), cluster_cols = FALSE, cluster_rows = FALSE, border_color = NA, legend = FALSE, labels_col = c("soil pathogens", "saprobes", "am", "ecm"), fontsize = 20)
  heat_soil_final <- ComplexHeatmap::pheatmap(heat_soil, scale = "column", treeheight_row = 0, treeheight_col = 0, col = brewer.pal(8, "YlOrBr"), cluster_cols = FALSE, cluster_rows = FALSE, border_color = NA, legend = FALSE, fontsize = 20)
  
  heat_soil_final
  
  ##Look at dataset richness and diversity
  env1 <- env1[match(row.names(new_otu), row.names(env1)),]
  per.table.soil <- env1[match(row.names(new_otu), row.names(per.table.soil)),]
  x <- specnumber(new_otu)
  y <- diversity(new_otu)
  per.table.soil$total.fungal.diversity <- y
  per.table.soil$richness <- x
  soil_d <- subset(per.table.soil, select = c(4,6,10))
  saveRDS(soil_d, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soil_richness_24.rds")
  
  ##look at diversity and tree host
  anova <- aov(total.fungal.diversity ~ elevation, data = per.table.soil)
  summary(anova)
  
  lm <- lm(pathogen ~ climate, data = cwm.per)
  summary(lm)
  plot1 <- ggplot(cwm.per, aes(x = climate, y = pathogen)) + geom_point()
  plot1
  plot2 <- ggplot(env1, aes(x = elevation, y = total.fungal.diversity)) + geom_point()
  ##grid.arrange(plot1,plot2)
  
