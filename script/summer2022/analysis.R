##### DATA ANALYSIS FOR THE SUMMER SCREEN ANALYSIS
##### Screens done by Aaron Zhang, Alexis Berger and Brandon Garcia
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 07/29/2021 

##### INITIALIZE
library(genefilter)
library(clusterProfiler)
library(org.Sc.sgd.db)

source('/home/sbp29/R/Projects/methionine/script/summer2022/data.R')
data.fit.clean$condition <- factor(data.fit.clean$condition, levels = cnd_limits)
data.stats$condition <- factor(data.stats$condition, levels = cnd_limits)
data.cnts.all$condition <- factor(data.cnts.all$condition, levels = cnd_limits)

##### REFERENCE COLONY SIZE DISTRIBUTION
plot.ref.dist <- data.fit.clean %>%
  filter(hours == saturation, orf_name == 'BF_control') %>%
  ggplot(aes(x = average, y = condition)) +
  geom_density_ridges(aes(height = ..density..), stat = 'density', trim = TRUE) +
  geom_text(data = sat.times, aes(x = 2000, y = condition,
                                  label = sprintf('t = %d hours',saturation)),
            size = 2, hjust = 0, nudge_y = 0.5) +
  labs(x = 'Colony Size (pixel counts)', y = 'Condition') +
  scale_x_continuous(trans = 'pseudo_log', breaks = c(10,500,2000,6000)) +
  facet_wrap(.~stage*bag) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.title.y = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(xlim = c(0,10000))
ggsave(sprintf("%s/REFERENCE_DISTRIBUTION.jpg",fig_path), plot.ref.dist,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600) 


##### FITNESS DISTRIBUTION 
plot.fit.dist <- data.stats %>%
  filter(hours == saturation, N > 4) %>%
  ggplot(aes(x = fitness_median, y = condition)) +
  geom_density_ridges(aes(height = ..density..), stat = 'density', trim = TRUE) +
  geom_text(data = sat.times, aes(x = 3, y = condition,
                                  label = sprintf('t = %d hours',saturation)),
            size = 2, hjust = 0, nudge_y = 0.5) +
  labs(x = 'Fitness', y = 'Condition') +
  scale_x_continuous(trans = 'pseudo_log') +
  facet_wrap(.~stage*bag) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.title.y = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/FITNESS_DISTRIBUTION.jpg",fig_path), plot.fit.dist,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600) 

##### PHENOTYPE COUNTS
plot.pheno.perc <- data.cnts.all %>%
  ggplot(aes(x = "", y = percentage, fill = phenotype)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",percentage)),
                   position = position_stack(vjust = 0.5),
                   colour ='white',
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#757575',
                               'Beneficial'='#FFC107')) +
  facet_grid(condition~stage*bag,
             labeller = labeller(condition = c('SC-Ura+Gal' = 'SC-Ura+Gal',
                                               'SD+Met-Cys-Ura+Gal' = 'SD+Met-Cys-Ura+Gal',
                                               'SD-Met-Cys-Ura+Gal'= 'SD-Met-Cys-Ura+Gal',
                                               'SD+Met-Cys-Ura+Gal+FeEDTA' = 'SD+Met-Cys-Ura+Gal\n+FeEDTA',
                                               'SD-Met-Cys-Ura+Gal+FeEDTA' = 'SD-Met-Cys-Ura+Gal\n+FeEDTA'))) +
  labs(title = 'All ORFs') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))
ggsave(sprintf("%s/PHENOTYPE_PERCENTAGES.jpg",fig_path), plot.pheno.perc,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600) 


plot.pheno.perc.trn <- data.cnts.trn %>%
  ggplot(aes(x = "", y = percentage, fill = phenotype)) +
  geom_bar(stat = 'identity') +
  coord_polar("y", start = 0) +
  geom_label_repel(aes(label = sprintf("%0.2f%%",percentage)),
                   position = position_stack(vjust = 0.5),
                   colour ='white',
                   label.size = 0.15,
                   show.legend = F) +
  scale_fill_manual(name = 'Phenotype',
                    breaks = c('Beneficial','Neutral','Deleterious'),
                    values = c('Deleterious'='#3F51B5',
                               'Neutral'='#757575',
                               'Beneficial'='#FFC107')) +
  facet_grid(condition~stage*bag,
             labeller = labeller(condition = c('SC-Ura+Gal' = 'SC-Ura+Gal',
                                               'SD+Met-Cys-Ura+Gal' = 'SD+Met-Cys-Ura+Gal',
                                               'SD-Met-Cys-Ura+Gal'= 'SD-Met-Cys-Ura+Gal',
                                               'SD+Met-Cys-Ura+Gal+FeEDTA' = 'SD+Met-Cys-Ura+Gal\n+FeEDTA',
                                               'SD-Met-Cys-Ura+Gal+FeEDTA' = 'SD-Met-Cys-Ura+Gal\n+FeEDTA'))) +
  labs(title = 'Transient ORFs') +
  theme_linedraw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.key.size = unit(3, "mm"),
        legend.position = 'bottom',
        strip.text = element_text(size = titles),
        legend.box.spacing = unit(0.5,"mm"))
ggsave(sprintf("%s/PHENOTYPE_PERCENTAGES_TRANSIENT_ORFS.jpg",fig_path), plot.pheno.perc.trn,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600) 


##### BENEFICIAL IN -Met
## Either MM or MMEDTA
orfs.minusmet <- unique(data.stats$orf_name[data.stats$hours == data.stats$saturation &
                                       data.stats$condition_id %in% c('MM','MMEDTA') &
                                       data.stats$stage == 'FS1']) 
orfs.plusmet <- unique(data.stats$orf_name[data.stats$hours == data.stats$saturation &
                                       data.stats$condition_id %in% c('PM','PMEDTA') &
                                       data.stats$stage == 'FS1']) 

orfs.mmedta.ben <- unique(data.stats$orf_name[data.stats$hours == data.stats$saturation &
                                            data.stats$condition_id %in% c('MMEDTA') &
                                            data.stats$phenotype == 'Beneficial' &
                                            data.stats$stage == 'FS1']) 
orfs.mm.ben <- unique(data.stats$orf_name[data.stats$hours == data.stats$saturation &
                                              data.stats$condition_id %in% c('MM') &
                                              data.stats$phenotype == 'Beneficial' &
                                              data.stats$stage == 'FS1']) 
orfs.minusmet.ben <- unique(c(orfs.mm.ben, orfs.mmedta.ben)) 


orfs.pmedta.ben <- unique(data.stats$orf_name[data.stats$hours == data.stats$saturation &
                                            data.stats$condition_id %in% c('PMEDTA') &
                                            data.stats$phenotype == 'Beneficial' &
                                            data.stats$stage == 'FS1']) 
orfs.pm.ben <- unique(data.stats$orf_name[data.stats$hours == data.stats$saturation &
                                              data.stats$condition_id %in% c('PM') &
                                              data.stats$phenotype == 'Beneficial' &
                                              data.stats$stage == 'FS1']) 
orfs.plusmet.ben <- unique(c(orfs.pm.ben, orfs.pmedta.ben)) 

orfs.minusmet.ben.unique <- orfs.minusmet.ben[orfs.minusmet.ben %notin% orfs.plusmet.ben]

# mm.all <- bitr(mm.all, fromType = "ORF",
#                toType = c("ENTREZID","GENENAME","ENSEMBL"),
#                OrgDb = org.Sc.sgd.db)
# mm.ben.unique <- bitr(mm.ben.unique, fromType = "ORF",
#                       toType = c("ENTREZID","GENENAME","ENSEMBL"),
#                       OrgDb = org.Sc.sgd.db)
# 
# mm.ben.goe <- enrichGO(gene          = mm.ben.unique$ENSEMBL,
#                        universe      = mm.all$ENSEMBL,
#                        OrgDb         = org.Sc.sgd.db,
#                        keyType       = "ENSEMBL",
#                        ont           = "ALL",
#                        pAdjustMethod = "BH",
#                        pvalueCutoff  = 0.01,
#                        qvalueCutoff  = 0.05) %>%
#   data.frame()
# 
# mm.ben.kegg <- enrichKEGG(gene         = mm.ben.unique$ENSEMBL,
#                           universe     = mm.all$ENSEMBL,
#                           organism     = 'sce',
#                           pvalueCutoff = 0.05) %>%
#   data.frame()
# 
# mm.ben.goe$GeneRatio <- as.numeric(str_split(mm.ben.goe$GeneRatio,'/',simplify = T)[,1])/
#   as.numeric(str_split(mm.ben.goe$GeneRatio,'/',simplify = T)[,2])
# mm.ben.goe$BgRatio <- as.numeric(str_split(mm.ben.goe$BgRatio,'/',simplify = T)[,1])/
#   as.numeric(str_split(mm.ben.goe$BgRatio,'/',simplify = T)[,2])
# mm.ben.goe$GO <- paste0(mm.ben.goe$ONTOLOGY, '_', mm.ben.goe$Description)
# 
# mm.ben.kegg$GeneRatio <- as.numeric(str_split(mm.ben.kegg$GeneRatio,'/',simplify = T)[,1])/
#   as.numeric(str_split(mm.ben.kegg$GeneRatio,'/',simplify = T)[,2])
# mm.ben.kegg$BgRatio <- as.numeric(str_split(mm.ben.kegg$BgRatio,'/',simplify = T)[,1])/
#   as.numeric(str_split(mm.ben.kegg$BgRatio,'/',simplify = T)[,2])
# 
# mm.ben.goe <- mm.ben.goe[order(-mm.ben.goe$GeneRatio,-mm.ben.goe$Count,mm.ben.goe$qvalue),]
# mm.ben.kegg <- mm.ben.kegg[order(-mm.ben.kegg$GeneRatio,-mm.ben.kegg$Count,mm.ben.kegg$qvalue),]


##### FITNESS PLOT

merge(data.sum, orf_cat_aaron[,c(1,7)], by = 'orf_name', all.x = T) %>%
  filter(hours == saturation, stage == 'FS1') %>%
  ggplot(aes(x = rep, y = fitness_median, color = orf_type)) +
  geom_point() +
  scale_y_continuous(trans = 'pseudo_log') +
  facet_grid(stage*condition~.) +
  coord_cartesian(ylim = c(0,5))



##### FITNESS SOM
library(kohonen)
library(lattice)
load('~/R/Projects/rnaseek/functions/som_functions.RData')

som.matrix <- fit.matrix[,str_detect(colnames(fit.matrix), 'FS1') & str_detect(colnames(fit.matrix), 'phenotype') & str_detect(colnames(fit.matrix), 'nobag')]
som.matrix[som.matrix == 'neutral' & !is.na(som.matrix)] <- 0
som.matrix[som.matrix == 'beneficial' & !is.na(som.matrix)] <- 1
som.matrix[som.matrix == 'deleterious' & !is.na(som.matrix)] <- -1
som.matrix <- sapply(som.matrix, as.numeric, simplify = 'list')
rownames(som.matrix) <- fit.matrix$strain_id

som_data <- run_soms_analysis( som.matrix, 3, 3, 7)


som_clusters <- som_data$pivot %>%
  group_by(cluster, id_a) %>%
  summarise(cluster = max(cluster), id_a = max(id_a)) %>%
  data.frame()
# som_clusters <- merge(som_clusters, bitr(som_clusters$id_a, fromType = "ORF",
#                                          toType = c("GENENAME","DESCRIPTION"),
#                                          OrgDb = org.Sc.sgd.db), by.x = 'id_a', by.y = 'ORF', all = T)


ggplot(som_data$pivot[som_data$pivot$cluster > 0,],
       aes(x = x,
           y = y)) +
  # geom_vline(xintercept = c(69,137), lwd = 0.2, linetype = 'dashed') +
  # geom_rect(mapping=aes(xmin=0, xmax=69, ymin=-20, ymax=25, fill= 'Pre-Screen #1', alpha= 'Pre-Screen #1')) +
  # geom_rect(mapping=aes(xmin=69, xmax=137, ymin=-20, ymax=25, fill= 'Pre-Screen #2', alpha= 'Pre-Screen #2')) +
  # geom_rect(mapping=aes(xmin=137, xmax=159, ymin=-20, ymax=25, fill= 'Final Screen', alpha= 'Final Screen')) +
  geom_line(aes(group = id_a), col = '#BDBDBD', alpha = 0.3) +
  # geom_line(data = som_data$pivot[som_data$pivot$cluster > 0 &
  #                                   som_data$pivot$id_a %in% unique(data.sul$orf_name),],
  #           aes(x = data.som2$hours[data.som2$arm == 'SD-Met-Cys+Gal' &
  #                                     data.som2$orf_name %in% unique(data.sul$orf_name)],
  #               y = y, group = id_a), col = 'blue', alpha = 0.3) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               aes(group = 1), geom="ribbon", alpha = 0.4) +
  stat_summary(fun=mean, aes(group=1), geom="line", colour="blue") +
  scale_y_continuous(trans = 'pseudo_log') +
  labs(title = 'Fitness Dynamics',
       x = 'Condition',
       y = 'Fitness') +
  facet_wrap(.~cluster, nrow = 4) +
  theme_linedraw() +
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        # axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position = 'bottom',
        legend.margin = margin(0.1,0.1,0.1,0.1, "mm"),
        strip.text = element_text(size = 5,
                                  margin = margin(0.1,0.1,0.1,0.1, "mm")))


