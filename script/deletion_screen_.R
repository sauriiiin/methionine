
load(sprintf("%s/%s/210728_environment.RData",out_path,expt.name))

for (a in unique(data.lim$arm)) {
  for (s in unique(data.lim$stage[data.lim$arm == a])) {
    data.lim$saturation[data.lim$arm == a & data.lim$stage == s] <- 
      max(data.lim$hours[data.lim$arm == a & data.lim$stage == s])
  }
}
head(data.lim)

diff.dist2 <- NULL
for (s in unique(diff.dist$stage)) {
  diff.dist2 <- rbind(diff.dist2, data.frame(stage = s,
                                        x = seq(0,100,0.01),
                                        y.ll = seq(0,100,0.01) - diff.dist$ll[diff.dist$stage == s],
                                        y.m = seq(0,100,0.01) - diff.dist$m[diff.dist$stage == s],
                                        y.ul = seq(0,100,0.01) - diff.dist$ul[diff.dist$stage == s]))
}
diff.dist2 <- data.frame(diff.dist2)
head(diff.dist2)

##### OVERLAP ANALYSIS
plot.olp <- data.diff %>%
  filter(stage == 'Final Screen') %>%
  ggplot(aes(x = fitness_MM, y = fitness_PM)) +
  geom_point(aes(col = phenotype), size = 0.8) +
  geom_hline(data = data.lim %>% filter(hours == saturation, arm == 'SD+Met-Cys+Gal', stage == 'Final Screen'),
             aes(yintercept = fitness_ll), linetype = 'dashed', col = 'red') +
  geom_hline(data = data.lim %>% filter(hours == saturation, arm == 'SD+Met-Cys+Gal', stage == 'Final Screen'),
             aes(yintercept = fitness_ul), linetype = 'dashed', col = 'red') +
  geom_vline(data = data.lim %>% filter(hours == saturation, arm == 'SD-Met-Cys+Gal', stage == 'Final Screen'),
             aes(xintercept = fitness_ll), linetype = 'dashed', col = 'red') +
  geom_vline(data = data.lim %>% filter(hours == saturation, arm == 'SD-Met-Cys+Gal', stage == 'Final Screen'),
             aes(xintercept = fitness_ul), linetype = 'dashed', col = 'red') +
  geom_line(data = diff.dist2 %>% filter(stage == 'Final Screen'), aes(x = x , y = y.ul)) +
  # geom_line(data = diff.dist2, aes(x = x , y = y.m), linetype = 'dashed') +
  geom_line(data = diff.dist2 %>% filter(stage == 'Final Screen'), aes(x = x , y = y.ll)) +
  geom_point(data = data.diff %>% filter(orf_name %in% unique(data.sul$orf_name), stage == 'Final Screen'),
             aes(x = fitness_MM, y = fitness_PM), shape = 0) +
  scale_x_log10(breaks = breaks_log(n = 10, base = 10),
                minor_breaks = NULL) +
  labs(x = 'Fitness in -Met',
      y = 'Fitness in +Met') +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Beneficial','Neutral','Deleterious','Dead'),
                     values = c('Beneficial' = '#4CAF50',
                                'Neutral'= '#757575',
                                'Deleterious' = '#FF5252',
                                'Dead' = 'black'),
                     na.translate = F) +
  coord_cartesian(xlim = c(0.01, 20),
                  ylim = c(0, 1.5)) +
  facet_wrap(.~stage) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/%s/FITNESS_OVERLAP_FS.jpg",fig_path, expt.name), plot.olp,
       height = one.5c, width = one.5c, units = 'mm',
       dpi = 600) 

##### OVERLAP INCLUDING DEAD
data.olp2 <- data.ded[,str_detect(colnames(data.ded), 'phenotype') | colnames(data.ded) %in% c('orf_name','orf_type')]

data.olp2$PS1 <- paste(data.olp2$phenotype_PS1_MM, data.olp2$phenotype_PS1_PM, sep = '/')
data.olp2$PS2 <- paste(data.olp2$phenotype_PS2_MM, data.olp2$phenotype_PS2_PM, sep = '/')
data.olp2$FS <- paste(data.olp2$phenotype_FS_MM, data.olp2$phenotype_FS_PM, sep = '/')

data.olp2 <- data.olp2[,c(1,2,12,13,14)]
data.olp2 <- melt(data.olp2, id.vars = c('orf_name','orf_type'), variable.name = 'stage', value.name = 'phenotype')
data.olp2 <- data.olp2[!str_detect(data.olp2$phenotype, 'NA'),]

head(data.olp2)
data.olp2 %>%
  group_by(stage, phenotype) %>%
  dplyr::count() %>%
  data.frame()

data.ded[!(data.ded$phenotype_FS_MM %in% c('Dead','Deleterious')) &
         (data.ded$phenotype_FS_PM %in% c('Dead','Deleterious')),] %>%
  dplyr::count()

data.diff %>%
  filter(
    # orf_name %in% unique(data.sul$orf_name),
    # phenotype_MM %in% c('Dead','Deleterious'),
    # phenotype_PM %in% c('Beneficial','Neutral'),
    # phenotype %in% c('Dead','Deleterious'),
    stage == 'Final Screen'
         ) %>%
  dplyr::count(phenotype_MM)


##### EFFECT SIZE ANALYSIS
## GO/KEGG ENRICHMENT
allgenes <- unique(data.diff$orf_name)
allgenes <- bitr(allgenes, fromType = "ORF",
                 toType = c("ENTREZID","GENENAME","ENSEMBL"),
                 OrgDb = org.Sc.sgd.db)
allgenes <- allgenes[!is.na(allgenes$ENSEMBL),]

goe.es <- data.frame()
kegg.es <- data.frame()
for (s in levels(data.diff$stage)) {
  for (es in seq(0,2,0.2)) {
    for (d in c(-1,1)) {
      temp.deg <- data.diff$orf_name[data.diff$stage == s & data.diff$fitness_diff*d >= es & !is.na(data.diff$fitness_diff)]
      if (length(temp.deg) != 0) {
        temp.deg <- bitr(temp.deg, fromType = "ORF",
                         toType = c("ENTREZID","GENENAME","ENSEMBL"),
                         OrgDb = org.Sc.sgd.db)
        temp.deg <- temp.deg[!is.na(temp.deg$ENSEMBL),]
        
        temp.goe.es <- enrichGO(gene          = temp.deg$ENSEMBL,
                                universe      = allgenes$ENSEMBL,
                                OrgDb         = org.Sc.sgd.db,
                                keyType       = "ENSEMBL",
                                ont           = "ALL",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05)
        
        if (length(temp.goe.es)[1] == 0) {
          # cat(sprintf('There are no GO term enrichment for %s ORFs in %s.\n',p,s))
        } else {
          if (dim(temp.goe.es)[1] != 0) {
            # cat(sprintf('GO term enrichment for %s ORFs in %s are:\n%s\n',
            #             p,s,
            #             paste(temp.goe.es$Description,collapse = ', ')))
            goe.es <- rbind(goe.es, data.frame(temp.goe.es, stage = s, effect_size = es, direction = d))
          }
        }
        
        temp.kegg.es <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                                      universe     = allgenes$ENSEMBL,
                                      organism     = 'sce',
                                      pvalueCutoff = 0.05)
        if (length(temp.kegg.es)[1] == 0) {
          # cat(sprintf('There are no kegg.es pathway enriched for %s ORFs in %s.\n',p,s))
        } else {
          if (dim(temp.kegg.es)[1] != 0) {
            # cat(sprintf('kegg.es pathway enrichment for %s ORFs in %s are:\n%s\n',
            #             p,s,
            #             paste(temp.kegg.es$Description,collapse = ', ')))
            kegg.es <- rbind(kegg.es, data.frame(temp.kegg.es, stage = s, effect_size = es, direction = d))
          }
        }
      }
    }
  }
}
goe.es$GeneRatio <- as.numeric(str_split(goe.es$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe.es$GeneRatio,'/',simplify = T)[,2])
goe.es$BgRatio <- as.numeric(str_split(goe.es$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe.es$BgRatio,'/',simplify = T)[,2])
goe.es$GO <- paste0(goe.es$ONTOLOGY, '_', goe.es$Description)

kegg.es$GeneRatio <- as.numeric(str_split(kegg.es$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg.es$GeneRatio,'/',simplify = T)[,2])
kegg.es$BgRatio <- as.numeric(str_split(kegg.es$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg.es$BgRatio,'/',simplify = T)[,2])

goe.es <- goe.es[order(goe.es$stage,goe.es$direction,goe.es$effect_size,-goe.es$GeneRatio,-goe.es$Count,goe.es$qvalue),]
kegg.es <- kegg.es[order(kegg.es$stage,kegg.es$direction,kegg.es$effect_size,-kegg.es$GeneRatio,-kegg.es$Count,kegg.es$qvalue),]

write.csv(goe.es, file = sprintf('%s/%s/go_es_enrichments.csv', res_path, expt.name))
write.csv(kegg.es, file = sprintf('%s/%s/kegg_es_enrichments.csv', res_path, expt.name))
  
##### SUL GENES IN BENEFICIAL CATEGORY
data.diff2 <- merge(data.diff, bitr(data.diff$orf_name, fromType = "ORF",
     toType = c("GENENAME","DESCRIPTION"),
     OrgDb = org.Sc.sgd.db), by.x = 'orf_name', by.y = 'ORF', all = T)

data.diff2 %>%
  filter(orf_name %in% unique(data.sul$orf_name), phenotype == 'Beneficial',
         stage == 'Final Screen')
data.diff2$sulfur_metabolism[data.diff2$orf_name %in% unique(data.sul$orf_name)] <- 'Yes'
data.diff2$sulfur_metabolism[is.na(data.diff2$sulfur_metabolism)] <- 'No'
data.diff2$sap[data.diff2$orf_name %in% unique(data.sap$orf_name)] <- 'Yes'
data.diff2$sap[is.na(data.diff2$sap)] <- 'No'
data.diff2$jm[data.diff2$orf_name %in% unique(data.jm$orf_name)] <- 'Yes'
data.diff2$jm[is.na(data.diff2$jm)] <- 'No'

head(data.diff2)
data.diff2 <- data.diff2[order(data.diff2$stage, data.diff2$orf_name),]

write.csv(data.diff2, file = sprintf('%s/%s/differential_analysis.csv',res_path,expt.name))

##### TRIPLE OVERLAP OF DELETERIOUS
allgenes <- unique(data.diff$orf_name)
allgenes <- bitr(allgenes, fromType = "ORF",
                 toType = c("ENTREZID","GENENAME","ENSEMBL"),
                 OrgDb = org.Sc.sgd.db)
allgenes <- allgenes[!is.na(allgenes$ENSEMBL),]

goe.ddd <- data.frame()
kegg.ddd <- data.frame()

temp.deg <- unique(data.diff$orf_name[data.diff$phenotype_MM  %in% c('Deleterious','Dead') &
                                 data.diff$phenotype_PM  %in% c('Deleterious','Dead') &
                                 data.diff$phenotype %in% c('Deleterious','Dead') &
                                 data.diff$stage == 'Final Screen'])
temp.deg <- bitr(temp.deg, fromType = "ORF",
                 toType = c("ENTREZID","GENENAME","ENSEMBL"),
                 OrgDb = org.Sc.sgd.db)
temp.deg <- temp.deg[!is.na(temp.deg$ENSEMBL),]

goe.ddd <- enrichGO(gene          = temp.deg$ENSEMBL,
                     universe      = allgenes$ENSEMBL,
                     OrgDb         = org.Sc.sgd.db,
                     keyType       = "ENSEMBL",
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)
goe.ddd <- data.frame(goe.ddd)

kegg.ddd <- enrichKEGG(gene          = temp.deg$ENSEMBL,
                        universe     = allgenes$ENSEMBL,
                        organism     = 'sce',
                        pvalueCutoff = 0.05)
kegg.ddd <- data.frame(kegg.ddd)


