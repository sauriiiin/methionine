##### DELETION SCREEN ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 07/13/2021 

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(ggforce)
library(plotly)
library(scales)
library(reshape2)
library(locfit)
library(growthcurver)
library(rstatix)
library(gtools)
library(locfit)
library(growthrates)
library(RMariaDB)
library(genefilter)
library(apeglm)
library(clusterProfiler)
library(org.Sc.sgd.db)

out_path <- "~/R/Projects/methionine/data"
fig_path <- "~/R/Projects/methionine/figures"
res_path <- "~/R/Projects/methionine/results"

expt.name <- "deletion"

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### GATHER DATA
stages <- c('Pre-Screen #1','Pre-Screen #2','Final Screen')
arms <- c('SD-Met-Cys+Gal','SD+Met-Cys+Gal')
phenotypes <- c('Beneficial','Neutral','Deleterious')

p2c <- dbGetQuery(conn, 'select a.*, b.orf_name from MET_DEL_pos2coor a, MET_DEL_pos2orf_name b
                  where a.pos = b.pos
                  order by density, plate, col, row')

info <- data.frame(rbind(c('PS1','MM',1536,'MET_DEL'), c('PS1','PM',1536,'MET_DEL'),
                         c('PS2','MM',1536,'MET_DEL'), c('PS2','PM',1536,'MET_DEL'),
                         c('FS','MM',6144,'MET_DEL_FS_MM'), c('FS','PM',6144,'MET_DEL')))
colnames(info) <- c('stage','arm','density','p2c')

orfs <- dbGetQuery(conn, 'select distinct(orf_name) from MET_DEL_FS_MM_6144_FITNESS_STATS
                   where hours = 22 and N >= 2
                   and orf_name in
                   (select distinct(orf_name) from MET_DEL_PS2_MM_1536_FITNESS_STATS
                   where hours = 68 and N >= 2
                   and orf_name in
                   (select distinct(orf_name) from MET_DEL_PS1_MM_1536_FITNESS_STATS
                   where hours = 69 and N >= 2))')

pgs <- dbGetQuery(conn, 'select orf_name from PROTOGENES where pg_2012 = 1')

data <- NULL
data.sum <- NULL
for (i in seq(1,dim(info)[1])) {
  temp <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate, b.row, b.col
                                   from MET_DEL_%s_%s_%s_FITNESS a,
                                   %s_pos2coor b
                                   where a.pos = b.pos
                                   order by a.hours, b.plate, b.col, b.row',
                                   info[i,1],info[i,2],info[i,3],
                                   info[i,4]))
  temp <- temp[!is.na(temp$orf_name) & temp$orf_name != 'NULL',]
  
  temp2 <- dbGetQuery(conn, sprintf('select a.*, b.p
                                    from MET_DEL_%s_%s_%s_FITNESS_STATS a
                                    left join
                                    MET_DEL_%s_%s_%s_PVALUE b
                                    on a.hours = b.hours and a.strain_id = b.strain_id
                                    order by a.hours',
                                    info[i,1],info[i,2],info[i,3],
                                    info[i,1],info[i,2],info[i,3]))

  # for (h in unique(temp$hours)) {
  #   for (o in unique(temp$orf_name)) {
  #     temp$average[temp$orf_name == o & temp$hours == h][isoutlier(temp$average[temp$orf_name == o & temp$hours == h], 2) |
  #                                                          isoutlier(temp$fitness[temp$orf_name == o & temp$hours == h], 2)] <- NA
  #     temp$fitness[temp$orf_name == o & temp$hours == h][isoutlier(temp$average[temp$orf_name == o & temp$hours == h], 2) |
  #                                                          isoutlier(temp$fitness[temp$orf_name == o & temp$hours == h], 2)] <- NA
  #   }
  # }
  temp$stage <- info[i,1]
  temp$arm <- info[i,2]
  
  temp2$stage <- info[i,1]
  temp2$arm <- info[i,2]
  
  data <- rbind(data,temp)
  data.sum <- rbind(data.sum,temp2)
}
data <- data.frame(data)
data$orf_type[data$orf_name %in% pgs$orf_name] <- 'Proto-gene'
data$orf_type[is.na(data$orf_type)] <- 'Gene'
data$arm <- as.character(data$arm)
data$stage <- as.character(data$stage)
data$arm[data$arm == 'MM'] <- 'SD-Met-Cys+Gal'
data$arm[data$arm == 'PM'] <- 'SD+Met-Cys+Gal'
data$stage[data$stage == 'PS1'] <- 'Pre-Screen #1'
data$stage[data$stage == 'PS2'] <- 'Pre-Screen #2'
data$stage[data$stage == 'FS'] <- 'Final Screen'
data$rep <- as.numeric(str_trunc(as.character(data$pos), 5, side = 'left', ellipsis = ''))

data$stage <- factor(data$stage, levels = stages)
data$arm <- factor(data$arm, levels = arms)
head(data)

data.sum <- data.frame(data.sum)
data.sum$orf_type[data.sum$orf_name %in% pgs$orf_name] <- 'Proto-gene'
data.sum$orf_type[is.na(data.sum$orf_type)] <- 'Gene'
data.sum$arm <- as.character(data.sum$arm)
data.sum$stage <- as.character(data.sum$stage)
data.sum$arm[data.sum$arm == 'MM'] <- 'SD-Met-Cys+Gal'
data.sum$arm[data.sum$arm == 'PM'] <- 'SD+Met-Cys+Gal'
data.sum$stage[data.sum$stage == 'PS1'] <- 'Pre-Screen #1'
data.sum$stage[data.sum$stage == 'PS2'] <- 'Pre-Screen #2'
data.sum$stage[data.sum$stage == 'FS'] <- 'Final Screen'
data.sum$phenotype[data.sum$p <= 0.05 & data.sum$cs_mean > 1] <- 'Beneficial'
data.sum$phenotype[data.sum$p <= 0.05 & data.sum$cs_mean < 1] <- 'Deleterious'
data.sum$phenotype[is.na(data.sum$phenotype)] <- 'Neutral'

data.sum$stage <- factor(data.sum$stage, levels = stages)
data.sum$arm <- factor(data.sum$arm, levels = arms)
data.sum$phenotype <- factor(data.sum$phenotype, levels = phenotypes)

for (a in arms) {
  for (s in stages) {
    data$saturation[data$arm == a & data$stage == s] <- 
      max(data$hours[data$arm == a & data$stage == s])
  }
}
head(data)

for (a in arms) {
  for (s in stages) {
    data.sum$saturation[data.sum$arm == a & data.sum$stage == s] <- 
      max(data.sum$hours[data.sum$arm == a & data.sum$stage == s])
  }
}
head(data.sum)

##### SUMMARIZE RESULTS
data.mad <- data %>%
  group_by(arm, stage, hours, rep) %>%
  summarize(fitness.median = median(fitness, na.rm = T), cs.median = median(average, na.rm = T),
            fitness.mad = mad(fitness, na.rm = T), cs.mad = mad(average, na.rm = T),
            .groups = 'keep') %>%
  data.frame()
data <- merge(data, data.mad, by = c('stage','arm','hours','rep'))
data$fitness[data$fitness < (data$fitness.median - 2*data$fitness.mad) |
               data$fitness > (data$fitness.median + 2*data$fitness.mad)] <- NA
data$average[data$average < (data$cs.median - 2*data$cs.mad) |
               data$average > (data$cs.median + 2*data$cs.mad)] <- NA

temp <- data[!(data$arm == 'MM' & data$stage == 'PS2' & data$plate %in% c(1,2,3,4,5,6,10,13)),] %>%
  # filter(fitness <= 50) %>%
  group_by(arm, stage, hours, orf_name) %>%
  summarise(fitness = median(fitness, na.rm = T), cs = median(average, na.rm = T), .groups = 'keep') %>%
  data.frame()
data.sum <- merge(data.sum, temp, by = c('stage','arm','hours','orf_name'))

data.sum <- merge(data.sum, data.sum %>%
                    group_by(stage, arm, orf_name) %>%
                    summarise(var = max(cs_std), .groups = 'keep') %>%
                    data.frame(), by = c('stage','arm','orf_name'), all = T)

data.lim <- data[!(data$arm == 'MM' & data$stage == 'PS2' & data$plate %in% c(1,2,3,4,5,6,10,13)),] %>%
  filter(orf_name == 'BY4741') %>%
  group_by(arm, stage, hours, orf_name) %>%
  summarize(cs_ll = quantile(average, 0.025, na.rm = T),
            cs_m = median(average, na.rm = T),
            cs_ul = quantile(average, 0.975, na.rm = T),
            fitness_ll = quantile(fitness, 0.025, na.rm = T),
            fitness_m = median(fitness, na.rm = T),
            fitness_ul = quantile(fitness, 0.975, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

for (a in arms) {
  for (s in stages) {
    data.lim$saturation[data.lim$arm == a & data.lim$stage == s] <- 
      max(data.lim$hours[data.lim$arm == a & data.lim$stage == s])
  }
}

data.cont.lim <- merge(data.lim[data.lim$arm == 'SD-Met-Cys+Gal' & data.lim$hours == data.lim$saturation,c(-3,-4,-11)],
                       data.lim[data.lim$arm == 'SD+Met-Cys+Gal' & data.lim$hours == data.lim$saturation,c(-3,-4,-11)],
                   by = c('stage'), suffixes = c('_MM','_PM'), all = T)
save(data.cont.lim, file = 'data/final/211123_del_cont_lims.RData')

###### PLATEMAPS 
plot.platemap <- p2c %>%
  filter(density > 384) %>%
  ggplot(aes(x = col, y = row)) +
  geom_tile(aes(fill = orf_name), col = 'black') +
  scale_fill_discrete(guide = F, na.value = 'white') +
  labs(x = 'Columns', y = 'Rows') +
  scale_y_reverse() +
  facet_wrap(~density*plate, scales = 'free', ncol = 10) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        panel.grid = element_blank(),
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
ggsave(sprintf("%s/%s/PLATEMAP.jpg",fig_path, expt.name), plot.platemap,
       height = one.5c*1.5, width = two.c*3, units = 'mm',
       dpi = 600)

##### PLOT GROWTH CURVES
plot.cs.gc <- data.sum %>%
  filter(orf_name %in% orfs$orf_name) %>%
  ggplot(aes(x = hours, y = cs)) +
  geom_line(aes(group = orf_name, col = orf_name), lwd = 0.7, alpha = 0.8) +
  # stat_summary(aes(group = orf_name, col = orf_name), fun=mean, geom="line", lwd =0.7) +
  scale_color_discrete(guide = F) +
  geom_line(data = data.lim, aes(x = hours, y = cs_ll), col = 'red', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = cs_m), col = 'black', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = cs_ul), col = 'red', linetype = 'dashed', lwd = 0.5) +
  labs(x = 'Time(hours)',
       y = 'Colony Size (pixels)') +
  facet_wrap(.~arm*stage, scales = 'free') +
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
ggsave(sprintf("%s/%s/ALL_GROWTH_CURVES.jpg",fig_path, expt.name), plot.cs.gc,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)

##### FIGURING OUT THE SUDDEN DROP
## ANS: plate 17 wasn't analyzed properly at later time points

# hi <- merge(data[data$hours == 50,], data[data$hours == 54,], 
#       by = c('arm','stage','orf_name','pos','density','plate','row','col'))
# hi[hi$average.x > hi$average.y & !is.na(hi$average.x) & hi$orf_name != 'BY4741' &
#      !(hi$plate %in% c(1,2,3,4,5,6,10,13)),]
# 
# hi$diff <- hi$average.x - hi$average.y
# hi[hi$diff > 500 & !is.na(hi$diff) & hi$orf_name != 'BY4741' &
#      !(hi$plate %in% c(1,2,3,4,5,6,10,13)),]
# 
# hello <- data[data$arm == 'MM' & data$stage == 'PS2' & data$orf_name %in% hi$orf_name[hi$cs.x > hi$cs.y] &
#                     data$hours == 54,]

###### SAP PATHWAY MEMBERS
orfs.sap <- read.csv(file = 'data/deletion/SAP_ORFs.csv')
orfs.sap <- orfs.sap[orfs.sap$standard_name %in% c('MET3','MET14','MET16','MET5','MET10'),]
data.sap <- merge(data.sum, orfs.sap, by = 'orf_name')

plot.sap.gc <- data.sap %>%
  ggplot(aes(x = hours, y = fitness)) +
  geom_line(aes(group = standard_name, col = standard_name), lwd = 0.7, alpha = 0.8) +
  geom_point(data = data.sap[data.sap$hours == data.sap$saturation,],
             aes(x = hours, y = fitness, fill = standard_name), size = 1.5, col = 'black', shape = 21) + 
  geom_text_repel(data = data.sap[data.sap$hours == data.sap$saturation,],
                  aes(x = hours, y = fitness, label = standard_name, col = standard_name), size = 1.2,
                  force = 2, max.overlaps = 30) +
  geom_line(data = data.lim, aes(x = hours, y = fitness_ll), col = 'red', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = fitness_m), col = 'black', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = fitness_ul), col = 'red', linetype = 'dashed', lwd = 0.5) +
  scale_color_discrete(guide = F) + scale_fill_discrete(guide = F) +
  labs(x = 'Time (hours)', y = 'Fitness') +
  # facet_zoom(ylim = c(0, 5), zoom.data = ifelse(fitness <= 5, NA, FALSE)) +
  facet_wrap(.~arm*stage) +
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
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(xlim = c(0,80),
                  ylim = c(0,2))
ggsave(sprintf("%s/%s/SAP_GROWTH_CURVES.jpg",fig_path, expt.name), plot.sap.gc,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)

plot.sap.cs.gc <- data.sap %>%
  ggplot(aes(x = hours, y = cs)) +
  geom_line(aes(group = standard_name, col = standard_name), lwd = 0.7, alpha = 0.8) +
  geom_point(data = data.sap[data.sap$hours == data.sap$saturation,],
             aes(x = hours, y = cs, fill = standard_name), size = 1.5, col = 'black', shape = 21) + 
  geom_text_repel(data = data.sap[data.sap$hours == data.sap$saturation,],
                  aes(x = hours, y = cs, label = standard_name, col = standard_name), size = 1.2,
                  force = 2, max.overlaps = 30) +
  geom_line(data = data.lim, aes(x = hours, y = cs_ll), col = 'red', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = cs_m), col = 'black', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = cs_ul), col = 'red', linetype = 'dashed', lwd = 0.5) +
  scale_color_discrete(guide = F) + scale_fill_discrete(guide = F) +
  labs(x = 'Time (hours)', y = 'Fitness') +
  # facet_zoom(ylim = c(0, 5), zoom.data = ifelse(fitness <= 5, NA, FALSE)) +
  facet_wrap(.~arm*stage) +
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
ggsave(sprintf("%s/%s/SAP_CS_GROWTH_CURVES.jpg",fig_path, expt.name), plot.sap.cs.gc,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)

##### JAKES MUTANTS
orfs.jm <- data.frame(orf_name = c('YLR303W','YNL277W','YJR010W','YER091C','YGL125W','YGR155W','YPL023C','YGL184C'),
                      standard_name = c('MET15','MET2','MET3','MET6','MET13','CYS4','MET12','STR3'))
data.jm <- merge(data.sum, orfs.jm, by = 'orf_name')

plot.jm.gc <- data.jm %>%
  ggplot(aes(x = hours, y = fitness)) +
  geom_line(aes(group = standard_name, col = standard_name), lwd = 0.7, alpha = 0.8) +
  geom_point(data = data.jm[data.jm$hours == data.jm$saturation,],
             aes(x = hours, y = fitness, fill = standard_name), size = 1.5, col = 'black', shape = 21) + 
  geom_text_repel(data = data.jm[data.jm$hours == data.jm$saturation,],
                  aes(x = hours, y = fitness, label = standard_name, col = standard_name), size = 1.2,
                  force = 2, max.overlaps = 30) +
  geom_line(data = data.lim, aes(x = hours, y = fitness_ll), col = 'red', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = fitness_m), col = 'black', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = fitness_ul), col = 'red', linetype = 'dashed', lwd = 0.5) +
  scale_color_discrete(guide = F) + scale_fill_discrete(guide = F) +
  labs(x = 'Time (hours)', y = 'Fitness') +
  # facet_zoom(ylim = c(0, 5), zoom.data = ifelse(fitness <= 5, NA, FALSE)) +
  facet_wrap(.~arm*stage) +
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
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(xlim = c(0,80),
                  ylim = c(0,2))
ggsave(sprintf("%s/%s/JM_GROWTH_CURVES.jpg",fig_path, expt.name), plot.jm.gc,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)


##### ALL SULFUR METABOLIC PROCESS GENES
orfs.sul <- read.table(file = 'data/deletion/sulfur_metabolic_process.tsv', sep = '\t', header = T)

data.sul <- merge(data.sum, data.frame(rbind(cbind(GENENAME = 'NULL', bitr(orfs.sul$SYMBOL,
                                                                           fromType = "ORF", toType = c("ORF","DESCRIPTION"),
                                                                           OrgDb = org.Sc.sgd.db)),
                                             bitr(orfs.sul$SYMBOL,
                                                  fromType = "GENENAME", toType = c("ORF","DESCRIPTION"),
                                                  OrgDb = org.Sc.sgd.db)), row.names = NULL),
                  by.x = 'orf_name', by.y = 'ORF')
data.sul$GENENAME <- as.character(data.sul$GENENAME)
data.sul$GENENAME[data.sul$GENENAME == 'NULL'] <- data.sul$orf_name[data.sul$GENENAME == 'NULL']
head(data.sul)

plot.sul.gc <- data.sul %>%
  # filter(var < 10) %>%
  ggplot(aes(x = hours, y = fitness)) +
  geom_line(aes(group = GENENAME, col = GENENAME), lwd = 0.7, alpha = 0.8) +
  # geom_smooth(method = 'loess', aes(group = GENENAME, col = GENENAME), lwd = 0.7, alpha = 0.8, se = F) +
  geom_point(data = data.sul[data.sul$hours == data.sul$saturation,],
             aes(x = hours, y = fitness, fill = GENENAME), size = 1.5, col = 'black', shape = 21) + 
  geom_text_repel(data = data.sul[data.sul$hours == data.sul$saturation,],
                  aes(x = hours, y = fitness, label = GENENAME, col = GENENAME), size = 1.2,
                  force = 2, max.overlaps = 30) +
  geom_line(data = data.lim, aes(x = hours, y = fitness_ll), col = 'red', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = fitness_m), col = 'black', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = fitness_ul), col = 'red', linetype = 'dashed', lwd = 0.5) +
  scale_color_discrete(guide = F) + scale_fill_discrete(guide = F) +
  labs(x = 'Time (hours)', y = 'Fitness') +
  facet_wrap(.~arm*stage) +
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
ggsave(sprintf("%s/%s/SUL_GROWTH_CURVES.jpg",fig_path, expt.name), plot.sul.gc,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)

###### DENSITY PLOTS
plot.fit.den <- data[data$hours == data$saturation,] %>%
  filter(orf_name == 'BY4741', fitness.mad <= 20) %>%
  ggplot(aes(x = 'Reference', y = fitness, fill = 'Reference')) +
  geom_violin(draw_quantiles = 0.5) +
  geom_violin(data = data[data$hours == data$saturation,] %>%
                filter(orf_name != 'BY4741', fitness.mad <= 20),
              aes(x = 'Mutants', y = fitness, fill = 'Mutants'), draw_quantiles = 0.5) +
  scale_y_log10() +
  # ggridges::geom_density_ridges(quantile_lines = F,
  #                     aes(y = 'Reference', fill = 'Reference'),
  #                     scale = 0.5, alpha = 0.7, size = 0.2,
  #                     vline_size = 0.2, vline_color = "black") +
  # ggridges::geom_density_ridges(data = data[data$hours == data$saturation,] %>%
  #             filter(orf_name != 'BY4741', fitness.mad <= 20), aes(x = fitness, y = 'Mutants', fill = 'Mutants'),
  #             scale = 0.5, alpha = 0.7, size = 0.2,
  #             vline_size = 0.2, vline_color = "black") +
  labs(y = 'Fitness (log10)') +
  scale_x_discrete(limits = c('Reference','Mutants')) +
  scale_fill_discrete(guide = F) +
  facet_wrap(.~arm * stage) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0.01,100)) #+
  # annotation_logticks(sides="b")
ggsave(sprintf("%s/%s/FITNESS_DENSITY.jpg",fig_path, expt.name), plot.fit.den,
       height = one.c, width = one.5c, units = 'mm',
       dpi = 600)

###### DIFFERENCE BETWEEN -MET AND +MET CONDITIONS
data.diff <- merge(data.sum[data.sum$arm == 'SD-Met-Cys+Gal' & data.sum$hours == data.sum$saturation,],
                   data.sum[data.sum$arm == 'SD+Met-Cys+Gal' & data.sum$hours == data.sum$saturation,],
                   by = c('stage','orf_name','strain_id','orf_type'), suffixes = c('_MM','_PM'), all = T)
data.diff$phenotype_MM <- as.character(data.diff$phenotype_MM)
data.diff$phenotype_MM[is.na(data.diff$phenotype_MM)] <- 'Dead'
data.diff <- data.diff %>%
  filter(strain_id %in% unique(data.diff$strain_id[data.diff$phenotype_MM != 'Dead' & data.diff$stage == 'Pre-Screen #1']))
# data.diff <- data.diff[data.diff$phenotype_PM[!is.na(data.diff$phenotype_PM)],]

data.diff$fitness_MM[data.diff$phenotype_MM == 'Dead'] <- 0
data.diff$cs_MM[data.diff$phenotype_MM == 'Dead'] <- 0

data.diff$fitness_diff <- data.diff$fitness_MM - data.diff$fitness_PM
data.diff$cs_diff <- data.diff$cs_MM - data.diff$cs_PM 
data.diff[data.diff$phenotype_MM == 'Dead',]

##### REFERENCE DIFF FOR EMPIRICAL DIFF DISTRIBUTION
diff.dist <- NULL
for (s in unique(data$stage)) {
  temp1 <- data$fitness[data$stage == s & data$arm == 'SD-Met-Cys+Gal' & data$orf_name == 'BY4741'
                        & data$hours == max(data$hours[data$stage == s & data$arm == 'SD-Met-Cys+Gal'])]
  temp2 <- data$fitness[data$stage == s & data$arm == 'SD+Met-Cys+Gal' & data$orf_name == 'BY4741'
                        & data$hours == max(data$hours[data$stage == s & data$arm == 'SD+Met-Cys+Gal'])]
  temp3 <- NULL
  for (i in 1:50000) {
    temp3 <- c(temp3, median(sample(temp1[!is.na(temp1)], 1)) - median(sample(temp2[!is.na(temp2)], 1)))
  }
  diff.dist <- rbind(diff.dist, cbind(stage = s, ll = quantile(temp3, 0.025)[[1]], m = quantile(temp3, 0.5)[[1]], ul = quantile(temp3, 0.975)[[1]]))
}
diff.dist <- data.frame(diff.dist)
diff.dist$ll <- as.numeric(as.character(diff.dist$ll))
diff.dist$m <- as.numeric(as.character(diff.dist$m))
diff.dist$ul <- as.numeric(as.character(diff.dist$ul))

data.diff <- merge(data.diff, diff.dist, by = 'stage')

data.diff$phenotype[data.diff$fitness_diff >= data.diff$ul] <- 'Beneficial'
data.diff$phenotype[data.diff$fitness_diff <= data.diff$ll] <- 'Deleterious'
# data.diff$phenotype[data.diff$phenotype_MM == 'Dead'] <- 'Dead'
# data.diff$phenotype[data.diff$phenotype_MM == 'Dead'] <- 'Deleterious'
data.diff$phenotype[is.na(data.diff$phenotype)] <- 'Neutral'

##### FITNESS PLOT FOR DIFFERENCE
plot.fit.diff <- data.diff %>%
  ggplot(aes(x = fitness_diff, y = '')) +
  ggridges::geom_density_ridges(scale = 1.8, alpha = 0.7, size = 0.2,
                                vline_size = 0.2, vline_color = "black") +
  geom_vline(data = diff.dist, aes(xintercept = ul)) +
  geom_vline(data = diff.dist, aes(xintercept = m)) +
  geom_vline(data = diff.dist, aes(xintercept = ll)) +
  labs(x = 'Fitness Differential (log10)', y = 'Density') +
  scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
  facet_wrap(.~stage, nrow = 1) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        # axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(xlim = c(-2,10))
ggsave(sprintf("%s/%s/FITNESS_DIFF_DENSITY.jpg",fig_path, expt.name), plot.fit.diff,
       height = 50, width = one.5c, units = 'mm',
       dpi = 600)

##### GO/KEGG ENRICHMENT
allgenes <- unique(data.diff$orf_name)
allgenes <- bitr(allgenes, fromType = "ORF",
                 toType = c("ENTREZID","GENENAME","ENSEMBL"),
                 OrgDb = org.Sc.sgd.db)
allgenes <- allgenes[!is.na(allgenes$ENSEMBL),]

goe <- data.frame()
kegg <- data.frame()
for (s in levels(data.diff$stage)) {
  for (p in unique(data.diff$phenotype[data.diff$stage == s])) {
    temp.deg <- data.diff$orf_name[data.diff$stage == s & data.diff$phenotype == p]
    temp.deg <- bitr(temp.deg, fromType = "ORF",
                     toType = c("ENTREZID","GENENAME","ENSEMBL"),
                     OrgDb = org.Sc.sgd.db)
    temp.deg <- temp.deg[!is.na(temp.deg$ENSEMBL),]
    
    temp.goe <- enrichGO(gene          = temp.deg$ENSEMBL,
                         universe      = allgenes$ENSEMBL,
                         OrgDb         = org.Sc.sgd.db,
                         keyType       = "ENSEMBL",
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)
    if (dim(temp.goe)[1] == 0) {
      cat(sprintf('There are no GO term enrichment for %s ORFs in %s.\n',p,s))
    } else{
      cat(sprintf('GO term enrichment for %s ORFs in %s are:\n%s\n',
                  p,s,
                  paste(temp.goe$Description,collapse = ', ')))
      goe <- rbind(goe, data.frame(temp.goe, stage = s, phenotype = p))
    }
    
    temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                            universe     = allgenes$ENSEMBL,
                            organism     = 'sce',
                            pvalueCutoff = 0.05)
    if (dim(temp.kegg)[1] == 0) {
      cat(sprintf('There are no KEGG pathway enriched for %s ORFs in %s.\n',p,s))
    } else{
      cat(sprintf('KEGG pathway enrichment for %s ORFs in %s are:\n%s\n',
                  p,s,
                  paste(temp.kegg$Description,collapse = ', ')))
      kegg <- rbind(kegg, data.frame(temp.kegg, stage = s, phenotype = p))
    }
  }
}
goe$GeneRatio <- as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,2])
goe$BgRatio <- as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,2])
goe$GO <- paste0(goe$ONTOLOGY, '_', goe$Description)

kegg$GeneRatio <- as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,2])
kegg$BgRatio <- as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,2])

goe <- goe[order(goe$stage,goe$phenotype,-goe$GeneRatio,-goe$Count,goe$qvalue),]
kegg <- kegg[order(kegg$stage,kegg$phenotype,-kegg$GeneRatio,-kegg$Count,kegg$qvalue),]

write.csv(goe, file = sprintf('%s/%s/go_enrichments.csv', res_path, expt.name))
write.csv(kegg, file = sprintf('%s/%s/kegg_enrichments.csv', res_path, expt.name))

##### GO-KEGG PLOTS
head(goe)
goe$phenotype <- factor(goe$phenotype, levels = c('Beneficial','Neutral','Deleterious','Dead'))
goe$stage <- factor(goe$stage, levels = c('Pre-Screen #1','Pre-Screen #2','Final Screen'))
kegg$phenotype <- factor(kegg$phenotype, levels = c('Beneficial','Neutral','Deleterious','Dead'))
kegg$stage <- factor(kegg$stage, levels = c('Pre-Screen #1','Pre-Screen #2','Final Screen'))

##### DEATH ANALYSIS
data.ded <- merge(merge(data.sum[data.sum$stage == 'Pre-Screen #1' & data.sum$hours == data.sum$saturation, c('arm','orf_name','strain_id','fitness','cs','orf_type','phenotype')],
                        data.sum[data.sum$stage == 'Pre-Screen #2' & data.sum$hours == data.sum$saturation, c('arm','orf_name','strain_id','fitness','cs','orf_type','phenotype')],
                        by = c('arm','orf_name','strain_id','orf_type'), suffixes = c('_PS1','_PS2'), all = T),
                  data.sum[data.sum$stage == 'Final Screen' & data.sum$hours == data.sum$saturation, c('arm','orf_name','strain_id','fitness','cs','orf_type','phenotype')],
                  by = c('arm','orf_name','strain_id','orf_type'), suffixes = c('','_FS'), all = T)
data.cnts$phenotype <- factor(data.cnts$phenotype, levels = c('Beneficial','Neutral','Deleterious','Dead'))
data.ded$phenotype_PS2 <- as.character(data.ded$phenotype_PS2)
data.ded$phenotype <- as.character(data.ded$phenotype)
data.ded$phenotype_PS2[is.na(data.ded$phenotype_PS2)] <- 'Dead'
data.ded$phenotype[is.na(data.ded$phenotype)] <- 'Dead'
data.ded$phenotype_PS2 <- factor(data.ded$phenotype_PS2, levels = c('Beneficial','Neutral','Deleterious','Dead'))
data.ded$phenotype <- factor(data.ded$phenotype, levels = c('Beneficial','Neutral','Deleterious','Dead'))

data.ded$fitness_PS2[data.ded$phenotype_PS2 == 'Dead'] <- 0
data.ded$cs_PS2[data.ded$phenotype_PS2 == 'Dead'] <- 0
data.ded$fitness[data.ded$phenotype == 'Dead'] <- 0
data.ded$cs[data.ded$phenotype == 'Dead'] <- 0

colnames(data.ded) <- c(colnames(data.ded)[1:10],'fitness_FS','cs_FS','phenotype_FS')

data.ded <- merge(data.ded[data.ded$arm == 'SD-Met-Cys+Gal',-1], data.ded[data.ded$arm == 'SD+Met-Cys+Gal',-1],
                  by = c('orf_name','strain_id','orf_type'), suffixes = c('_MM','_PM'), all = T)

data.ded$diff_PS1 <- data.ded$fitness_PS1_MM - data.ded$fitness_PS1_PM
data.ded$diff_PS2 <- data.ded$fitness_PS2_MM - data.ded$fitness_PS2_PM
data.ded$diff_FS <- data.ded$fitness_FS_MM - data.ded$fitness_FS_PM
head(data.ded)

data.ded %>%
  group_by(orf_type, phenotype_PS1_MM) %>%
  dplyr::count()

##### COUNTS OF PGS BENEFICIAL, DELETERIOUS, NEUTRAL AND DEAD
head(data.diff)
data.ded <- merge(merge(merge(data.ded, data.diff[,c('stage','orf_name','strain_id','orf_type','phenotype')] %>% filter(stage == 'Pre-Screen #1'),
                  by = c('orf_name','strain_id','orf_type'), all = T),
            data.diff[,c('stage','orf_name','strain_id','orf_type','phenotype')] %>% filter(stage == 'Pre-Screen #2'),
            by = c('orf_name','strain_id','orf_type'), all = T),
      data.diff[,c('stage','orf_name','strain_id','orf_type','phenotype')] %>% filter(stage == 'Final Screen'),
      by = c('orf_name','strain_id','orf_type'), all = T)
data.ded <- data.ded[,c(-29,-27,-25)]
colnames(data.ded) <- c(colnames(data.ded)[1:24],'phenotype_PS1_diff','phenotype_PS2_diff','phenotype_FS_diff')

data.cnts <- data.ded[,str_detect(colnames(data.ded), 'phenotype') | colnames(data.ded) %in% c('orf_name','orf_type')]
data.cnts <- melt(data.cnts, id.vars = c('orf_name','orf_type'), variable.name = 'level', value.name = 'phenotype')
data.cnts$phenotype[data.cnts$phenotype == 'Dead'] <- 'Deleterious'

data.cnts <- merge(data.cnts[!is.na(data.cnts$phenotype),] %>%
                     group_by(level,orf_type,phenotype) %>%
                     dplyr::count(), data.cnts[!is.na(data.cnts$phenotype),] %>%
                     group_by(level,orf_type) %>%
                     dplyr::count(), by = c('level','orf_type'), suffixes = c('','_total'))
data.cnts$percentage <- data.cnts$n/data.cnts$n_total * 100 

data.cnts$stage[str_detect(data.cnts$level,'PS1')] <- 'Pre-Screen #1'
data.cnts$stage[str_detect(data.cnts$level,'PS2')] <- 'Pre-Screen #2'
data.cnts$stage[str_detect(data.cnts$level,'FS')] <- 'Final Screen'
data.cnts$stage <- factor(data.cnts$stage, levels = c('Pre-Screen #1','Pre-Screen #2','Final Screen'))

data.cnts$arm[str_detect(data.cnts$level,'MM')] <- 'SD-Met-Cys+Gal'
data.cnts$arm[str_detect(data.cnts$level,'PM')] <- 'SD+Met-Cys+Gal'
data.cnts$arm[str_detect(data.cnts$level,'diff')] <- 'Differential'
data.cnts$arm <- factor(data.cnts$arm, levels = c('SD+Met-Cys+Gal','SD-Met-Cys+Gal','Differential'))

data.cnts$phenotype <- factor(data.cnts$phenotype, levels = c('Beneficial','Neutral','Deleterious','Dead'))

plot.pheno.cnts <- data.cnts %>%
  ggplot(aes(x = stage, y = percentage/100)) +
  geom_point(aes(col = phenotype)) +
  geom_line(data = data.cnts %>% filter(orf_type == 'Gene'),
                                        aes(col = phenotype, group = phenotype, linetype = orf_type)) +
  geom_line(data = data.cnts %>% filter(orf_type == 'Proto-gene'),
            aes(col = phenotype, group = phenotype, linetype = orf_type)) +
  geom_label_repel(data = data.cnts %>% filter(stage == 'Final Screen'), 
                  aes(label = n, fill = phenotype), size = 1.5, xlim = c('Final Screen', NA)) +
  scale_y_continuous(breaks = seq(-10,200,10)/100, minor_breaks = seq(-10,200,2)/100) +
  labs(x = 'Screen Stage', y = 'ORF Proportion with Phenotype') +
  scale_linetype_discrete(name = 'ORF Type') +
  scale_color_manual(name = 'Phenotype',
                     values = c('Beneficial' = '#4CAF50','Neutral'= '#757575','Deleterious' = '#FF5252')) +
  scale_fill_manual(guide = F,
                    values = c('Beneficial' = '#4CAF50','Neutral'= '#757575','Deleterious' = '#FF5252')) +
  facet_wrap(.~arm) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,.7))
ggsave(sprintf("%s/%s/PHENOTYPE_COUNTS.jpg",fig_path, expt.name), plot.pheno.cnts,
       height = one.c, width = two.c, units = 'mm',
       dpi = 600)  

write.csv(data.cnts, file = sprintf('%s/%s/phenotype_counts.csv',res_path,expt.name))

##### ODDS RATIOS OF PGS BENEFICIAL, DELETERIOUS, NEUTRAL AND DEAD
head(data.cnts)

data.or <- NULL
for (a in unique(data.cnts$arm)) {
  for (s in unique(data.cnts$stage[data.cnts$arm == a])) {
    for (p in data.cnts$phenotype[data.cnts$arm == a & data.cnts$stage == s]) {
      n_pg_effect <- data.cnts$n[data.cnts$arm == a & data.cnts$stage == s & data.cnts$phenotype == p & data.cnts$orf_type == 'Proto-gene']
      n_pg_noeffect <- data.cnts$n_total[data.cnts$arm == a & data.cnts$stage == s & data.cnts$phenotype == p & data.cnts$orf_type == 'Proto-gene'] -
        data.cnts$n[data.cnts$arm == a & data.cnts$stage == s & data.cnts$phenotype == p & data.cnts$orf_type == 'Proto-gene']
      n_g_effect <- data.cnts$n[data.cnts$arm == a & data.cnts$stage == s & data.cnts$phenotype == p & data.cnts$orf_type == 'Gene']
      n_g_noeffect <- data.cnts$n_total[data.cnts$arm == a & data.cnts$stage == s & data.cnts$phenotype == p & data.cnts$orf_type == 'Gene'] -
        data.cnts$n[data.cnts$arm == a & data.cnts$stage == s & data.cnts$phenotype == p & data.cnts$orf_type == 'Gene']
      
      ftest <- fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))
      data.or <- rbind(data.or, c(a,s,p,n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect,
                                  ftest$conf.int[1],ftest$estimate[[1]],ftest$conf.int[2],ftest$p.value))
    }
  }
}
data.or <- data.frame(data.or, stringsAsFactors = F)
colnames(data.or) <- c('arm','stage','phenotype','pg_w_effect','pg_wo_effect','g_w_effect','g_wo_effect','top','or','bottom','p')
head(data.or)
data.or$significant[data.or$p <= 0.05] <- 'Yes'
data.or$significant[is.na(data.or$significant)] <- 'No'

data.or$phenotype <- factor(data.or$phenotype, levels = c('Beneficial','Neutral','Deleterious','Dead'))
data.or$stage <- factor(data.or$stage, levels = c('Pre-Screen #1','Pre-Screen #2','Final Screen'))
data.or$arm <- factor(data.or$arm, levels = c('SD+Met-Cys+Gal','SD-Met-Cys+Gal','Differential'))

data.or$pg_w_effect <- as.numeric(data.or$pg_w_effect)
data.or$pg_wo_effect <- as.numeric(data.or$pg_wo_effect)
data.or$g_w_effect <- as.numeric(data.or$g_w_effect)
data.or$g_wo_effect <- as.numeric(data.or$g_wo_effect)
data.or$top <- as.numeric(data.or$top)
data.or$or <- as.numeric(data.or$or)
data.or$bottom <- as.numeric(data.or$bottom)
data.or$p <- as.numeric(data.or$p)

data.or$label[data.or$p > 0.05] <- 'ns'
data.or$label[data.or$p <= 0.05] <- '*'
data.or$label[data.or$p <= 0.01] <- '**'
data.or$label[data.or$p <= 0.001] <- '***'
data.or$label[data.or$p <= 0.0001] <- '****'

plot.or <- ggplot(data.or,aes(x=phenotype,y=or,ymin=bottom,ymax=top))+
  geom_point(stat="identity",  shape=21, size=2, stroke=1, fill = "white")+
  geom_errorbar(width=0.25)+
  geom_hline(yintercept=1,linetype="dashed",color="red")+
  geom_text(aes(x = phenotype, y = 2, label = label), size = 2, col = 'red') +
  labs(x='Phenotype',y='Odds Ratio') +
  # scale_y_continuous(trans="log10", breaks=c(0.2,0.5,1,2,5,10)) + 
  # annotation_logticks(sides="l") +
  facet_wrap(.~stage*arm) +
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
ggsave(sprintf("%s/%s/PHENOTYPE_ODDS.jpg",fig_path, expt.name), plot.or,
       height = one.5c, width = one.5c, units = 'mm',
       dpi = 600) 

write.csv(data.or, file = sprintf('%s/%s/phenotype_odds.csv',res_path,expt.name))

##### PHENOTYPE OVERLAP BETWEEN -MET AND +MET CONDITIONS
data.olp <- data.ded[,str_detect(colnames(data.ded), 'phenotype') | colnames(data.ded) %in% c('orf_name','orf_type')]
data.olp[data.olp == 'Dead' & !is.na(data.olp)] <- 'Deleterious'

data.olp$PS1 <- paste(data.olp$phenotype_PS1_MM, data.olp$phenotype_PS1_PM, sep = '/')
data.olp$PS2 <- paste(data.olp$phenotype_PS2_MM, data.olp$phenotype_PS2_PM, sep = '/')
data.olp$FS <- paste(data.olp$phenotype_FS_MM, data.olp$phenotype_FS_PM, sep = '/')

data.olp <- data.olp[,c(1,2,12,13,14)]
data.olp <- melt(data.olp, id.vars = c('orf_name','orf_type'), variable.name = 'stage', value.name = 'phenotype')
data.olp <- data.olp[!str_detect(data.olp$phenotype, 'NA'),]

data.olp %>%
  group_by(stage, phenotype) %>%
  dplyr::count() %>%
  data.frame()

allgenes <- unique(data.olp$orf_name)
allgenes <- bitr(allgenes, fromType = "ORF",
                 toType = c("ENTREZID","GENENAME","ENSEMBL"),
                 OrgDb = org.Sc.sgd.db)
# allgenes <- allgenes[!is.na(allgenes$ENSEMBL),]

goe.olp <- data.frame()
kegg.olp <- data.frame()
for (s in unique(data.olp$stage)) {
  for (p in unique(data.olp$phenotype[data.olp$stage == s])) {
    temp.deg <- data.olp$orf_name[data.olp$stage == s & data.olp$phenotype == p]
    temp.deg <- bitr(temp.deg, fromType = "ORF",
                     toType = c("ENTREZID","GENENAME","ENSEMBL"),
                     OrgDb = org.Sc.sgd.db)
    temp.deg <- temp.deg[!is.na(temp.deg$ENSEMBL),]
    
    temp.goe <- enrichGO(gene          = temp.deg$ENSEMBL,
                         universe      = allgenes$ENSEMBL,
                         OrgDb         = org.Sc.sgd.db,
                         keyType       = "ENSEMBL",
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)
    if (dim(temp.goe)[1] == 0) {
      cat(sprintf('There are no GO term enrichment for %s ORFs in %s.\n',p,s))
    } else{
      cat(sprintf('GO term enrichment for %s ORFs in %s are:\n%s\n',
                  p,s,
                  paste(temp.goe$Description,collapse = ', ')))
      goe.olp <- rbind(goe.olp, data.frame(temp.goe, stage = s, phenotype = p))
    }
    
    temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                            universe     = allgenes$ENSEMBL,
                            organism     = 'sce',
                            pvalueCutoff = 0.05)
    if (dim(temp.kegg)[1] == 0) {
      cat(sprintf('There are no KEGG pathway enriched for %s ORFs in %s.\n',p,s))
    } else{
      cat(sprintf('KEGG pathway enrichment for %s ORFs in %s are:\n%s\n',
                  p,s,
                  paste(temp.kegg$Description,collapse = ', ')))
      kegg.olp <- rbind(kegg.olp, data.frame(temp.kegg, stage = s, phenotype = p))
    }
  }
}
goe.olp$GeneRatio <- as.numeric(str_split(goe.olp$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe.olp$GeneRatio,'/',simplify = T)[,2])
goe.olp$BgRatio <- as.numeric(str_split(goe.olp$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe.olp$BgRatio,'/',simplify = T)[,2])
goe.olp$GO <- paste0(goe.olp$ONTOLOGY, '_', goe.olp$Description)

kegg.olp$GeneRatio <- as.numeric(str_split(kegg.olp$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg.olp$GeneRatio,'/',simplify = T)[,2])
kegg.olp$BgRatio <- as.numeric(str_split(kegg.olp$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg.olp$BgRatio,'/',simplify = T)[,2])

goe.olp <- goe.olp[order(goe.olp$stage,goe.olp$phenotype,-goe.olp$GeneRatio,-goe.olp$Count,goe.olp$qvalue),]
kegg.olp <- kegg.olp[order(kegg.olp$stage,kegg.olp$phenotype,-kegg.olp$GeneRatio,-kegg.olp$Count,kegg.olp$qvalue),]

write.csv(goe.olp, file = sprintf('%s/%s/go_overlap_enrichments.csv', res_path, expt.name))
write.csv(kegg.olp, file = sprintf('%s/%s/kegg_overlap_enrichments.csv', res_path, expt.name))

##### SAVING ALL DATA
save(data, data.cnts, data.ded, data.diff, data.jm, data.lim, data.mad,
     data.or, data.sap, data.sul, data.sum, diff.dist, goe, kegg,
     orfs, orfs.jm, orfs.sap, orfs.sul, pgs,
     file = sprintf('%s/%s/all_data.RData',res_path,expt.name))
