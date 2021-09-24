##### OVEREXPRESSION SCREEN ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 07/21/2021 

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
library(growthrates)
library(RMariaDB)
library(genefilter)
library(apeglm)
library(clusterProfiler)
library(org.Sc.sgd.db)

out_path <- "~/R/Projects/methionine/data"
fig_path <- "~/R/Projects/methionine/figures"
res_path <- "~/R/Projects/methionine/results"

expt.name <- "overexpression"

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### GATHER DATA
stages <- c('Pre-Screen #1','Pre-Screen #2','Final Screen')
arms <- c('SD-Met-Cys-Ura+Gal','SC-Ura+Gal')
phenotypes <- c('Beneficial','Neutral','Deleterious','Dead')
orfs <- c('Proto-gene','Gene','Reference')

p2c <- dbGetQuery(conn, 'select a.*, b.orf_name from MET_OE_pos2coor a, MET_OE_pos2orf_name b
                  where a.pos = b.pos
                  order by density, plate, col, row')

info <- data.frame(rbind(c('PS1','MM',1536,'MET_OE'), c('PS1','PM',1536,'MET_OE'),
                         c('PS2','MM',1536,'MET_OE'), c('PS2','PM',1536,'MET_OE'),
                         c('FS','PM',6144,'MET_OE'),c('FS','MM',6144,'MET_OE')))
colnames(info) <- c('stage','arm','density','p2c')

pgs <- dbGetQuery(conn, 'select orf_name from PROTOGENES where pg_2012 = 1')

data <- NULL
data.sum <- NULL
for (i in seq(1,dim(info)[1])) {
  temp <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate, b.row, b.col
                                   from MET_OE_%s_%s_%s_FITNESS a,
                                   %s_pos2coor b
                                   where a.pos = b.pos
                                   order by a.hours, b.plate, b.col, b.row',
                                   info[i,1],info[i,2],info[i,3],
                                   info[i,4]))
  temp <- temp[!is.na(temp$orf_name) & temp$orf_name != 'NULL',]
  
  temp2 <- dbGetQuery(conn, sprintf('select a.*, b.p
                                    from MET_OE_%s_%s_%s_FITNESS_STATS a
                                    left join
                                    MET_OE_%s_%s_%s_PVALUE b
                                    on a.hours = b.hours and a.strain_id = b.strain_id
                                    order by a.hours',
                                    info[i,1],info[i,2],info[i,3],
                                    info[i,1],info[i,2],info[i,3]))
  
  temp$stage <- info[i,1]
  temp$arm <- info[i,2]
  
  temp2$stage <- info[i,1]
  temp2$arm <- info[i,2]
  
  data <- rbind(data,temp)
  data.sum <- rbind(data.sum,temp2)
}
data <- data.frame(data)
data$orf_type[data$orf_name %in% pgs$orf_name] <- 'Proto-gene'
data$orf_type[data$orf_name == 'BF_control'] <- 'Reference'
data$orf_type[is.na(data$orf_type)] <- 'Gene'
data$arm <- as.character(data$arm)
data$stage <- as.character(data$stage)
data$arm[data$arm == 'MM'] <- 'SD-Met-Cys-Ura+Gal'
data$arm[data$arm == 'PM'] <- 'SC-Ura+Gal'
data$stage[data$stage == 'PS1'] <- 'Pre-Screen #1'
data$stage[data$stage == 'PS2'] <- 'Pre-Screen #2'
data$stage[data$stage == 'FS'] <- 'Final Screen'
data$rep <- as.numeric(str_trunc(as.character(data$pos), 5, side = 'left', ellipsis = ''))

data$stage <- factor(data$stage, levels = stages)
data$arm <- factor(data$arm, levels = arms)
data$orf_type <- factor(data$orf_type, levels = orfs)
head(data)

data.sum <- data.frame(data.sum)
data.sum$orf_type[data.sum$orf_name %in% pgs$orf_name] <- 'Proto-gene'
data.sum$orf_type[data.sum$orf_name == 'BF_control'] <- 'Reference'
data.sum$orf_type[is.na(data.sum$orf_type)] <- 'Gene'
data.sum$arm <- as.character(data.sum$arm)
data.sum$stage <- as.character(data.sum$stage)
data.sum$arm[data.sum$arm == 'MM'] <- 'SD-Met-Cys-Ura+Gal'
data.sum$arm[data.sum$arm == 'PM'] <- 'SC-Ura+Gal'
data.sum$stage[data.sum$stage == 'PS1'] <- 'Pre-Screen #1'
data.sum$stage[data.sum$stage == 'PS2'] <- 'Pre-Screen #2'
data.sum$stage[data.sum$stage == 'FS'] <- 'Final Screen'
data.sum$phenotype[data.sum$p <= 0.05 & data.sum$cs_mean > 1] <- 'Beneficial'
data.sum$phenotype[data.sum$p <= 0.05 & data.sum$cs_mean < 1] <- 'Deleterious'
data.sum$phenotype[is.na(data.sum$phenotype)] <- 'Neutral'

data.sum$stage <- factor(data.sum$stage, levels = stages)
data.sum$arm <- factor(data.sum$arm, levels = arms)
data.sum$phenotype <- factor(data.sum$phenotype, levels = phenotypes)
data.sum$orf_type <- factor(data.sum$orf_type, levels = orfs)

data <- data[data$hours != 180,]
data.sum <- data.sum[data.sum$hours != 180,]

for (a in unique(data$arm)) {
  for (s in unique(data$stage[data$arm == a])) {
    data$saturation[data$arm == a & data$stage == s] <- 
      max(data$hours[data$arm == a & data$stage == s])
  }
}
head(data)

for (a in unique(data.sum$arm)) {
  for (s in unique(data.sum$stage[data.sum$arm == a])) {
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

temp <- data %>%
  # filter(fitness <= 50) %>%
  group_by(arm, stage, hours, orf_name) %>%
  summarise(fitness = median(fitness, na.rm = T), cs = median(average, na.rm = T), .groups = 'keep') %>%
  data.frame()
data.sum <- merge(data.sum, temp, by = c('stage','arm','hours','orf_name'))

data.sum <- merge(data.sum, data.sum %>%
                    group_by(stage, arm, orf_name) %>%
                    summarise(var = max(cs_std), .groups = 'keep') %>%
                    data.frame(), by = c('stage','arm','orf_name'), all = T)

data.lim <- data %>%
  filter(orf_name == 'BF_control') %>%
  group_by(arm, stage, hours, orf_name, rep) %>%
  summarize(average = median(average, na.rm = T),
            fitness = median(fitness, na.rm = T),
            .groups = 'keep') %>%
  group_by(arm, stage, hours, orf_name) %>%
  summarize(cs_ll = quantile(average, 0.025, na.rm = T),
            cs_m = median(average, na.rm = T),
            cs_ul = quantile(average, 0.975, na.rm = T),
            fitness_ll = quantile(fitness, 0.025, na.rm = T),
            fitness_m = median(fitness, na.rm = T),
            fitness_ul = quantile(fitness, 0.975, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

data.sum <- merge(data.sum, data.lim[,-4], by = c('arm','stage','hours'))
data.sum$es <- abs(data.sum$fitness - data.sum$fitness_m)/data.sum$fitness_m

##### PLOT GROWTH CURVES
plot.cs.gc <- data.sum %>%
  ggplot(aes(x = hours, y = cs)) +
  geom_line(aes(group = orf_name, col = orf_name), lwd = 0.7, alpha = 0.8) +
  # stat_summary(aes(group = orf_name, col = orf_name), fun=mean, geom="line", lwd =0.7) +
  scale_color_discrete(guide = F) +
  geom_line(data = data.lim, aes(x = hours, y = cs_ll), col = 'red', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = cs_m), col = 'black', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = data.lim, aes(x = hours, y = cs_ul), col = 'red', linetype = 'dashed', lwd = 0.5) +
  labs(x = 'Time(hours)',
       y = 'Colony Size (pixels)') +
  facet_wrap(.~arm*stage, scales = 'free_x') +
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
ggsave(sprintf("%s/%s/ALL_GROWTH_CURVES_CS.jpg",fig_path, expt.name), plot.cs.gc,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)

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
  facet_wrap(.~arm*stage, scale = 'free_x') +
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
  facet_wrap(.~arm*stage, scale = 'free_x') +
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
  facet_wrap(.~arm*stage, scale = 'free_x') +
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
ggsave(sprintf("%s/%s/JM_GROWTH_CURVES.jpg",fig_path, expt.name), plot.jm.gc,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)

###### DENSITY PLOTS
plot.fit.vio <- data[data$hours == data$saturation,] %>%
  filter(orf_name == 'BF_control', fitness.mad <= 20) %>%
  ggplot(aes(x = 'Reference', y = fitness, fill = 'Reference')) +
  geom_violin(draw_quantiles = 0.5) +
  geom_violin(data = data[data$hours == data$saturation,] %>%
                filter(orf_name != 'BF_control', fitness.mad <= 20),
              aes(x = 'Mutants', y = fitness, fill = 'Mutants'), draw_quantiles = 0.5) +
  scale_y_log10() +
  # ggridges::geom_density_ridges(quantile_lines = F,
  #                     aes(y = 'Reference', fill = 'Reference'),
  #                     scale = 0.5, alpha = 0.7, size = 0.2,
  #                     vline_size = 0.2, vline_color = "black") +
  # ggridges::geom_density_ridges(data = data[data$hours == data$saturation,] %>%
  #             filter(orf_name != 'BF_control', fitness.mad <= 20), aes(x = fitness, y = 'Mutants', fill = 'Mutants'),
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
  coord_cartesian(ylim = c(0.01,100))
ggsave(sprintf("%s/%s/FITNESS_DENSITY_VIOLIN.jpg",fig_path, expt.name), plot.fit.vio,
       height = one.c, width = one.5c, units = 'mm',
       dpi = 600)

plot.fit.den <- data.sum %>%
  filter(hours == saturation, orf_type != 'Reference') %>%
  ggplot(aes(x = fitness, y = orf_type, fill = orf_type)) +
  geom_density_ridges(quantile_lines = TRUE,
                      scale = 2, alpha = 0.7, size = 0.2,
                      vline_size = 0.2, vline_color = "black") +
  facet_wrap(.~arm*stage, nrow = 3) +
  scale_fill_discrete(name = 'ORF Type') +
  labs(x = 'Fitness') +
  coord_cartesian(xlim = c(0,2)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
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
ggsave(sprintf("%s/%s/FITNESS_DENSITY.jpg",fig_path, expt.name), plot.fit.den,
       height = one.c, width = one.5c, units = 'mm',
       dpi = 600)

# ####### CRISPR TARGETS
# data.crsp2 <- merge(data.sum %>% filter(hours == saturation, orf_type == 'Proto-gene',
#                                         arm == 'SD-Met-Cys-Ura+Gal', phenotype == 'Beneficial'),
#       data.crsp, by = 'orf_name')
# 
# unique(hi$orf_name)[!(unique(hi$orf_name) %in% unique(data.crsp2$orf_name))]
# 
# length(unique(data.crsp2$orf_name[data.crsp2$Hsu2013 >= 95 & data.crsp2$Doench2014OnTarget >= 0.65 & !is.na(data.crsp2$Doench2014OnTarget)]))
# 
# plot.crspr.good <- ggplot(data = data.sum %>%
#          filter(orf_name %in%
#                   unique(data.crsp$orf_name[data.crsp$Hsu2013 >= 95 & data.crsp$Doench2014OnTarget >= 0.7 &
#                                               !is.na(data.crsp$Doench2014OnTarget)]),
#                 hours == saturation, orf_type == 'Proto-gene',
#                 arm == 'SD-Met-Cys-Ura+Gal'),
#        aes(y = cs_mean, x = 'Proto-genes')) +
#   geom_jitter(aes(col = phenotype), size = 1) + 
#   geom_violin(fill = 'transparent') +
#   geom_jitter(data = data.sum %>%
#                 filter(orf_name %in% unique(data.crsp$orf_name[data.crsp$Hsu2013 >= 95 & data.crsp$Doench2014OnTarget >= 0.7 &
#                                                                  !is.na(data.crsp$Doench2014OnTarget)]),
#                        hours == saturation, orf_type == 'Gene',
#                        arm == 'SD-Met-Cys-Ura+Gal', cs_mean <= 4),
#               aes(y = cs_mean, x = 'Gene', col = phenotype), size = 1) +
#   geom_violin(data = data.sum %>%
#                 filter(orf_name %in% unique(data.crsp$orf_name[data.crsp$Hsu2013 >= 95 & data.crsp$Doench2014OnTarget >= 0.7 &
#                                                                  !is.na(data.crsp$Doench2014OnTarget)]),
#                        hours == saturation, orf_type == 'Gene',
#                        arm == 'SD-Met-Cys-Ura+Gal', cs_mean <= 4),
#               aes(y = cs_mean, x = 'Gene'), fill = 'transparent') +
#   scale_y_continuous(minor_breaks = seq(-2,10,0.1)) +
#   facet_wrap(.~arm*stage, nrow = 3) +
#   scale_color_discrete(name = 'Phenotype') +
#   labs(title = 'ORFs with at least one good CRISPR Target',
#        x = 'ORF Type', y = 'Fitness') +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%s/%s/GOOD_CRISPR_FITNESS.jpg",fig_path, expt.name), plot.crspr.good,
#        height = one.c, width = one.c, units = 'mm',
#        dpi = 600)
# 
# # 107 PGs are Beneficial >> 24 PGs have no predicted CRISPR Target
# data.sum %>%
#   filter(orf_name %in% unique(data.crsp$orf_name[data.crsp$Hsu2013 >= 95 & data.crsp$Doench2014OnTarget >= 0.7 &
#                                                    !is.na(data.crsp$Doench2014OnTarget)]),
#          hours == saturation,
#          arm == 'SD-Met-Cys-Ura+Gal', phenotype == 'Beneficial', orf_type == 'Proto-gene') %>%
#   group_by(orf_type) %>%
#   data.frame()
# 
# data.sum %>%
#   filter(orf_name %in% unique(data.crsp$orf_name[data.crsp$Hsu2013 >= 95 & data.crsp$Doench2014OnTarget >= 0.7 &
#                                                    !is.na(data.crsp$Doench2014OnTarget)]),
#          hours == saturation) %>%
#   group_by(arm, stage, orf_type, phenotype) %>%
#   dplyr::count()
# 
# data.sum %>%
#   filter(hours == saturation,
#          arm == 'SD-Met-Cys-Ura+Gal') %>%
#   group_by(orf_type, phenotype) %>%
#   dplyr::count()
# 
# data.sum %>%
#   filter(hours == saturation) %>%
#   group_by(arm, stage, orf_type, phenotype) %>%
#   dplyr::count() %>%
#   data.frame()
# 
# pgs %>%
#   filter(orf_name %in% unique(data.crsp$orf_name[data.crsp$Hsu2013 >= 95 & data.crsp$Doench2014OnTarget >= 0.7 &
#                                                    !is.na(data.crsp$Doench2014OnTarget)])) %>%
#   dplyr::count()
# 
# 428-372
# 
# 
# plot.theme <- theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# 
# plot.onandoff <- data.crsp2 %>%
#   ggplot(aes(x = Doench2014OnTarget, y = Hsu2013)) +
#   geom_point(size = 0.2) +
#   plot.theme
# 
# plot.on <- data.crsp2 %>%
#   ggplot() +
#   geom_line(aes(x = Doench2014OnTarget), stat = 'density') +
#   plot.theme + theme_minimal() +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank())
# 
# plot.off <- data.crsp2 %>%
#   ggplot() +
#   geom_line(aes(y = Hsu2013), stat = 'density') +
#   plot.theme + theme_minimal() +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank())
# 
# ggpubr::ggarrange(plot.on,NULL,plot.onandoff,plot.off,nrow = 2,ncol = 2,align = 'hv',
#                   widths = c(3, 1), heights = c(1, 3))
# 


##### GO/KEGG ENRICHMENT
goe <- data.frame()
kegg <- data.frame()
for (a in unique(data.sum$arm)) {
  allgenes <- unique(data.sum$orf_name[data.sum$arm == a])
  allgenes <- bitr(allgenes, fromType = "ORF",
                   toType = c("ENTREZID","GENENAME","ENSEMBL"),
                   OrgDb = org.Sc.sgd.db)
  allgenes <- allgenes[!is.na(allgenes$ENSEMBL),]
  for (s in unique(data.sum$stage[data.sum$arm == a])) {
    for (p in unique(data.sum$phenotype[data.sum$stage == s & data.sum$arm == a])) {
      temp.deg <- data.sum$orf_name[data.sum$stage == s & data.sum$phenotype == p & data.sum$arm == a]
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
        cat(sprintf('There are no GO term enrichment for %s ORFs in %s %s.\n',p,a,s))
      } else{
        cat(sprintf('GO term enrichment for %s ORFs in %s %s are:\n%s\n',
                    p,a,s,
                    paste(temp.goe$Description,collapse = ', ')))
        goe <- rbind(goe, data.frame(temp.goe, arm = a, stage = s, phenotype = p))
      }
      
      temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                              universe     = allgenes$ENSEMBL,
                              organism     = 'sce',
                              pvalueCutoff = 0.05)
      if (dim(temp.kegg)[1] == 0) {
        cat(sprintf('There are no KEGG pathway enriched for %s ORFs in %s %s.\n',p,a,s))
      } else{
        cat(sprintf('KEGG pathway enrichment for %s ORFs in %s %s are:\n%s\n',
                    p,a,s,
                    paste(temp.kegg$Description,collapse = ', ')))
        kegg <- rbind(kegg, data.frame(temp.kegg, arm = a, stage = s, phenotype = p))
      }
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

goe <- goe[order(goe$arm,goe$stage,goe$phenotype,-goe$GeneRatio,-goe$Count,goe$qvalue),]
kegg <- kegg[order(kegg$arm,kegg$stage,kegg$phenotype,-kegg$GeneRatio,-kegg$Count,kegg$qvalue),]

write.csv(goe, file = sprintf('%s/%s/go_enrichments.csv', res_path, expt.name))
write.csv(kegg, file = sprintf('%s/%s/kegg_enrichments.csv', res_path, expt.name))

##### GO-KEGG PLOTS
# head(goe)
# goe$phenotype <- factor(goe$phenotype, levels = c('Beneficial','Neutral','Deleterious','Dead'))
# goe$stage <- factor(goe$stage, levels = c('Pre-Screen #1','Pre-Screen #2','Final Screen'))
# kegg$phenotype <- factor(kegg$phenotype, levels = c('Beneficial','Neutral','Deleterious','Dead'))
# kegg$stage <- factor(kegg$stage, levels = c('Pre-Screen #1','Pre-Screen #2','Final Screen'))

##### DEATH ANALYSIS
data.ded <- merge(merge(data.sum[data.sum$stage == 'Pre-Screen #1' & data.sum$hours == data.sum$saturation, c('arm','orf_name','strain_id','fitness','cs','orf_type','phenotype')],
                        data.sum[data.sum$stage == 'Pre-Screen #2' & data.sum$hours == data.sum$saturation, c('arm','orf_name','strain_id','fitness','cs','orf_type','phenotype')],
                        by = c('arm','orf_name','strain_id','orf_type'), suffixes = c('_PS1','_PS2'), all = T),
                  data.sum[data.sum$stage == 'Final Screen' & data.sum$hours == data.sum$saturation, c('arm','orf_name','strain_id','fitness','cs','orf_type','phenotype')],
                  by = c('arm','orf_name','strain_id','orf_type'), suffixes = c('','_FS'), all = T)
head(data.ded)
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

data.ded <- merge(data.ded[data.ded$arm == 'SD-Met-Cys-Ura+Gal',-1], data.ded[data.ded$arm == 'SC-Ura+Gal',-1],
                  by = c('orf_name','strain_id','orf_type'), suffixes = c('_MM','_PM'), all = T)
# data.ded$diff_PS1 <- data.ded$fitness_PS1_MM - data.ded$fitness_PS1_PM
# data.ded$diff_PS2 <- data.ded$fitness_PS2_MM - data.ded$fitness_PS2_PM
# data.ded$diff_FS <- data.ded$fitness_FS_MM - data.ded$fitness_FS_PM
head(data.ded)

data.ded %>%
  group_by(orf_type, phenotype_PS2_MM) %>%
  dplyr::count()

##### COUNTS OF PGS BENEFICIAL, DELETERIOUS, NEUTRAL AND DEAD
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
data.cnts$stage <- factor(data.cnts$stage, levels = stages)

data.cnts$arm[str_detect(data.cnts$level,'MM')] <- 'SD-Met-Cys-Ura+Gal'
data.cnts$arm[str_detect(data.cnts$level,'PM')] <- 'SC-Ura+Gal'
data.cnts$arm <- factor(data.cnts$arm, levels = arms)

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
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/%s/PHENOTYPE_COUNTS.jpg",fig_path, expt.name), plot.pheno.cnts,
       height = one.c, width = two.c, units = 'mm',
       dpi = 600)  

# write.csv(data.cnts, file = sprintf('%s/%s/phenotype_counts.csv',res_path,expt.name))

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
data.or$arm <- factor(data.or$arm, levels = c('SC-Ura+Gal','SD-Met-Cys-Ura+Gal','Differential'))

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
  geom_text(aes(x = phenotype, y = 8, label = label), size = 2, col = 'red') +
  labs(x='Phenotype',y='Odds Ratio') +
  # scale_y_continuous(trans="log10", breaks=c(0.2,0.5,1,2,5,10)) + 
  # annotation_logticks(sides="l") +
  facet_wrap(.~arm*stage, nrow = 2) +
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

# write.csv(data.or, file = sprintf('%s/%s/phenotype_odds.csv',res_path,expt.name))

##### SAVING ALL DATA
# save(data, data.cnts, data.ded, data.diff, data.jm, data.lim, data.mad,
#      data.or, data.sap, data.sul, data.sum, diff.dist, goe, kegg,
#      orfs, orfs.jm, orfs.sap, orfs.sul, pgs,
#      file = sprintf('%s/%s/all_data.RData',res_path,expt.name))
