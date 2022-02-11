##### JAKES MUTANT EXPERIMENT - FINAL FIGURES
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 05/19/2021 

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(ggridges)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(ggtext)
library(plotly)
library(scales)
library(reshape2)
library(locfit)
library(growthcurver)
library(rstatix)
library(gtools)
library(effsize)

out_path <- "~/R/Projects/methionine/data"
fig_path <- "~/R/Projects/methionine/figures"
res_path <- "~/R/Projects/methionine/results"

expt.name <- "jakes"

source("~/R/Projects/methionine/functions/colorstrip.R")

##### LOAD DATA
load(file = sprintf('%s/%s/colonysizes2.RData', out_path, expt.name))
load(file = sprintf('%s/%s/anova_results.RData', out_path, expt.name))
# load(file = sprintf('%s/%s/anova_results2.RData', out_path, expt.name))
load(file = sprintf('%s/%s/anova_results2_nopilot.RData', out_path, expt.name))

##### LEVELS
attempt.levels <- c('pilot', 'copy1', 'copy2')
strain.levels <- c('BY4742', 'BY4741', 'FY4', 'met15', 'met3', 'met2', 'met6', 'met13', 'cys4', 'met12', 'str3', 'yll')
condition.levels <- c('YPDA', 'SD-Met+Glu')

# STRAIN LABELES
strain.labs <- c('BY4742', 'BY4741', 'FY4', 'FY4-*met15Δ*', 'FY4-*met3Δ*', 'FY4-*met2Δ*', 'FY4-*met6Δ*', 'FY4-*met13Δ*','FY4-*cys4Δ*',
                 'FY4-*met12Δ*', 'FY4-*str3Δ*', 'FY4-*yll058wΔ*')
strain.labs.parse <- c('BY4742', 'BY4741', 'FY4', 'FY4-italic(met15Δ)', 'FY4-italic(met3Δ)', 'FY4-italic(met2Δ)',
                       'FY4-italic(met6Δ)', 'FY4-italic(met13Δ)','FY4-italic(cys4Δ)',
                       'FY4-italic(met12Δ)', 'FY4-italic(str3Δ)', 'FY4-italic(yll058wΔ)')
strain.labs2 <- data.frame(orf_name = strain.levels, labels = strain.labs, parsed = strain.labs.parse)
strain.labs2$auxotrophy[strain.labs2$orf_name %in% c('BY4742','FY4')] <- 'Prototroph'
strain.labs2$auxotrophy[strain.labs2$orf_name %in% c('BY4741','met15','met3','met2','met6','met13','cys4')] <- 'Presumed Auxotroph'
strain.labs2$auxotrophy[strain.labs2$orf_name %in% c('met12')] <- 'Prototroph'
strain.labs2$auxotrophy[strain.labs2$orf_name %in% c('str3','yll')] <- 'Unknown'
strain.labs2$auxotrophy <- factor(strain.labs2$auxotrophy, levels = c('Prototroph', 'Presumed Auxotroph', 'Unknown'))

data$condition <- as.character(data$condition)
data$condition[data$condition == 'SD-MET+Glucose'] <- 'SD-Met+Glu'
anova.res$condition <- as.character(anova.res$condition)
anova.res$condition[anova.res$condition == 'SD-MET+Glucose'] <- 'SD-Met+Glu'
anova.res2$condition <- as.character(anova.res2$condition)
anova.res2$condition[anova.res2$condition == 'SD-MET+Glucose'] <- 'SD-Met+Glu'

data$attempt <- factor(data$attempt, levels = attempt.levels)
data$orf_name <- factor(data$orf_name, levels = strain.levels)
data$condition <- factor(data$condition, levels = condition.levels)

anova.res$strain1 <- factor(anova.res$strain1, levels = strain.levels)
anova.res$strain2 <- factor(anova.res$strain2, levels = strain.levels)
anova.res$condition <- factor(anova.res$condition, levels = condition.levels)

anova.res2$strain1 <- factor(anova.res2$strain1, levels = strain.levels)
anova.res2$strain2 <- factor(anova.res2$strain2, levels = strain.levels)
anova.res2$condition <- factor(anova.res2$condition, levels = condition.levels)

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### STAT LABELS
anova.res$label[anova.res$rcs_between > 0.05] <- 'ns'
anova.res$label[anova.res$rcs_between <= 0.05] <- '*'
anova.res$label[anova.res$rcs_between <= 0.01] <- '**'
anova.res$label[anova.res$rcs_between <= 0.001] <- '***'
anova.res$label[anova.res$rcs_between <= 0.0001] <- '****'

anova.res2$label[anova.res2$rcs_kw > 0.05] <- 'ns'
anova.res2$label[anova.res2$rcs_kw <= 0.05] <- '*'
anova.res2$label[anova.res2$rcs_kw <= 0.01] <- '**'
anova.res2$label[anova.res2$rcs_kw <= 0.001] <- '***'
anova.res2$label[anova.res2$rcs_kw <= 0.0001] <- '****'

data$attempt_label[data$attempt == 'pilot'] <- 'Replicate 0'
data$attempt_label[data$attempt == 'copy1'] <- 'Replicate 1'
data$attempt_label[data$attempt == 'copy2'] <- 'Replicate 2'

anova.res$attempt_label[anova.res$attempt == 'pilot'] <- 'Replicate 0'
anova.res$attempt_label[anova.res$attempt == 'copy1'] <- 'Replicate 1'
anova.res$attempt_label[anova.res$attempt == 'copy2'] <- 'Replicate 2'

##### WITHOUT YLL
temp <- data[data$time == 't_final' & data$attempt != 'pilot',] %>%
  group_by(orf_name, condition) %>%
  summarize(f = median(relative_cs, na.rm = T), .groups = 'keep') %>% data.frame()
temp.es <- data.frame()
# data.es <- data.frame()
for (o in unique(temp$orf_name)) {
  es <- (temp$f[temp$orf_name == o & temp$condition == 'SD-Met+Glu'] - temp$f[temp$orf_name == o & temp$condition == 'YPDA'])/
    temp$f[temp$orf_name == o & temp$condition == 'YPDA']
  temp.es <- rbind(temp.es, data.frame(orf_name = o, effsize = es, PM = temp$f[temp$orf_name == o & temp$condition == 'YPDA'],
                                       MM = temp$f[temp$orf_name == o & temp$condition == 'SD-Met+Glu']))
  # data.es <- rbind(data.es, data.frame(orf_name = o, label = strain.labs2$labels[strain.labs2$orf_name == o],
  #                                      parsed = strain.labs2$parsed[strain.labs2$orf_name == o],
  #                                      PM = temp$f[temp$orf_name == o & temp$condition == 'YPDA'],
  #                                      MM = temp$f[temp$orf_name == o & temp$condition == 'SD-Met-Cys+Glu'], effsize = es))
}
data.es <- merge(temp.es, strain.labs2, by = 'orf_name')

plot.rcs.box.kw <- data[data$time == 't_final' & data$attempt != 'pilot' & data$orf_name != 'yll',] %>%
  group_by(attempt, condition, time, orf_name, bio_rep) %>%
  # summarise(average = mean(average, na.rm = T), relative_cs = mean(relative_cs, na.rm = T)) %>%
  data.frame() %>%
  ggplot(aes(x = condition, y = relative_cs)) +
  geom_boxplot(#aes(fill = bio_rep), 
               size = 0.3, outlier.shape = NA) +
  # geom_jitter(aes(col = bio_rep, shape = attempt), size = 1) +
  # geom_jitter(size = 1) +
  stat_compare_means(aes(label = ..p.signif..), method = 'kruskal',
                      label.x = 1.4, label.y = 1.45, size = 2.2, col = 'red') +
  geom_text(data = temp.es[temp.es$orf_name != 'yll',], aes(x = 1.45, y = 1.35, label = sprintf('%0.2f%%',effsize*100)), size = 2) +
  # geom_text(data = anova.res2[anova.res2$strain2 == 'FY4' & 
  #                              anova.res2$time == 't_final' & anova.res2$strain1 != 'yll',],
  #           aes(x = strain1, y = 1.5, label = label), size = 2.2, col = 'red') +
  # geom_text(data = anova.res2[anova.res2$strain2 == 'FY4' &
  #                               anova.res2$time == 't_final' & anova.res2$strain1 != 'yll',],
  #           aes(x = strain1, y = 1.4, label = sprintf('%0.2f%%',effect_size*100)),
  #           size = 2.2) +
  # scale_x_discrete(labels = strain.labs[-10]) +
  scale_y_continuous(minor_breaks = seq(-2,2,0.05)) +
  facet_wrap(.~orf_name, nrow = 2,
             labeller = labeller(orf_name = c('FY4'='FY4',
                                              'met12'='FY4-*met12Δ*',
                                              'str3'='FY4-*str3Δ*',
                                              'met3'='FY4-*met3Δ*',
                                              'met15'='FY4-*met15Δ*',
                                              'met2'='FY4-*met2Δ*',
                                              'met6'='FY4-*met6Δ*', 
                                              'met13'='FY4-*met13Δ*',
                                              'cys4'='FY4-*cys4Δ*',
                                              'yll'='FY4-*yll058wΔ*',
                                              'BY4742'='BY4742',
                                              'BY4741'='BY4741'))) +
  labs(x = 'Media Type', y = 'Relative Colony Size') +
  scale_fill_manual(name = 'Biological Replicate',
                    values = c("1" = "#673AB7",
                               "2" = "#009688",
                               "3" = "#607D8B",
                               "4" = "#FFC107")) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.x = ggtext::element_markdown(angle = 40, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))
# plot.rcs.box.kw <- colorstrip(plot.rcs.box.kw,c("#00796B","#212121"))
fig3A <- plot.rcs.box.kw
save(fig3A, file = 'figures/final/fig3A.RData')
ggsave(sprintf("%s/%s/final/RELATIVE_COLONY_SIZE_BOX_KW.jpg",fig_path, expt.name),
       plot.rcs.box.kw,
       height = one.5c, width = one.5c, units = 'mm',
       dpi = 600)
# write.csv(anova.res2[anova.res2$strain2 == 'FY4' & anova.res2$time == 't_final',],
#           file = sprintf('%s/%s/final/RELATIVE_COLONY_SIZE_KW.csv', res_path, expt.name))


##### EFFECT SIZE PLOT
fig3A.2 <- data.es[data.es$orf_name != 'yll',]  %>%
  ggplot(aes(x = PM, y = MM)) +
  geom_abline(linetype = 'dashed', col = 'black') +
  annotate("rect", xmin = 0.9, xmax = 1, ymin = 0.3, ymax = 0.6,
           alpha = .1, fill = "blue", col = "blue", size = 0.2, linetype = 'dashed') +
  # annotate("text", x = 1, y = 0.6,
  #           label = 'Unexpected', hjust = 0.5, vjust = -0.5, size = 2, col = 'red') +
  geom_point(size = 3) +
  # geom_richtext(aes(label = label), size = 1.4, nudge_y = 0.05) +
  geom_label_repel(aes(label = parsed, fill = auxotrophy),
                   parse = T, size = 2.5, col = 'black') +
  labs(x = 'Median Relative Fitness in YPDA',
       y = 'Median Relative Fitness in SD-Met-Cys+Glu') +
  scale_fill_manual(name = 'Methionine Auxotrophy',
                    values = c('Prototroph' = '#4CAF50',
                               'Presumed Auxotroph' = '#F44336',
                               'Unknown' = '#03A9F4'),
                    labels = c('Prototroph' = 'Prototroph',
                               'Presumed Auxotroph' = 'Presumed Auxotroph',
                               'Unknown' = 'Unknown')) +
  scale_x_continuous(breaks = seq(-1,1,0.2)) + 
  scale_y_continuous(breaks = seq(-1,1,0.2)) +
  coord_cartesian(xlim = c(0, 1.1),
                  ylim = c(0, 1.1)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              face = 'bold',
                                              margin = margin(0.1,0,0.1,0, "mm")))
save(fig3A.2, file = 'figures/final/fig3A.2.RData')

##### SUPPLEMENTARY FIGURES #####
##### DENSITY PLOT WITH ANOVA RESULTS
plot.rcs.den.aov <- data[data$attempt != 'pilot' & data$time == 't_final' & data$orf_name != 'yll',] %>%
  ggplot(aes(x = relative_cs, y = orf_name)) +
  geom_density_ridges(quantile_lines = TRUE,
                      aes(fill = bio_rep),
                      scale = 2, alpha = 0.7, size = 0.2,
                      vline_size = 0.2, vline_color = "black") +
  geom_text(data = anova.res[anova.res$attempt != 'pilot' &
                               anova.res$strain2 == 'FY4' & 
                               anova.res$time == 't_final' &
                               anova.res$strain1 != 'yll',],
            aes(y = strain1, x = 1.5, label = label), size = 2, col = 'red') +
  labs(x = 'Relative Colony Size',
       y = 'Strain') +
  scale_y_discrete(labels = strain.labs[-12]) +
  scale_fill_manual(name = 'Biological Replicate',
                    values = c("1" = "#673AB7",
                               "2" = "#009688",
                               "3" = "#607D8B",
                               "4" = "#FFC107")) +
  facet_wrap(.~condition*attempt_label, nrow = 1) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = ggtext::element_markdown(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
# plot.rcs.den.aov <- colorstrip(plot.rcs.den.aov,c("#C2185B","#C2185B","#448AFF","#448AFF"))
ggsave(sprintf("%s/%s/final/RELATIVE_COLONY_SIZE_DENSITY_ANOVA.jpg",fig_path, expt.name), plot.rcs.den.aov,
       height = 70, width = two.c, units = 'mm',
       dpi = 600)


##### ONLY YLL
plot.rcs.box.kw <- data[data$time == 't_final' & data$attempt != 'pilot' & data$orf_name %in% c('FY4','yll'),] %>%
  group_by(attempt, condition, time, orf_name, bio_rep) %>%
  # summarise(average = mean(average, na.rm = T), relative_cs = mean(relative_cs, na.rm = T)) %>%
  data.frame() %>%
  ggplot(aes(x = condition, y = relative_cs)) +
  geom_boxplot(aes(fill = bio_rep), size = 0.3, outlier.shape = NA) +
  # geom_jitter(aes(col = bio_rep, shape = attempt), size = 1) +
  # geom_jitter(size = 1) +
  stat_compare_means(aes(label = ..p.signif..), method = 'kruskal',
                     label.x = 1.4, label.y = 1.45, size = 2.2, col = 'red') +
  geom_text(data = temp.es[temp.es$orf_name %in% c('FY4','yll'),], aes(x = 1.45, y = 1.35, label = sprintf('%0.2f%%',effsize*100)), size = 2) +
  # geom_text(data = anova.res2[anova.res2$strain2 == 'FY4' & 
  #                              anova.res2$time == 't_final' & anova.res2$strain1 != 'yll',],
  #           aes(x = strain1, y = 1.5, label = label), size = 2.2, col = 'red') +
  # geom_text(data = anova.res2[anova.res2$strain2 == 'FY4' &
  #                               anova.res2$time == 't_final' & anova.res2$strain1 != 'yll',],
  #           aes(x = strain1, y = 1.4, label = sprintf('%0.2f%%',effect_size*100)),
  #           size = 2.2) +
  # scale_x_discrete(labels = strain.labs[-10]) +
  scale_y_continuous(minor_breaks = seq(-2,2,0.05)) +
  facet_wrap(.~orf_name, nrow = 1,
             labeller = labeller(orf_name = c('FY4'='FY4',
                                              'met12'='FY4-*met12Δ*',
                                              'str3'='FY4-*str3Δ*',
                                              'met3'='FY4-*met3Δ*',
                                              'met15'='FY4-*met15Δ*',
                                              'met2'='FY4-*met2Δ*',
                                              'met6'='FY4-*met6Δ*', 
                                              'met13'='FY4-*met13Δ*',
                                              'cys4'='FY4-*cys4Δ*',
                                              'yll'='FY4-*yll058wΔ*',
                                              'BY4742'='BY4742',
                                              'BY4741'='BY4741'))) +
  labs(x = 'Media Type', y = 'Relative Colony Size') +
  scale_fill_manual(name = 'Biological Replicate',
                    values = c("1" = "#673AB7",
                               "2" = "#009688",
                               "3" = "#607D8B",
                               "4" = "#FFC107")) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.x = ggtext::element_markdown(angle = 40, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              face = 'bold',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))
# plot.rcs.box.kw <- colorstrip(plot.rcs.box.kw,c("#212121","#00796B"))
fig6B <- plot.rcs.box.kw
save(fig6B, file = 'figures/final/fig6B.RData')
# ggsave(sprintf("%s/%s/final/RELATIVE_COLONY_SIZE_BOX_KW_YLL.jpg",fig_path, expt.name),
#        plot.rcs.box.kw,
#        height = 70, width = one.c, units = 'mm',
#        dpi = 600)
# write.csv(anova.res2[anova.res2$strain2 == 'FY4' & anova.res2$time == 't_final',],
#           file = sprintf('%s/%s/final/RELATIVE_COLONY_SIZE_KW.csv', res_path, expt.name))

## DENSITY PLOT WITH ANOVA RESULTS WITH YLL
plot.rcs.den.aov <- data[data$attempt != 'pilot' & data$time == 't_final' & data$orf_name %in% c('FY4','yll'),] %>%
  ggplot(aes(x = relative_cs, y = orf_name)) +
  geom_density_ridges(quantile_lines = TRUE,
                      aes(fill = bio_rep),
                      scale = 2, alpha = 0.7, size = 0.2,
                      vline_size = 0.2, vline_color = "black") +
  geom_text(data = anova.res[anova.res$attempt != 'pilot' &
                               anova.res$strain2 == 'FY4' & 
                               anova.res$time == 't_final' &
                               anova.res$strain1 == 'yll',],
            aes(y = strain1, x = 1.5, label = label), size = 2, col = 'red') +
  labs(x = 'Relative Colony Size',
       y = 'Strain') +
  scale_y_discrete(labels = strain.labs[c(1,10)]) +
  scale_fill_manual(name = 'Biological Replicate',
                    values = c("1" = "#673AB7",
                               "2" = "#009688",
                               "3" = "#607D8B",
                               "4" = "#FFC107")) +
  facet_wrap(.~condition*attempt_label, nrow = 1) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = ggtext::element_markdown(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
plot.rcs.den.aov <- colorstrip(plot.rcs.den.aov,c("#C2185B","#C2185B","#448AFF","#448AFF"))
ggsave(sprintf("%s/%s/final/RELATIVE_COLONY_SIZE_DENSITY_ANOVA_YLL.jpg",fig_path, expt.name), plot.rcs.den.aov,
       height = 70, width = two.c, units = 'mm',
       dpi = 600)

##### GROWTH CURVES
head(data)
data %>%
  filter(orf_name %in% c('FY4','met15','BY4741'), data$attempt != 'pilot') %>%
  ggplot(aes(x = hours, y = average, col = condition, group = condition)) +
  geom_line() +
  facet_wrap(.~orf_name)
