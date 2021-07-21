##### JAKES MUTANT EXPERIMENT - FIGURES
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 05/05/2021 

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

load(sprintf("%s/%s/colonysizes.RData", out_path, expt.name))

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### REMOVING BORDER RELATED DATA
data <- data[data$orf_name != 'BOR',]

##### VARIABLE LEVELS
attempt.levels <- c('pilot', 'copy1', 'copy2')
strain.levels <- c('FY4', 'met15', 'BY4742', 'BY4741', 'met2', 'met3', 'met6', 'met12', 'met13', 'yll', 'str3', 'cys4')
condition.levels <- c('YPDA', 'SD-MET+Glucose')

data$attempt <- factor(data$attempt, levels = attempt.levels)
data$orf_name <- factor(data$orf_name, levels = strain.levels)
data$condition <- factor(data$condition, levels = condition.levels)

##### ALL GROWTH CURVES
# all.gc <- data %>%
#   ggplot(aes(x = hours, y = average, col = bio_rep)) +
#   geom_point() +
#   geom_smooth(method = 'loess') +
#   # scale_x_continuous(breaks = unique(data$hours), minor_breaks = NULL) +
#   scale_color_discrete(name = 'Biological Replicate') +
#   labs(y = 'Colony Size (pixels)',
#        x = 'Time (hours)') +
#   facet_wrap(.~attempt*condition*orf_name,
#              ncol = 12) +
#   theme_linedraw() +
#   theme(axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%s/%s/ALL_GROWTH_CURVES_POINTS.jpg",fig_path, expt.name), all.gc,
#        height = 300, width = 600, units = 'mm',
#        dpi = 600)

##### CLEAN DATA
data <- data[!(data$condition == 'YPDA' & data$hours > 100),]
data$time[data$hours < 10] <- 't_0'
data$time[data$hours > 10] <- 't_final'

##### T_0 and T_Final CORR
plot.rcs.t0tf.cor <- merge(data[data$time == 't_0',], data[data$time == 't_final',],
      by = c('attempt','condition','orf_name','bio_rep','pos'),
      suffixes = c('.t_0','.t_final')) %>%
  ggplot(aes(x = relative_cs.t_final, y = relative_cs.t_0)) +
  geom_point(aes(shape = condition, col = attempt), size = 1) +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', size = 1.2, label.x = 0.4) +
  labs(x = 'T_Final Relative Colony Size',
       y = 'T_0 Relative Colony Size') +
  scale_color_discrete(name = 'Experimental\nReplicate') +
  scale_shape_discrete(name = 'Media') +
  # facet_wrap(.~condition, ncol = 12) +
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
ggsave(sprintf("%s/%s/COLONY_SIZE_T0TF_CORR.jpg",fig_path, expt.name), plot.rcs.t0tf.cor,
       height = 100, width = 100, units = 'mm',
       dpi = 600)

##### ANOVA AND EFFECT SIZE
strain.pairs <- permutations(n = length(unique(data$orf_name)), r = 2,
                             v = as.character(unique(data$orf_name)),
                             repeats.allowed = F) %>%
  data.frame(stringsAsFactors = F)

anova.res <- NULL
for (a in unique(data$attempt)) {
  for (c in unique(data$condition[data$attempt == a])) {
    for (h in unique(data$hours[data$attempt == a & data$condition == c])) {
      data$relative_cs[data$attempt == a & data$hours == h & data$condition == c] <- 
        data$average[data$attempt == a & data$hours == h & data$condition == c]/
        median(data$average[data$attempt == a & data$hours == h & data$condition == c & data$orf_name == 'FY4'], na.rm = T)
      for (i in seq(1,dim(strain.pairs)[1])) {
        o1 <- strain.pairs[i,1]
        o2 <- strain.pairs[i,2]
        
        res.aov <- data[data$attempt == a & data$hours == h & data$condition == c & data$orf_name %in% c(o1, o2),] %>%
          anova_test(average ~ orf_name * bio_rep)
        res.rcs.aov <- data[data$attempt == a & data$hours == h & data$condition == c & data$orf_name %in% c(o1, o2),] %>%
          anova_test(relative_cs ~ orf_name * bio_rep)
        
        cs1 <- data$average[data$attempt == a & data$hours == h & data$condition == c & data$orf_name == o1]
        cs2 <- data$average[data$attempt == a & data$hours == h & data$condition == c & data$orf_name == o2]
        effsize.temp <- wilcox.test(cs1, cs2, conf.int = T)
        emp_effsize <- (median(cs1, na.rm = T) - median(cs2, na.rm = T))/median(cs2, na.rm = T)
        delta.temp <- cliff.delta(cs1, cs2)
        relative_cs <- median(cs1, na.rm = T)/median(cs2, na.rm = T)
        
        anova.res <- rbind(anova.res, cbind(a, c, h, o1, o2,
                                            res.aov$p[res.aov$Effect == 'orf_name'],
                                            res.aov$p[res.aov$Effect == 'orf_name:bio_rep'],
                                            res.rcs.aov$p[res.aov$Effect == 'orf_name'],
                                            res.rcs.aov$p[res.aov$Effect == 'orf_name:bio_rep'],
                                            effsize.temp$estimate,
                                            delta.temp$estimate, delta.temp$magnitude,
                                            emp_effsize, relative_cs))
      }
    }
  }
}
colnames(anova.res) <- c('attempt','condition','hours','strain1','strain2','between','within','rcs_between','rcs_within',
                         'wilcox_estimate','cliffs_estimate','cliffs_magnitude','effect_size','relative_cs')
anova.res <- data.frame(anova.res, stringsAsFactors = F)
anova.res$hours <- as.numeric(anova.res$hours)
anova.res$between <- as.numeric(anova.res$between)
anova.res$within <- as.numeric(anova.res$within)
anova.res$rcs_between <- as.numeric(anova.res$rcs_between)
anova.res$rcs_within <- as.numeric(anova.res$rcs_within)
anova.res$wilcox_estimate <- as.numeric(anova.res$wilcox_estimate)
anova.res$cliffs_estimate <- as.numeric(anova.res$cliffs_estimate)
anova.res$effect_size <- as.numeric(anova.res$effect_size)
anova.res$relative_cs <- as.numeric(anova.res$relative_cs)
anova.res$between <- p.adjust(anova.res$between, method = 'BH')
anova.res$within <- p.adjust(anova.res$within, method = 'BH')
anova.res$rcs_between <- p.adjust(anova.res$rcs_between, method = 'BH')
anova.res$rcs_within <- p.adjust(anova.res$rcs_within, method = 'BH')

anova.res$attempt <- factor(anova.res$attempt, levels = attempt.levels)
anova.res$strain1 <- factor(anova.res$strain1, levels = strain.levels)
anova.res$strain2 <- factor(anova.res$strain2, levels = strain.levels)
anova.res$condition <- factor(anova.res$condition, levels = condition.levels)
anova.res$cliffs_magnitude <- factor(anova.res$cliffs_magnitude, levels = c('negligible','small','medium','large'))

anova.res$time[anova.res$hours < 10] <- 't_0'
anova.res$time[anova.res$hours > 10] <- 't_final'

# write.csv(anova.res, file = sprintf('%s/%s/COLONY_SIZE_ANOVA_RESULTS.csv', res_path, expt.name))
# save(anova.res, file = sprintf('%s/%s/anova_results.RData', out_path, expt.name))
# save(data, file = sprintf('%s/%s/colonysizes2.RData', out_path, expt.name))
load(file = sprintf('%s/%s/colonysizes2.RData', out_path, expt.name))
load(file = sprintf('%s/%s/anova_results.RData', out_path, expt.name))

plot.anova.res <- ggplot(anova.res) +
  geom_tile(aes(x = strain1, y = strain2, fill = between <= 0.05),
            col = 'black') +
  facet_wrap(.~attempt*condition*time, ncol = 2) +
  scale_fill_manual(name = 'Significant Anova',
                    breaks = c('TRUE', 'FALSE'),
                    values = c('TRUE' = '#4CAF50', 'FALSE' = '#536DFE')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")),
        panel.grid = element_blank())
ggsave(sprintf("%s/%s/COLONY_SIZE_ANOVA.jpg",fig_path, expt.name), plot.anova.res,
       height = 200, width = 150, units = 'mm',
       dpi = 600) 

##### COLONY SIZE COMPARISONS
data %>%
  group_by(attempt, condition, time, orf_name, bio_rep) %>%
  count() %>%
  data.frame()

# plot.cs.vio.1 <- data[!(data$condition == 'YPDA' & data$hours > 100) & data$hours == 0,] %>%
#   ggplot(aes(x = orf_name, y = average)) +
#   geom_violin(aes(fill = bio_rep),lwd = 0.1) +
#   geom_point(data = anova.res[anova.res$strain2 == 'FY4' &
#                                 anova.res$time == 't_0' & anova.res$between <= 0.05,],
#              aes(x = strain1, y = 610),
#              shape = 8, size = 1, col = 'red') +
#   # geom_jitter() +
#   scale_fill_discrete(name = 'Biological Replicate') +
#   facet_wrap(.~attempt*condition*time, ncol = 1) +
#   labs(x = 'Strain',
#        y = 'Colony Size (pixels)') +
#   # scale_y_log10() +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm")))

plot.cs.vio.2 <- data[!(data$condition == 'YPDA' & data$hours > 100) & data$hours != 0,] %>%
  ggplot(aes(x = orf_name, y = average)) +
  geom_violin(aes(fill = bio_rep),lwd = 0.1) +
  geom_point(data = anova.res[anova.res$strain2 == 'FY4' &
                                anova.res$time == 't_final' & anova.res$between <= 0.05,],
             aes(x = strain1, y = 10000),
             shape = 8, size = 1, col = 'red') +
  # geom_jitter() +
  scale_fill_discrete(name = 'Biological Replicate') +
  facet_wrap(.~attempt*condition*time, nrow = 1) +
  labs(x = 'Strain',
       y = 'Colony Size (pixels)') +
  # scale_y_log10() +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/%s/COLONY_SIZE_VIOLIN_ANOVA.jpg",fig_path, expt.name), plot.cs.vio.2,
       height = 80, width = 200, units = 'mm',
       dpi = 600) 

plot.den.anova <- data[data$time == 't_final',] %>%
  ggplot(aes(x = average, y = orf_name)) +
  geom_density_ridges(quantile_lines = TRUE,
                      aes(fill = bio_rep),
                      scale = 2, alpha = 0.7, size = 0.2,
                      vline_size = 0.2, vline_color = "black") +
  geom_point(data = anova.res[anova.res$strain2 == 'FY4' & 
                                anova.res$time == 't_final' & anova.res$between < 0.05,],
             aes(y = strain1, x = 10000), shape = 8, size = 1, col = 'red') +
  labs(x = 'Colony Size',
       y = 'Strain') +
  scale_fill_discrete(name = 'Biological Replicate') +
  facet_wrap(.~attempt*condition*time, nrow = 1) +
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
ggsave(sprintf("%s/%s/COLONY_SIZE_DENSITY_ANOVA.jpg",fig_path, expt.name), plot.den.anova,
       height = 80, width = 200, units = 'mm',
       dpi = 600)

# plot.cs.vio <- ggpubr::ggarrange(plot.cs.vio.1, plot.cs.vio.2,
#                                  nrow = 1,
#                                  common.legend = T, legend = 'bottom')
# ggsave(sprintf("%s/%s/COLONY_SIZE_VIOLIN.jpg",fig_path, expt.name), plot.cs.vio,
#        height = 200, width = 150, units = 'mm',
#        dpi = 600) 

# plot.cs.den.1 <- data[!(data$condition == 'YPDA' & data$hours > 100) & data$hours == 0,] %>%
#   ggplot(aes(x = average, y = orf_name, fill = bio_rep)) +
#   geom_density_ridges(quantile_lines = TRUE,
#                       scale = 3, alpha = 0.7, size = 0.3,
#                       vline_size = 0.2, vline_color = "black") +
#   labs(x = 'Colony Size (pixels)',
#        y = 'Strains') +
#   scale_fill_discrete(name = 'Biological Replicate') +
#   facet_wrap(.~attempt*condition*time, ncol = 1) +
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
# 
# plot.cs.den.2 <- data[!(data$condition == 'YPDA' & data$hours > 100) & data$hours != 0,] %>%
#   ggplot(aes(x = average, y = orf_name, fill = bio_rep)) +
#   geom_density_ridges(quantile_lines = TRUE,
#                       scale = 3, alpha = 0.7, size = 0.3,
#                       vline_size = 0.2, vline_color = "black") +
#   labs(x = 'Colony Size (pixels)',
#        y = 'Strains') +
#   scale_fill_discrete(name = 'Biological Replicate') +
#   facet_wrap(.~attempt*condition*time, ncol = 1) +
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
# 
# plot.cs.den <- ggpubr::ggarrange(plot.cs.den.1, plot.cs.den.2,
#                                  nrow = 1,
#                                  common.legend = T, legend = 'bottom')
# ggsave(sprintf("%s/%s/COLONY_SIZE_DENSITY.jpg",fig_path, expt.name), plot.cs.den,
#        height = 200, width = 150, units = 'mm',
#        dpi = 600)  

##### EFFECT SIZE OF CS DIFFERENCES
# plot.rcs.change <- ggplot(anova.res[anova.res$strain2 == 'FY4',],
#        aes(x = time, y = relative_cs, col = strain1)) +
#   geom_point() + geom_line(aes(group = strain1)) +
#   labs(x = 'Time Point', y = 'Relative Colony Size') +
#   scale_color_discrete(name = 'Strains') +
#   facet_wrap(.~attempt*condition, nrow = 1) +
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
# ggsave(sprintf("%s/%s/COLONY_SIZE_RELATIVE_COLONY_SIZE_CHANGE.jpg",fig_path, expt.name), plot.rcs.change,
#        height = 100, width = 150, units = 'mm',
#        dpi = 600)  
# 
# 
# plot.effsize.res <- ggplot(anova.res[anova.res$strain2 == 'FY4',],
#                               aes(x = strain1, y = effect_size)) +
#   geom_point() +
#   scale_x_discrete(breaks = unique(anova.res$strain1), drop = F) +
#   labs(x = 'Strain', y = 'Effect Size') +
#   facet_wrap(.~attempt*condition*time, ncol = 2) +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%s/%s/COLONY_SIZE_EFFECT_SIZE.jpg",fig_path, expt.name), plot.effsize.res,
#        height = 200, width = 150, units = 'mm',
#         dpi = 600)
# 
##### RELATIVE COLONY SIZE ANALYSIS (DOESN'T INCLUDE BIOLOGICAL REPLICATES)
# plot.rcs.den.1 <- data[data$time == 't_0',] %>%
#   ggplot(aes(x = relative_cs, y = orf_name, fill = bio_rep, col = attempt)) +
#   geom_density_ridges(scale = 3, alpha = 0.7, size = 1.2,
#                       vline_size = 0.2, vline_color = "black") +
#   labs(x = 'Relative Colony Size',
#        y = 'Strains') +
#   scale_fill_discrete(name = 'Biological Replicate') +
#   scale_color_discrete(name = 'Attempt') +
#   facet_wrap(.~condition*time, ncol = 1) +
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

plot.rcs.den.2 <- data[data$time == 't_final',] %>%
  ggplot(aes(x = relative_cs, y = attempt, fill = bio_rep)) +
  geom_density_ridges(quantile_lines = TRUE,
                      scale = 3, alpha = 0.7, size = 0.5,
                      vline_size = 0.2, vline_color = "black") +
  labs(x = 'Relative Colony Size',
       y = 'Attempt') +
  scale_fill_discrete(name = 'Biological Replicate') +
  facet_wrap(.~condition*orf_name*time, nrow = 2) +
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
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_DENSITY.jpg",fig_path, expt.name), plot.rcs.den.2,
       height = 100, width = 400, units = 'mm',
       dpi = 600)
# plot.rcs.den <- ggpubr::ggarrange(plot.rcs.den.1, plot.rcs.den.2,
#                                  nrow = 1,
#                                  common.legend = T, legend = 'bottom')
# ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_DENSITY.jpg",fig_path, expt.name), plot.rcs.den,
#        height = 200, width = 150, units = 'mm',
#        dpi = 600) 

plot.rcs.den.anova <- data[data$time == 't_final',] %>%
  ggplot(aes(x = relative_cs, y = orf_name)) +
  geom_density_ridges(quantile_lines = TRUE,
                      aes(fill = bio_rep),
                      scale = 2, alpha = 0.7, size = 0.2,
                      vline_size = 0.2, vline_color = "black") +
  geom_point(data = anova.res[anova.res$strain2 == 'FY4' & 
                                anova.res$time == 't_final' & anova.res$rcs_between < 0.05,],
             aes(y = strain1, x = 1.5), shape = 8, size = 1, col = 'red') +
  labs(x = 'Relative Colony Size',
       y = 'Strain') +
  scale_fill_discrete(name = 'Biological Replicate') +
  facet_wrap(.~attempt*condition*time, nrow = 1) +
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
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_DENSITY_ANOVA.jpg",fig_path, expt.name), plot.rcs.den.anova,
       height = 80, width = 200, units = 'mm',
       dpi = 600)

plot.rcs.anova.res <- ggplot(data[data$time == 't_final',],
       aes(x = orf_name, y = relative_cs)) +
  # geom_boxplot(aes(fill = orf_name), outlier.shape = NA) +
  geom_violin(aes(fill = bio_rep), lwd = 0.1) +
  geom_point(data = anova.res[anova.res$strain2 == 'FY4' & 
                                anova.res$time == 't_final' & anova.res$rcs_between < 0.05,],
             aes(x = strain1, y = 1.5), shape = 8, size = 1, col = 'red') +
  labs(x = 'Strains', y = 'Relative Colony Size') +
  scale_fill_discrete(name = 'Biological Replicate') +
  facet_wrap(.~attempt*condition*time, nrow = 1) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_VIOLIN_ANOVA.jpg",fig_path, expt.name),
       plot.rcs.anova.res,
       height = 80, width = 200, units = 'mm',
       dpi = 600)

plot.rcs.anova.res2 <- ggplot(data[data$time == 't_final',],
                             aes(x = attempt, y = relative_cs)) +
  # geom_boxplot(aes(fill = orf_name), outlier.shape = NA) +
  geom_violin(aes(fill = bio_rep), lwd = 0.1) +
  labs(x = 'Strains', y = 'Relative Colony Size') +
  scale_fill_discrete(name = 'Biological Replicate') +
  stat_compare_means(method = 'anova', size = 1.2, label.y = 1.5) +
  facet_wrap(.~condition*orf_name*time, nrow = 2) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_VIOLIN_ANOVA2.jpg",fig_path, expt.name),
       plot.rcs.anova.res2,
       height = 100, width = 240, units = 'mm',
       dpi = 600)

plot.rcs.anova.res3 <- data[data$time == 't_final',] %>%
  group_by(attempt, condition, time, orf_name, bio_rep) %>%
  summarise(average = mean(average, na.rm = T), relative_cs = mean(relative_cs, na.rm = T)) %>%
  data.frame() %>%
  ggplot(aes(x = attempt, y = relative_cs)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col = bio_rep, shape = 'attempt')) +
  labs(x = 'Strains', y = 'Relative Colony Size') +
  scale_y_continuous(minor_breaks = seq(-2,2,0.02)) +
  scale_color_discrete(name = 'Biological Replicate') +
  scale_shape_discrete(name = 'Experimental Replicate') +
  stat_compare_means(method = 'kruskal', size = 1.2, label.y = 1.5) +
  facet_wrap(.~condition*orf_name*time, nrow = 2) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX_ANOVA.jpg",fig_path, expt.name),
       plot.rcs.anova.res3,
       height = 100, width = 240, units = 'mm',
       dpi = 600)

plot.rcs.ttest.res <- data[data$time == 't_final',] %>%
  group_by(attempt, condition, time, orf_name, bio_rep) %>%
  summarise(average = mean(average, na.rm = T), relative_cs = mean(relative_cs, na.rm = T)) %>%
  data.frame() %>%
  ggplot(aes(x = orf_name, y = relative_cs)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col = bio_rep, shape = attempt), size = 1) +
  labs(x = 'Strains', y = 'Relative Colony Size') +
  scale_color_discrete(name = 'Biological Replicate') +
  scale_shape_discrete(name = 'Experimental Replicate') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))

ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX_TTEST.jpg",fig_path, expt.name),
       plot.rcs.ttest.res + stat_compare_means(method = 't.test', ref.group = 'FY4',
                                               size = 1.2, label.y = 1.5, label = 'p.signif') +
         facet_wrap(.~attempt*condition*time, nrow = 1),
       height = 80, width = 200, units = 'mm',
       dpi = 600)
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX_TTEST2.jpg",fig_path, expt.name),
       plot.rcs.ttest.res + stat_compare_means(method = 't.test', ref.group = 'FY4',
                                               size = 1.2, label.y = 1.5, label = 'p.signif') +
         facet_wrap(.~condition*time, nrow = 1),
       height = 80, width = 200, units = 'mm',
       dpi = 600)

##### OVERALL RESULTS (COMBINING EVERYTHING)
anova.res2 <- NULL
for (c in unique(data$condition)) {
  for (t in unique(data$time[data$condition == c])) {
    for (i in seq(1,dim(strain.pairs)[1])) {
      o1 <- strain.pairs[i,1]
      o2 <- strain.pairs[i,2]

      res.aov <- data[data$time == t & data$condition == c & data$attempt != 'pilot' & data$orf_name %in% c(o1, o2),] %>%
        anova_test(average ~ orf_name * bio_rep)
      res.rcs.aov <- data[data$time == t & data$condition == c & data$attempt != 'pilot' & data$orf_name %in% c(o1, o2),] %>%
        anova_test(relative_cs ~ orf_name * bio_rep)

      res.rcs.kw <- compare_means(data = data[data$time == t & data$condition == c & data$attempt != 'pilot' & data$orf_name %in% c(o1, o2),] %>%
                                    group_by(time, condition, attempt, orf_name, bio_rep) %>%
                                    summarize(average = median(average, na.rm = T), relative_cs = median(relative_cs, na.rm = T),
                                              .groups = 'keep') %>%
                                    data.frame(), method = 'kruskal', formula = relative_cs ~ orf_name) %>%
        data.frame()
      
      cs1 <- data$average[data$time == t & data$condition == c & data$attempt != 'pilot' & data$orf_name == o1]
      cs2 <- data$average[data$time == t & data$condition == c & data$attempt != 'pilot' & data$orf_name == o2]
      effsize.temp <- wilcox.test(cs1, cs2, conf.int = T)
      emp_effsize <- (median(cs1, na.rm = T) - median(cs2, na.rm = T))/median(cs2, na.rm = T)
      delta.temp <- cliff.delta(cs1, cs2)
      relative_cs <- median(cs1, na.rm = T)/median(cs2, na.rm = T)

      anova.res2 <- rbind(anova.res2, cbind(c, t, o1, o2,
                                            res.aov$p[res.aov$Effect == 'orf_name'],
                                            res.aov$p[res.aov$Effect == 'orf_name:bio_rep'],
                                            res.rcs.aov$p[res.aov$Effect == 'orf_name'],
                                            res.rcs.aov$p[res.aov$Effect == 'orf_name:bio_rep'],
                                            effsize.temp$estimate,
                                            delta.temp$estimate, delta.temp$magnitude,
                                            emp_effsize, relative_cs,
                                            res.rcs.kw$p))
    }
  }
}
colnames(anova.res2) <- c('condition','time','strain1','strain2','between','within','rcs_between','rcs_within',
                         'wilcox_estimate','cliffs_estimate','cliffs_magnitude','effect_size','relative_cs','rcs_kw')
anova.res2 <- data.frame(anova.res2, stringsAsFactors = F)
# anova.res2$hours <- as.numeric(anova.res2$hours)
anova.res2$between <- as.numeric(anova.res2$between)
anova.res2$within <- as.numeric(anova.res2$within)
anova.res2$rcs_between <- as.numeric(anova.res2$rcs_between)
anova.res2$rcs_within <- as.numeric(anova.res2$rcs_within)
anova.res2$wilcox_estimate <- as.numeric(anova.res2$wilcox_estimate)
anova.res2$cliffs_estimate <- as.numeric(anova.res2$cliffs_estimate)
anova.res2$effect_size <- as.numeric(anova.res2$effect_size)
anova.res2$relative_cs <- as.numeric(anova.res2$relative_cs)
anova.res2$between <- p.adjust(anova.res2$between, method = 'BH')
anova.res2$within <- p.adjust(anova.res2$within, method = 'BH')
anova.res2$rcs_between <- p.adjust(anova.res2$rcs_between, method = 'BH')
anova.res2$rcs_within <- p.adjust(anova.res2$rcs_within, method = 'BH')
anova.res2$rcs_kw <- p.adjust(anova.res2$rcs_kw, method = 'BH')

anova.res2$strain1 <- factor(anova.res2$strain1, levels = strain.levels)
anova.res2$strain2 <- factor(anova.res2$strain2, levels = strain.levels)
anova.res2$condition <- factor(anova.res2$condition, levels = condition.levels)
anova.res2$cliffs_magnitude <- factor(anova.res2$cliffs_magnitude, levels = c('negligible','small','medium','large'))

# anova.res2$time[anova.res2$hours < 10] <- 't_0'
# anova.res2$time[anova.res2$hours > 10] <- 't_final'

save(anova.res2, file = sprintf('%s/%s/anova_results2_nopilot.RData', out_path, expt.name))

ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX_KW.jpg",fig_path, expt.name),
       plot.rcs.ttest.res +
         geom_point(data = anova.res2[anova.res2$strain2 == 'FY4' &
                                        anova.res2$time == 't_final' & anova.res2$rcs_kw < 0.05,],
                    aes(x = strain1, y = 1.5), shape = 8, size = 1, col = 'red') +
         geom_text(data = anova.res2[anova.res2$strain2 == 'FY4' & anova.res2$time == 't_final',],
                   aes(x = strain1, y = 1.4, label = sprintf('%0.2f%%',effect_size*100)),
                   size = 1.5) +
         facet_wrap(.~condition*time, nrow = 1),
       height = 80, width = 200, units = 'mm',
       dpi = 600)

plot.rcs.anova.res4 <- ggplot(data[data$time == 't_final',],
                             aes(x = orf_name, y = relative_cs)) +
  geom_boxplot(aes(fill = orf_name), outlier.shape = NA) +
  geom_violin(fill = 'transparent') +
  geom_point(data = anova.res2[anova.res2$strain2 == 'FY4' & anova.res2$time == 't_final' & anova.res2$rcs_between < 0.05,],
             aes(x = strain1, y = 1.5), shape = 8, size = 1, col = 'red') +
  labs(x = 'Strains', y = 'Relative Colony Size') +
  scale_fill_discrete(guide = F) +
  facet_wrap(.~condition*time, nrow = 1) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_ALLCOMBINED_ANOVA.jpg",fig_path, expt.name),
       plot.rcs.anova.res4,
       height = 100, width = 120, units = 'mm',
       dpi = 600)

##### ALL x ALL EFFECT SIZE (INCLUDES BIOLOGICAL REPLICATES)
# strain.pairs <- permutations(n = length(unique(data$orf_name)), r = 2,
#                              v = as.character(unique(data$orf_name)),
#                              repeats.allowed = T) %>%
#   data.frame(stringsAsFactors = F)
# 
# effsize.res <- NULL
# for (a in unique(data$attempt)) {
#   for (c in unique(data$condition[data$attempt == a])) {
#     for (h in unique(data$hours[data$attempt == a & data$condition == c])) {
#       for (i in seq(1,dim(strain.pairs)[1])) {
#         o1 <- strain.pairs[i,1]
#         o2 <- strain.pairs[i,2]
#         for (b1 in unique(data$bio_rep[data$attempt == a & data$hours == h & data$condition == c & data$orf_name == o1])) {
#           for (b2 in unique(data$bio_rep[data$attempt == a & data$hours == h & data$condition == c & data$orf_name == o2])) {
#             cs1 <- data$average[data$attempt == a & data$hours == h & data$condition == c & data$orf_name == o1 & data$bio_rep == b1]
#             cs2 <- data$average[data$attempt == a & data$hours == h & data$condition == c & data$orf_name == o2 & data$bio_rep == b2]
# 
#             effsize.temp <- wilcox.test(cs1, cs2, conf.int = T)
#             emp_effsize <- (median(cs1, na.rm = T) - median(cs2, na.rm = T))/median(cs2, na.rm = T)
#             delta.temp <- cliff.delta(cs1, cs2)
#             relative_cs <- median(cs1, na.rm = T)/median(cs2, na.rm = T)
# 
#             effsize.res <- rbind(effsize.res, cbind(a, c, h, o1, b1, o2, b2,
#                                                     effsize.temp$estimate[[1]],
#                                                     delta.temp$estimate, as.character(delta.temp$magnitude),
#                                                     emp_effsize, relative_cs))
#           }
#         }
#       }
#     }
#   }
# }
# colnames(effsize.res) <- c('attempt','condition','hours','strain1','bio_rep1','strain2','bio_rep2',
#                            'wilcox_estimate','cliffs_estimate','cliffs_magnitude','effect_size','relative_cs')
# effsize.res <- data.frame(effsize.res, stringsAsFactors = F)
# effsize.res$hours <- as.numeric(effsize.res$hours)
# effsize.res$wilcox_estimate <- as.numeric(effsize.res$wilcox_estimate)
# effsize.res$cliffs_estimate <- as.numeric(effsize.res$cliffs_estimate)
# effsize.res$effect_size <- as.numeric(effsize.res$effect_size)
# effsize.res$relative_cs <- as.numeric(effsize.res$relative_cs)
# 
# effsize.res$attempt <- factor(effsize.res$attempt, levels = attempt.levels)
# effsize.res$strain1 <- factor(effsize.res$strain1,
#                             levels = strain.levels)
# effsize.res$strain2 <- factor(effsize.res$strain2,
#                             levels = strain.levels)
# effsize.res$condition <- factor(effsize.res$condition, levels = condition.levels)
# effsize.res$cliffs_magnitude <- factor(effsize.res$cliffs_magnitude, levels = c('negligible','small','medium','large'))
# effsize.res <- effsize.res[!(effsize.res$strain1 == effsize.res$strain2 & effsize.res$bio_rep1 == effsize.res$bio_rep2),]
# 
# write.csv(effsize.res, file = sprintf('%s/%s/COLONY_SIZE_EFFECT_SIZE_RESULTS.csv', res_path, expt.name))
# save(effsize.res, file = sprintf("%s/%s/effectsizes.RData", out_path, expt.name))

load(sprintf("%s/%s/effectsizes.RData", out_path, expt.name))

effsize.res$time[effsize.res$hours == 0] <- 't_0'
effsize.res$time[effsize.res$hours > 0] <- 't_final'
head(effsize.res)

# plot.effsize.all.res <- effsize.res[!(effsize.res$condition == 'YPDA' & effsize.res$hours > 100),] %>%
#   filter(strain2 == 'FY4') %>%
#   ggplot(aes(x = strain1, y = effect_size)) +
#   geom_boxplot(fill = 'transparent', outlier.shape = NA) +
#   geom_jitter(aes(col = attempt), size = 0.2) +
#   stat_compare_means(method = 't.test', ref.group = 'FY4',
#                      label = 'p.signif', label.y = 0.5,
#                      size = 1.2) +
#   scale_x_discrete(breaks = unique(effsize.res$strain1), drop = F) +
#   scale_y_continuous(minor_breaks = seq(-1,2,0.02)) +
#   labs(x = 'Strains', y = 'Effect Size') +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm"))) +
#   coord_cartesian(ylim = c(0,1.5))
# 
# ggsave(sprintf("%s/%s/COLONY_SIZE_ALLbyALL_EFFECT_SIZE.jpg",fig_path, expt.name),
#        plot.effsize.all.res + facet_wrap(.~attempt*condition*time, ncol = 2),
#        height = 200, width = 150, units = 'mm',
#        dpi = 600)
# ggsave(sprintf("%s/%s/COLONY_SIZE_ALLbyALL_EFFECT_SIZE_2.jpg",fig_path, expt.name),
#        plot.effsize.all.res + facet_wrap(.~condition*time, nrow = 2),
#        height = 80, width = 200, units = 'mm',
#        dpi = 600) 

plot.rcs.all.res <- effsize.res[!(effsize.res$condition == 'YPDA' & effsize.res$hours > 100) &
                                  effsize.res$time == 't_final',] %>%
  filter(strain2 == 'FY4') %>%
  ggplot(aes(x = strain1, y = relative_cs)) +
  geom_boxplot(fill = 'transparent', outlier.shape = NA) +
  geom_jitter(aes(col = attempt), size = 0.2) +
  stat_compare_means(method = 't.test', ref.group = 'FY4',
                     label = 'p.signif', label.y = 1.4,
                     size = 1.2) +
  scale_x_discrete(breaks = unique(effsize.res$strain1), drop = F) +
  scale_y_continuous(minor_breaks = seq(-2,2,0.02)) +
  labs(x = 'Strains', y = 'Relative Colony Size') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))

ggsave(sprintf("%s/%s/COLONY_SIZE_ALLbyALL_RELATIVE_CS.jpg",fig_path, expt.name),
       plot.rcs.all.res + facet_wrap(.~attempt*condition*time, nrow = 1),
       height = 80, width = 200, units = 'mm',
       dpi = 600)
ggsave(sprintf("%s/%s/COLONY_SIZE_ALLbyALL_RELATIVE_CS_2.jpg",fig_path, expt.name),
       plot.rcs.all.res + facet_wrap(.~condition*time, nrow = 1),
       height = 80, width = 100, units = 'mm',
       dpi = 600)


effsize.res[!(effsize.res$condition == 'YPDA' & effsize.res$hours > 100),] %>%
  filter(strain2 == 'FY4') %>%
  group_by(attempt, condition, hours, strain1) %>%
  count()


plot.effsize.all.res.media <- effsize.res[!(effsize.res$condition == 'YPDA' & effsize.res$hours > 100) &
                                            effsize.res$hours != 0 &
                                            effsize.res$attempt != 'pilot',] %>%
  filter(strain2 == 'FY4') %>%
  ggplot(aes(x = strain1, y = effect_size, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  # stat_summary(fun = mean, group = 'condition', geom = "line") +
  stat_compare_means(method = 't.test', label = 'p.signif',
                     size = 1.2, label.y = 0.5) +
  facet_wrap(.~time) +
  scale_x_discrete(breaks = unique(effsize.res$strain1), drop = F) +
  scale_y_continuous(minor_breaks = seq(-2,1,0.02)) +
  labs(x = 'Strains', y = 'Effect Size') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(-1,0.5))

ggsave(sprintf("%s/%s/COLONY_SIZE_ALLbyALL_EFFECT_SIZE_3.jpg",fig_path, expt.name),
       plot.effsize.all.res.media,
       height = 100, width = 150, units = 'mm',
       dpi = 600) 


plot.rcs.all.res.media <- effsize.res[!(effsize.res$condition == 'YPDA' & effsize.res$hours > 100) &
                                            effsize.res$hours != 0 &
                                            effsize.res$attempt != 'pilot',] %>%
  filter(strain2 == 'FY4') %>%
  ggplot(aes(x = strain1, y = relative_cs, fill = condition)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_violin(alpha = 0.5) +
  # stat_summary(fun = mean, group = 'condition', geom = "line") +
  stat_compare_means(method = 't.test', label = 'p.signif',
                     size = 1.2, label.y = 1.2) +
  facet_wrap(.~time) +
  scale_x_discrete(breaks = unique(effsize.res$strain1), drop = F) +
  scale_y_continuous(minor_breaks = seq(-2,2,0.01)) +
  scale_fill_discrete(name = 'Media Type') +
  labs(x = 'Strains', y = 'Relative Colony Size') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.2))
ggsave(sprintf("%s/%s/COLONY_SIZE_ALLbyALL_RELATIVE_CS_3.jpg",fig_path, expt.name),
       plot.rcs.all.res.media,
       height = 100, width = 150, units = 'mm',
       dpi = 600)


head(effsize.res)
plot.rcs.t0tf.cor <- merge(effsize.res[effsize.res$time == 't_0',], effsize.res[!(effsize.res$condition == 'YPDA' & effsize.res$hours > 100) & 
                                                                                  effsize.res$time == 't_final',],
      by = c('attempt','condition','strain1','bio_rep1','strain2','bio_rep2'),
      suffixes = c('.t_0','.t_final')) %>%
  filter(strain2 == 'FY4') %>%
  ggplot(aes(x = relative_cs.t_final, y = relative_cs.t_0)) +
  geom_point(size = 1) +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'spearman', size = 1.2, label.x = 0.4) +
  labs(x = 'T_Final Relative Colony Size',
       y = 'T_0 Relative Colony Size') +
  # facet_wrap(.~condition, ncol = 12) +
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
ggsave(sprintf("%s/%s/COLONY_SIZE_T0TF_CORR2.jpg",fig_path, expt.name), plot.rcs.t0tf.cor,
       height = 60, width = 60, units = 'mm',
       dpi = 600)





