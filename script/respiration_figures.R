##### RESPIRATION EXPERIMENT - FIGURES
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 05/03/2021 

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
library(plotly)
library(scales)
library(reshape2)
library(locfit)
library(growthcurver)
library(rstatix)
library(gtools)

load("~/R/Projects/methionine/data/methionine_respiration_colonysizes.RData")
load("~/R/Projects/methionine/data/methionine_respiration_growthcurve.RData")
fig_path <- "~/R/Projects/methionine/figures"
res_path <- "~/R/Projects/methionine/results"

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### ALL GROWTH CURVES
data <- data[data$orf_name != 'BOR',]
data$orf_name <- factor(data$orf_name,
                          levels = c('FY4', 'FY4_pet', 'FY4-met15del', 'FY4-met15del_pet', 'BY4742', 'BY4742_pet', 'BY4741', 'BY4741_pet'))
data$condition <- factor(data$condition,
                           levels = c('SD+MET+Glucose','SD-MET+Glucose','SD+MET+EtOH','SD-MET+EtOH','SD+MET-URA+Glucose'))

all.gc <- data %>%
  ggplot(aes(x = hours, y = average, col = bio_rep)) +
  geom_point() +
  geom_smooth(method = 'loess') +
  scale_x_continuous(breaks = unique(data$hours), minor_breaks = NULL) +
  scale_color_discrete(name = 'Biological Replicate') +
  labs(y = 'Colony Size (pixels)',
       x = 'Time (hours)') +
  facet_wrap(.~condition*orf_name,
             ncol = 8) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(size = txt - 2, angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/RESPIRATION_ALL_GROWTH_CURVES_PROBLEM.jpg",fig_path), all.gc,
       height = 300, width = 600, units = 'mm',
       dpi = 600)

## CLEANED OUT PROBLEM HOURS
all.gc <- data[!(data$hours %in% c(46.23, 87.35, 25.88, 16.75, 54.36, 37.13)),] %>%
  ggplot(aes(x = hours, y = average, col = bio_rep)) +
  geom_smooth(method = 'loess') +
  scale_color_discrete(name = 'Biological Replicate') +
  labs(y = 'Colony Size (pixels)',
       x = 'Time (hours)') +
  facet_wrap(.~condition*orf_name,
             ncol = 8) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(size = txt - 2, angle = 90),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/RESPIRATION_ALL_GROWTH_CURVES.jpg",fig_path), all.gc,
       height = 300, width = 600, units = 'mm',
       dpi = 600)

##### COLONY SIZE COMPARISONS
data[data$hours %in% c(174.43, 327.88),] %>%
  group_by(condition, hours, orf_name, bio_rep) %>%
  count()

plot.cs.vio <- data[data$hours %in% c(174.43, 327.88),] %>%
  ggplot(aes(x = orf_name, y = average)) +
  geom_violin(aes(fill = bio_rep),lwd = 0.1) +
  # geom_jitter() +
  scale_fill_discrete(name = 'Biological Replicate') +
  facet_wrap(.~condition*hours, nrow = 5) +
  labs(x = 'Strain',
       y = 'Colony Size (pixels)') +
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
ggsave(sprintf("%s/RESPIRATION_COLONY_SIZE_VIOLIN.jpg",fig_path), plot.cs.vio,
       height = 200, width = 150, units = 'mm',
       dpi = 600)  

##### ANOVA
strain.pairs <- permutations(n = length(unique(data$orf_name)), r = 2,
                             v = as.character(unique(data$orf_name)),
                             repeats.allowed = F) %>%
  data.frame(stringsAsFactors = F)

hrs <- c(174.43, 327.88)
anova.res <- NULL
for (h in hrs) {
  for (c in unique(data$condition)) {
    for (i in seq(1,dim(strain.pairs)[1])) {
      o1 <- strain.pairs[i,1]
      o2 <- strain.pairs[i,2]
      res.aov <- data[data$hours == h & data$condition == c & data$orf_name %in% c(o1, o2),] %>%
        anova_test(average ~ orf_name * bio_rep)
      anova.res <- rbind(anova.res, cbind(c, h, o1, o2,
                                          res.aov$p[res.aov$Effect == 'orf_name'],
                                          res.aov$p[res.aov$Effect == 'orf_name:bio_rep']))
    }
  }
}
colnames(anova.res) <- c('condition','hours','strain1','strain2','between','within')
anova.res <- data.frame(anova.res, stringsAsFactors = F)
anova.res$between <- as.numeric(anova.res$between)
anova.res$within <- as.numeric(anova.res$within)
# anova.res$between <- p.adjust(anova.res$between, method = 'BH')
# anova.res$within <- p.adjust(anova.res$within, method = 'BH')

anova.res$strain1 <- factor(anova.res$strain1,
                            levels = c('FY4', 'FY4_pet', 'FY4-met15del', 'FY4-met15del_pet', 'BY4742', 'BY4742_pet', 'BY4741', 'BY4741_pet'))
anova.res$strain2 <- factor(anova.res$strain2,
                            levels = c('FY4', 'FY4_pet', 'FY4-met15del', 'FY4-met15del_pet', 'BY4742', 'BY4742_pet', 'BY4741', 'BY4741_pet'))
anova.res$condition <- factor(anova.res$condition,
                              levels = c('SD+MET+Glucose','SD-MET+Glucose','SD+MET+EtOH','SD-MET+EtOH','SD+MET-URA+Glucose'))

write.csv(anova.res, file = sprintf('%s/RESPIRATION_COLONY_SIZE_ANOVA_RESULTS.csv', res_path))

plot.anova.res <- ggplot(anova.res) +
  geom_tile(aes(x = strain1, y = strain2, fill = between <= 0.05),
            col = 'black') +
  facet_wrap(.~condition*hours, nrow = 5) +
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
ggsave(sprintf("%s/RESPIRATION_COLONY_SIZE_ANOVA.jpg",fig_path), plot.anova.res,
       height = 200, width = 150, units = 'mm',
       dpi = 600)  


##### GROWTH CURVE ANALYSIS
data.pred <- data.pred[,colSums(log(data.pred[,-1]) < 5) <= 300]
# data.pred[,-1] <- log(data.pred[,-1])
gc.res <- SummarizeGrowthByPlate(data.pred)
temp <- str_split(gc.res$sample, ',', simplify = T)
colnames(temp) <- c('orf_name','condition','bio_rep')
gc.res <- cbind(temp, gc.res)
head(gc.res)

temp.count <- gc.res %>%
  group_by(orf_name, condition) %>%
  count() %>% 
  data.frame()
gc.res <- merge(gc.res, temp.count, by = c('orf_name', 'condition'))

gc.res$orf_name <- as.character(gc.res$orf_name)
gc.res$condition <- as.character(gc.res$condition)
gc.res$orf_name <- factor(gc.res$orf_name,
                          levels = c('FY4', 'FY4_pet', 'FY4-met15del', 'FY4-met15del_pet', 'BY4742', 'BY4742_pet', 'BY4741', 'BY4741_pet'))
gc.res$condition <- factor(gc.res$condition,
                           levels = c('SD+MET+Glucose','SD-MET+Glucose','SD+MET+EtOH','SD-MET+EtOH','SD+MET-URA+Glucose'))

write.csv(gc.res, file = sprintf('%s/RESPIRATION_GROWTH_CURVE_RESULTS.csv', res_path))

plot.gc.auc <- gc.res[gc.res$n > 1,] %>%
  ggplot(aes(x = orf_name, y = auc_e)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col = bio_rep)) +
  # stat_compare_means(method = 't.test', ref.group = 'FY4', label = "p.signif") +
  scale_color_discrete(name = 'Biological Replicate') +
  facet_wrap(.~condition, nrow = 5) +
  labs(x = 'Strain',
       y = 'AUC',
       title = 'Area Under the Curve (AUC)') +
  # coord_flip() +
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

plot.gc.tgen <- gc.res[gc.res$n > 1,] %>%
  ggplot(aes(x = orf_name, y = t_gen)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col = bio_rep)) +
  scale_color_discrete(name = 'Biological Replicate') +
  facet_wrap(.~condition, nrow = 5) +
  labs(x = 'Strain',
       y = 'DT (hours)',
       title = 'Doubling Time (DT)') +
  # coord_flip() +
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

plot.gc.res <- ggpubr::ggarrange(plot.gc.tgen, plot.gc.auc,
          ncol = 2,
          common.legend = T, legend = 'bottom')

ggsave(sprintf("%s/RESPIRATION_GROWTH_CURVES_RESULTS.jpg",fig_path), plot.gc.res,
       height = 200, width = 150, units = 'mm',
       dpi = 600)

## AUC AND DT T.TESTS
strain.pairs <- permutations(n = length(unique(data$orf_name)), r = 2,
                             v = as.character(unique(data$orf_name)),
                             repeats.allowed = F) %>%
  data.frame(stringsAsFactors = F)

gc.ttest.res <- NULL
for (c in unique(gc.res$condition)) {
  for (i in seq(1,dim(strain.pairs)[1])) {
    o1 <- strain.pairs[i,1]
    o2 <- strain.pairs[i,2]
    
    if (length(gc.res$t_gen[gc.res$condition == c & gc.res$orf_name == o1]) > 1 &
        length(gc.res$t_gen[gc.res$condition == c & gc.res$orf_name == o2]) > 1) {
      res.ttest.t_gen <- t.test(gc.res$t_gen[gc.res$condition == c & gc.res$orf_name == o1],
                          gc.res$t_gen[gc.res$condition == c & gc.res$orf_name == o2])
      res.ttest.auc_e <- t.test(gc.res$auc_e[gc.res$condition == c & gc.res$orf_name == o1],
                                gc.res$auc_e[gc.res$condition == c & gc.res$orf_name == o2])
      gc.ttest.res <- rbind(gc.ttest.res, cbind(c, o1, o2,
                                                res.ttest.t_gen$p.value,
                                                res.ttest.auc_e$p.value))

    }
  }
}

colnames(gc.ttest.res) <- c('condition','strain1','strain2','dt_p','auc_p')
gc.ttest.res <- data.frame(gc.ttest.res, stringsAsFactors = F)
gc.ttest.res$dt_p <- as.numeric(gc.ttest.res$dt_p)
gc.ttest.res$auc_p <- as.numeric(gc.ttest.res$auc_p)

gc.ttest.res$strain1 <- factor(gc.ttest.res$strain1,
                            levels = c('FY4', 'FY4_pet', 'FY4-met15del', 'FY4-met15del_pet', 'BY4742', 'BY4742_pet', 'BY4741', 'BY4741_pet'))
gc.ttest.res$strain2 <- factor(gc.ttest.res$strain2,
                            levels = c('FY4', 'FY4_pet', 'FY4-met15del', 'FY4-met15del_pet', 'BY4742', 'BY4742_pet', 'BY4741', 'BY4741_pet'))
gc.ttest.res$condition <- factor(gc.ttest.res$condition,
                              levels = c('SD+MET+Glucose','SD-MET+Glucose','SD+MET+EtOH','SD-MET+EtOH','SD+MET-URA+Glucose'))

write.csv(gc.ttest.res, file = sprintf('%s/RESPIRATION_GROWTH_CURVE_TTEST_RESULTS.csv', res_path))

plot.ttest.res.dt <- ggplot(gc.ttest.res) +
  geom_tile(aes(x = strain1, y = strain2, fill = dt_p <= 0.05),
            col = 'black') +
  labs(title = 'Doubling Time (DT)') +
  facet_wrap(.~condition, nrow = 5) +
  scale_fill_manual(name = 'Significant T-Test',
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
  
plot.ttest.res.auc <- ggplot(gc.ttest.res) +
  geom_tile(aes(x = strain1, y = strain2, fill = auc_p <= 0.05),
            col = 'black') +
  labs(title = 'Area Under the Curve (AUC)') +
  facet_wrap(.~condition, nrow = 5) +
  scale_fill_manual(name = 'Significant T-Test',
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

plot.gc.ttest.res <- ggpubr::ggarrange(plot.ttest.res.dt, plot.ttest.res.auc,
                                 ncol = 2,
                                 common.legend = T, legend = 'bottom')
ggsave(sprintf("%s/RESPIRATION_GROWTH_CURVES_TTEST.jpg",fig_path), plot.gc.ttest.res,
       height = 200, width = 150, units = 'mm',
       dpi = 600)

