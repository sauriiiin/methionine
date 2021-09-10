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
library(locfit)
library(growthrates)

out_path <- "~/R/Projects/methionine/data"
fig_path <- "~/R/Projects/methionine/figures"
res_path <- "~/R/Projects/methionine/results"

expt.name <- "respiration"

load(sprintf("%s/%s/colonysizes.RData", out_path, expt.name))

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### LEVELS
strain.levels = c('FY4', 'FY4_pet', 'FY4-met15del', 'FY4-met15del_pet', 'BY4742', 'BY4742_pet', 'BY4741', 'BY4741_pet')
condition.levels = c('SD+Met-Cys+Glu','SD-Met-Cys+Glu','SD+Met-Cys+EtOH','SD-Met-Cys+EtOH','SD+Met-Ura+Glu')

##### ALL GROWTH CURVES
data <- data[data$orf_name != 'BOR',]
data$orf_name <- factor(data$orf_name, levels = strain.levels)
data$condition <- factor(data$condition, levels = condition.levels)

# for (rr in unique(data$expt_rep)) {
#   all.gc <- data %>%
#     filter(expt_rep == rr) %>% 
#     ggplot(aes(x = hours, y = average, col = bio_rep)) +
#     geom_point(size = 0.5) +
#     geom_smooth(method = 'loess') +
#     scale_x_continuous(breaks = unique(data$hours), minor_breaks = NULL) +
#     scale_color_discrete(name = 'Biological Replicate') +
#     labs(title = sprintf('Experimental Replicate #%s', rr),
#          y = 'Colony Size (pixels)',
#          x = 'Time (hours)') +
#     facet_wrap(.~condition*orf_name,
#                ncol = 8) +
#     theme_linedraw() +
#     theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#           axis.title = element_text(size = titles),
#           axis.text = element_text(size = txt),
#           axis.text.x = element_text(size = 2, angle = 90, hjust = 0.5, vjust = 0.5),
#           legend.title = element_text(size = titles),
#           legend.text = element_text(size = txt),
#           legend.position = 'bottom',
#           legend.key.size = unit(3, "mm"),
#           legend.box.spacing = unit(0.5,"mm"),
#           strip.text = element_text(size = txt,
#                                     face = 'bold',
#                                     margin = margin(0.1,0,0.1,0, "mm")))
#   ggsave(sprintf("%s/%s/ALL_GROWTH_CURVES_%s.jpg",fig_path, expt.name, rr), all.gc,
#          height = 200, width = 400, units = 'mm',
#          dpi = 600)
# }


# ## CLEANED OUT PROBLEM HOURS
# all.gc <- data[!(data$hours %in% c(46.23, 87.35, 25.88, 16.75, 54.36, 37.13)),] %>%
#   ggplot(aes(x = hours, y = average, col = bio_rep)) +
#   geom_smooth(method = 'loess') +
#   scale_color_discrete(name = 'Biological Replicate') +
#   labs(y = 'Colony Size (pixels)',
#        x = 'Time (hours)') +
#   facet_wrap(.~condition*orf_name,
#              ncol = 8) +
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
# ggsave(sprintf("%s/%s/ALL_GROWTH_CURVES.jpg",fig_path, expt.name), all.gc,
#        height = 200, width = 400, units = 'mm',
#        dpi = 600)
# 
# ##### COLONY SIZE COMPARISONS
# data[data$hours %in% c(174.43, 327.88),] %>%
#   group_by(condition, hours, orf_name, bio_rep) %>%
#   count() %>%
#   data.frame()
# 
# plot.cs.vio <- data[data$hours %in% c(174.43, 327.88),] %>%
#   ggplot(aes(x = orf_name, y = average)) +
#   geom_violin(aes(fill = bio_rep),lwd = 0.1) +
#   # geom_jitter() +
#   scale_fill_discrete(name = 'Biological Replicate') +
#   facet_wrap(.~condition*hours, nrow = 5) +
#   labs(x = 'Strain',
#        y = 'Colony Size (pixels)') +
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
# ggsave(sprintf("%s/%s/COLONY_SIZE_VIOLIN.jpg",fig_path, expt.name), plot.cs.vio,
#        height = 200, width = 150, units = 'mm',
#        dpi = 600)  
# 
# ##### ANOVA
# strain.pairs <- permutations(n = length(unique(data$orf_name)), r = 2,
#                              v = as.character(unique(data$orf_name)),
#                              repeats.allowed = F) %>%
#   data.frame(stringsAsFactors = F)
# 
# hrs <- c(174.43, 327.88)
# anova.res <- NULL
# for (h in hrs) {
#   for (c in unique(data$condition)) {
#     for (i in seq(1,dim(strain.pairs)[1])) {
#       o1 <- strain.pairs[i,1]
#       o2 <- strain.pairs[i,2]
#       res.aov <- data[data$hours == h & data$condition == c & data$orf_name %in% c(o1, o2),] %>%
#         anova_test(average ~ orf_name * bio_rep)
#       anova.res <- rbind(anova.res, cbind(c, h, o1, o2,
#                                           res.aov$p[res.aov$Effect == 'orf_name'],
#                                           res.aov$p[res.aov$Effect == 'orf_name:bio_rep']))
#     }
#   }
# }
# colnames(anova.res) <- c('condition','hours','strain1','strain2','between','within')
# anova.res <- data.frame(anova.res, stringsAsFactors = F)
# anova.res$between <- as.numeric(anova.res$between)
# anova.res$within <- as.numeric(anova.res$within)
# anova.res$between <- p.adjust(anova.res$between, method = 'BH')
# anova.res$within <- p.adjust(anova.res$within, method = 'BH')
# 
# anova.res$strain1 <- factor(anova.res$strain1, levels = strain.levels)
# anova.res$strain2 <- factor(anova.res$strain2, levels = strain.levels)
# anova.res$condition <- factor(anova.res$condition, levels = condition.levels)
# 
# write.csv(anova.res, file = sprintf('%s/%s/COLONY_SIZE_ANOVA_RESULTS.csv', res_path, expt.name))
# 
# plot.anova.res <- ggplot(anova.res) +
#   geom_tile(aes(x = strain1, y = strain2, fill = between <= 0.05),
#             col = 'black') +
#   facet_wrap(.~condition*hours, nrow = 5) +
#   scale_fill_manual(name = 'Significant Anova',
#                     breaks = c('TRUE', 'FALSE'),
#                     values = c('TRUE' = '#4CAF50', 'FALSE' = '#536DFE')) +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_blank(),
#         axis.text = element_text(size = txt),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm")),
#         panel.grid = element_blank())
# ggsave(sprintf("%s/%s/COLONY_SIZE_ANOVA.jpg",fig_path, expt.name), plot.anova.res,
#        height = 200, width = 150, units = 'mm',
#        dpi = 600)  
# 
# ## VIOLIN + ANOVA
# plot.cs.vio.anova <- data[data$hours == 327.88,] %>%
#   ggplot(aes(x = orf_name, y = average)) +
#   geom_violin(aes(fill = bio_rep),lwd = 0.1) +
#   geom_point(data = anova.res[anova.res$strain2 == 'FY4' &
#                                 anova.res$between <= 0.05 &
#                                 anova.res$hours == 327.88,],
#              aes(x = strain1, y = 5000), shape = 8, size = 1, col = 'red') +
#   scale_fill_discrete(name = 'Biological Replicate') +
#   facet_wrap(.~condition*hours, nrow = 1) +
#   labs(x = 'Strain',
#        y = 'Colony Size (pixels)') +
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
# ggsave(sprintf("%s/%s/COLONY_SIZE_VIOLIN_ANOVA.jpg",fig_path, expt.name), plot.cs.vio.anova,
#        height = 80, width = 200, units = 'mm',
#        dpi = 600) 
# 
# ## COLONY SIZE DATA SUMMARY
# data.sum <- data %>%
#   group_by(condition, hours, orf_name, bio_rep) %>%
#   summarise(average = median(average, na.rm = T)) %>%
#   data.frame()
# 
# plot.cs.sum.ttest <- data.sum %>%
#   filter(hours == 327.88) %>%
#   ggplot(aes(x = orf_name, y = average)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(aes(col = bio_rep), size = 1.2) +
#   stat_compare_means(method = 't.test', ref.group = 'FY4',
#                      size = 1, label = 'p.signif', label.y = 4000) +
#   scale_color_discrete(name = 'Biological Replicate') +
#   facet_wrap(.~condition*hours, nrow = 1) +
#   labs(x = 'Strain',
#        y = 'Median Colony Size (pixels)') +
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
# ggsave(sprintf("%s/%s/COLONY_SIZE_TTEST.jpg",fig_path, expt.name), plot.cs.sum.ttest,
#        height = 80, width = 200, units = 'mm',
#        dpi = 600) 


##### DATA FOR GROWTH CURVE ANALYSIS
data <- data[!(data$hours %in% c(46.23, 87.35, 25.88, 16.75, 54.36, 37.13)) &
               data$orf_name != 'BOR',]

data.gc <- data %>%
  group_by(expt_rep, condition, orf_name, bio_rep, pos, hours) %>%
  summarise(average = median(average, na.rm = T)) %>%
  data.frame()

data.pred <- NULL
col.names <- NULL
for (e in unique(data.gc$expt_rep)) {
  for (o in unique(data.gc$orf_name[data.gc$expt_rep == e])) {
    for (c in unique(data.gc$condition[data.gc$expt_rep == e & data.gc$orf_name == o])) {
      for (b in unique(data.gc$bio_rep[data.gc$expt_rep == e & data.gc$orf_name == o & data.gc$condition == c])) {
        for (p in unique(data.gc$pos[data.gc$expt_rep == e & data.gc$orf_name == o & data.gc$condition == c & data.gc$bio_rep == b])) {
          temp <- data.gc[data.gc$expt_rep == e &
                            data.gc$orf_name == o & data.gc$condition == c &
                            data.gc$bio_rep == b & data.gc$pos == p,]
          if (sum(is.na(temp$average)) <= 5) {
            lo <- loess.smooth(temp$hours, log(temp$average),
                               span = 0.6, evaluation = 100, degree = 2,
                               family = 'gaussian')
            data.pred <- cbind(data.pred,exp(lo$y))
            col.names <- cbind(col.names,paste(e,o,c,b,p,sep = ','))
            
            # temp.plot <- ggplot() +
            #   geom_line(data = data.frame(lo), aes(x = x, y = y)) +
            #   geom_point(data = temp, aes(x = hours, y = log(average)), shape = 1) +
            #   labs(y = 'Log( Colony Size (pixels) )',
            #        x = 'Time (hours)',
            #        title = paste(o,c,b,sep = ' | ')) +
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
            #                                   margin = margin(0.1,0,0.1,0, "mm"))) +
            #   coord_cartesian(ylim = c(3,9))
            # ggsave(sprintf("%s/%s/growthcurves/%s_%s_%s.jpg",fig_path, expt.name, o, c, b), temp.plot,
            #        height = 100, width = 100, units = 'mm',
            #        dpi = 600)
          }
        }
      }
    }
  }
}
data.pred <- cbind(lo$x, data.pred)
data.pred <- data.frame(data.pred)
colnames(data.pred) <- c('Time',col.names)
head(data.pred)

##### SAVE GROWTH CURVE DATA
save(data.pred, file = sprintf("%s/%s/growthcurve.RData", out_path, expt.name))

##### GROWTH CURVE ANALYSIS
load(file = sprintf("%s/%s/growthcurve.RData", out_path, expt.name))
data.pred <- data.pred[,colSums(log(data.pred[,-1]) < 5) <= 70]
# data.pred[,-1] <- log(data.pred[,-1])
gc.gr.res <- NULL
for (i in 2:dim(data.pred)[2]) {
  fit0 <- fit_easylinear(data.pred$Time[2:dim(data.pred)[1]], data.pred[2:dim(data.pred)[1],i], h = 20, quota = 1);
  
  temp_res <- data.frame(colnames(data.pred[i]), maxgr = coef(fit0)[[3]],
                         dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]])
  gc.gr.res <- rbind(gc.gr.res, temp_res)
}
gc.gr.res <- data.frame(gc.gr.res)
colnames(gc.gr.res) <- c('sample','gr','dtime','lag')
temp <- str_split(gc.gr.res$sample, ',', simplify = T)
colnames(temp) <- c('expt_rep','orf_name','condition','bio_rep','pos')
gc.gr.res <- cbind(temp, gc.gr.res)
head(gc.gr.res)

gc.res <- SummarizeGrowthByPlate(data.pred)
temp <- str_split(gc.res$sample, ',', simplify = T)
colnames(temp) <- c('expt_rep','orf_name','condition','bio_rep','pos')
gc.res <- cbind(temp, gc.res)
head(gc.res)

gc.res <- merge(gc.res, gc.gr.res, by = c('expt_rep','sample','orf_name','condition','bio_rep','pos'))
head(gc.res)

temp.count <- gc.res %>%
  group_by(expt_rep ,orf_name, condition) %>%
  count() %>% 
  data.frame()
gc.res <- merge(gc.res, temp.count, by = c('expt_rep', 'orf_name', 'condition'))
head(gc.res)

gc.res$pos <- as.numeric(as.character(gc.res$pos))
gc.res$orf_name <- as.character(gc.res$orf_name)
gc.res$condition <- as.character(gc.res$condition)
gc.res$expt_rep <- factor(gc.res$expt_rep, levels = c("1","2","3"))
gc.res$orf_name <- factor(gc.res$orf_name, levels = strain.levels)
gc.res$condition <- factor(gc.res$condition, levels = condition.levels)

write.csv(gc.res, file = sprintf('%s/%s/GROWTH_CURVE_RESULTS.csv', res_path, expt.name))

# plot.gc.auc <- gc.res[gc.res$n > 200,] %>%
#   ggplot(aes(x = orf_name, y = auc_e)) +
#   geom_jitter(aes(col = bio_rep)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_violin(fill = 'transparent') +
#   # stat_compare_means(method = 't.test', ref.group = 'FY4', label = "p.signif") +
#   scale_color_discrete(name = 'Biological Replicate') +
#   facet_wrap(.~condition, nrow = 5) +
#   labs(x = 'Strain',
#        y = 'AUC',
#        title = 'Area Under the Curve (AUC)') +
#   # coord_flip() +
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
# 
# plot.gc.tgen <- gc.res[gc.res$n > 200 & gc.res$t_gen <= 500,] %>%
#   ggplot(aes(x = orf_name, y = gr)) +
#   geom_jitter(aes(col = bio_rep)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_violin(fill = 'transparent') +
#   scale_color_discrete(name = 'Biological Replicate') +
#   facet_wrap(.~condition, nrow = 5, scales = 'free_y') +
#   labs(x = 'Strain',
#        y = 'Rate of Change in Colony Size (pix/hours)',
#        title = 'Growth Rate') +
#   # coord_flip() +
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
# 
# plot.gc.res <- ggpubr::ggarrange(plot.gc.tgen, plot.gc.auc,
#           ncol = 2,
#           common.legend = T, legend = 'bottom')
# ggsave(sprintf("%s/%s/GROWTH_CURVES_RESULTS2.jpg", fig_path, expt.name), plot.gc.res,
#        height = 200, width = 150, units = 'mm',
#        dpi = 600)
# 
# ## CORR BETWEEN AUC AND GROWTH RATE
# plot.auc.gr.corr <- ggplot(gc.res[gc.res$n > 200 & abs(gc.res$dtime) <= 500,], 
#        aes(x = auc_e, y = log(2)/dtime)) + 
#   geom_point(aes(col = orf_name, shape = condition)) +
#   stat_cor(method = 'spearman') +
#   geom_smooth(method = 'lm') +
#   labs(x = 'Area Under the Curve',
#        y = 'Growth Rate') +
#   scale_shape_discrete(name = 'Condition') +
#   scale_color_discrete(name = 'Strain') +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'right',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%s/%s/GROWTH_CURVES_AUC_GR_CORR.jpg", fig_path, expt.name),
#        plot.auc.gr.corr,
#        height = 120, width = 150, units = 'mm',
#        dpi = 600)

## AUC AND DT T.TESTS
strain.pairs <- permutations(n = length(unique(data$orf_name)), r = 2,
                             v = as.character(unique(data$orf_name)),
                             repeats.allowed = F) %>%
  data.frame(stringsAsFactors = F)

gc.anova.res <- NULL
for (e in unique(gc.res$expt_rep)) {
  for (c in unique(gc.res$condition[gc.res$expt_rep == e])) {
    for (i in seq(1,dim(strain.pairs)[1])) {
      o1 <- strain.pairs[i,1]
      o2 <- strain.pairs[i,2]
      
      if (length(gc.res$t_gen[gc.res$expt_rep == e & gc.res$condition == c & gc.res$orf_name == o1]) > 1 &
          length(gc.res$t_gen[gc.res$expt_rep == e & gc.res$condition == c & gc.res$orf_name == o2]) > 1) {
        res.aov.gr <- gc.res[gc.res$expt_rep == e & gc.res$condition == c & gc.res$orf_name %in% c(o1, o2),] %>%
          anova_test(gr ~ orf_name * bio_rep)
        res.aov.t_gen <- gc.res[gc.res$expt_rep == e & gc.res$condition == c & gc.res$orf_name %in% c(o1, o2),] %>%
          anova_test(t_gen ~ orf_name * bio_rep)
        res.aov.auc_e <- gc.res[gc.res$expt_rep == e & gc.res$condition == c & gc.res$orf_name %in% c(o1, o2),] %>%
          anova_test(auc_e ~ orf_name * bio_rep)
        
        gc.anova.res <- rbind(gc.anova.res, cbind(e, c, o1, o2,
                                                  res.aov.gr$p[res.aov.gr$Effect == 'orf_name'],
                                                  res.aov.gr$p[res.aov.gr$Effect == 'orf_name:bio_rep'],
                                                  res.aov.t_gen$p[res.aov.t_gen$Effect == 'orf_name'],
                                                  res.aov.t_gen$p[res.aov.t_gen$Effect == 'orf_name:bio_rep'],
                                                  res.aov.auc_e$p[res.aov.auc_e$Effect == 'orf_name'],
                                                  res.aov.auc_e$p[res.aov.auc_e$Effect == 'orf_name:bio_rep']))
      }
    }
  }
}
colnames(gc.anova.res) <- c('expt_rep','condition','strain1','strain2',
                         'gr_between','gr_within','t_gen_between','t_gen_within','auc_e_between','auc_e_within')
gc.anova.res <- data.frame(gc.anova.res, stringsAsFactors = F)
gc.anova.res$gr_between <- as.numeric(gc.anova.res$gr_between)
gc.anova.res$gr_within <- as.numeric(gc.anova.res$gr_within)
gc.anova.res$t_gen_between <- as.numeric(gc.anova.res$t_gen_between)
gc.anova.res$t_gen_within <- as.numeric(gc.anova.res$t_gen_within)
gc.anova.res$auc_e_between <- as.numeric(gc.anova.res$auc_e_between)
gc.anova.res$auc_e_within <- as.numeric(gc.anova.res$auc_e_within)
gc.anova.res$gr_between <- p.adjust(gc.anova.res$gr_between, method = 'BH')
gc.anova.res$gr_within <- p.adjust(gc.anova.res$gr_within, method = 'BH')
gc.anova.res$t_gen_between <- p.adjust(gc.anova.res$t_gen_between, method = 'BH')
gc.anova.res$t_gen_within <- p.adjust(gc.anova.res$t_gen_within, method = 'BH')
gc.anova.res$auc_e_between <- p.adjust(gc.anova.res$auc_e_between, method = 'BH')
gc.anova.res$auc_e_within <- p.adjust(gc.anova.res$auc_e_within, method = 'BH')

gc.anova.res$strain1 <- factor(gc.anova.res$strain1, levels = strain.levels)
gc.anova.res$strain2 <- factor(gc.anova.res$strain2, levels = strain.levels)
gc.anova.res$condition <- factor(gc.anova.res$condition, levels = condition.levels)
gc.anova.res$expt_rep <- factor(gc.anova.res$expt_rep, levels = c("1","2","3"))

write.csv(gc.anova.res, file = sprintf('%s/%s/GROWTH_CURVE_ANOVA_RESULTS.csv', res_path, expt.name))

for (e in unique(gc.anova.res$expt_rep)) {
  plot.anova.res.gr <- ggplot(gc.anova.res) +
    geom_tile(aes(x = strain1, y = strain2, fill = gr_between <= 0.05),
              col = 'black') +
    # labs(title = sprintf('Experimental Replicate #%s',e)) +
    facet_wrap(.~condition, nrow = 5) +
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
  
  plot.anova.res.auc <- ggplot(gc.anova.res) +
    geom_tile(aes(x = strain1, y = strain2, fill = auc_e_between <= 0.05),
              col = 'black') +
    # labs(title = 'Area Under the Curve (AUC)') +
    facet_wrap(.~condition, nrow = 5) +
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
  
  plot.gc.anova.res <- ggpubr::ggarrange(plot.anova.res.gr, plot.anova.res.auc,
                                         ncol = 2,
                                         common.legend = T, legend = 'bottom')
  ggsave(sprintf("%s/%s/GROWTH_CURVES_ANOVA_%s.jpg",fig_path, expt.name,e), plot.gc.anova.res,
         height = 200, width = 150, units = 'mm',
         dpi = 600)
}

## VIOLIN + ANOVA
plot.gc.gr.anova <- gc.res %>%
  filter(expt_rep == e) %>%
  ggplot(aes(x = orf_name, y = gr)) +
  geom_violin(aes(fill = bio_rep),lwd = 0.1) +
  geom_point(data = gc.anova.res[gc.anova.res$strain2 == 'FY4' &
                                gc.anova.res$gr_between <= 0.05,],
             aes(x = strain1, y = 0.08), shape = 8, size = 1, col = 'red') +
  scale_fill_discrete(name = 'Biological Replicate') +
  facet_wrap(.~condition, nrow = 1) +
  labs(x = 'Strain',
       y = 'Growth Rate') +
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
ggsave(sprintf("%s/%s/GROWTH_CURVES_GROWTHRATE_ANOVA.jpg",fig_path, expt.name), plot.gc.gr.anova,
       height = 80, width = 200, units = 'mm',
       dpi = 600)

plot.gc.auc.anova <- gc.res %>%
  ggplot(aes(x = orf_name, y = auc_e)) +
  geom_violin(aes(fill = bio_rep),lwd = 0.1) +
  geom_point(data = gc.anova.res[gc.anova.res$strain2 == 'FY4' &
                                   gc.anova.res$auc_e_between <= 0.05,],
             aes(x = strain1, y = 1200000), shape = 8, size = 1, col = 'red') +
  scale_fill_discrete(name = 'Biological Replicate') +
  facet_wrap(.~condition, nrow = 1) +
  labs(x = 'Strain',
       y = 'Area Under the Curve') +
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
ggsave(sprintf("%s/%s/GROWTH_CURVES_AUC_ANOVA.jpg",fig_path, expt.name), plot.gc.auc.anova,
       height = 80, width = 200, units = 'mm',
       dpi = 600)


##### GROWTH CURVE ANALYSIS SUMMARY
head(gc.res)
gc.res.sum <- gc.res %>%
  filter(n > 10) %>%
  group_by(expt_rep, condition, orf_name, bio_rep) %>%
  summarise(t_gen = median(t_gen, na.rm = T), dtime = median(dtime, na.rm = T),
            gr = median(gr, na.rm = T), auc_e = median(auc_e, na.rm = T)) %>%
  data.frame()

plot.gc.auc <- gc.res.sum %>%
  ggplot(aes(x = orf_name, y = auc_e)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col = bio_rep)) +
  # geom_violin(fill = 'transparent') +
  stat_compare_means(method = 't.test', ref.group = 'FY4', label = "p.signif",
                     size = 1.5, label.y = 1000000) +
  scale_color_discrete(name = 'Biological Replicate') +
  facet_wrap(.~expt_rep*condition, nrow = 3) +
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
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1000000))

plot.gc.tgen <- gc.res.sum %>%
  ggplot(aes(x = orf_name, y = gr)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col = bio_rep)) +
  # geom_violin(fill = 'transparent') +
  scale_color_discrete(name = 'Biological Replicate') +
  stat_compare_means(method = 't.test', ref.group = 'FY4', label = "p.signif",
                     size = 1.5, label.y = 0.06) +
  facet_wrap(.~expt_rep*condition, nrow = 3) +
  labs(x = 'Strain',
       y = 'Rate of Change in Colony Size (pix/hours)',
       title = 'Growth Rate') +
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
                                 ncol = 2, align = 'hv',
                                 common.legend = T, legend = 'bottom')
ggsave(sprintf("%s/%s/GROWTH_CURVES_RESULTS.jpg", fig_path, expt.name), plot.gc.res,
       height = 160, width = 300, units = 'mm',
       dpi = 600)


##### GROWTH CURVE AND COLONY SIZE SUMMARY
all.res.sum <- merge(gc.res.sum, data.sum %>% filter(hours == 327.88), by = c('condition','orf_name','bio_rep'), all = T)
write.csv(all.res.sum, file = sprintf('%s/%s/ALL_RESULTS_SUMMARY.csv', res_path, expt.name))

