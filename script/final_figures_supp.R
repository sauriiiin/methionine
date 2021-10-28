##### FINAL SUPPLEMENTARY FIGURES FOR METHIONINE PAPER
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 09/27/2021 

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
library(rstatix)
library(gtools)
library(effsize)
library(png)
library(ggstatsplot)
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
load(file = 'figures/final/data.RData')

fig_path <- "~/R/Projects/methionine/figures/final/final/"

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### LEGEND GRABBER
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

##### COLOR LEGEND
fill.col <- data.frame(Auxotrophy = c('None','Methionine','Uracil-Leucine','Methionine-Uracil-Leucine','Unknown'),
                       Color = c('#FFC107','#536DFE','#E040FB','#FF5722','#BDBDBD'))
legend.aux <- fill.col %>%
  ggplot(aes(x = seq(1,5), y = 1, fill = Auxotrophy)) +
  geom_tile() +
  scale_fill_manual(limits = c('None',
                               'Methionine',
                               'Uracil-Leucine',
                               'Methionine-Uracil-Leucine',
                               'Unknown'),
                    values = c('None' = '#FFC107',
                               'Methionine' = '#536DFE',
                               'Uracil-Leucine' = '#E040FB',
                               'Methionine-Uracil-Leucine' = '#FF5722',
                               'Unknown' = '#BDBDBD'),
                    position = 'bottom') +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(fill = guide_legend(nrow=2, byrow=F, order = 1))

legend.aux <- g_legend(legend.aux)

##### CARBON SOURCE - SD
data.cbn.leu$orf_name[data.cbn.leu$orf_name == 'met15'] <- 'FY4-met15del'
figScbn1 <- rbind(data.cbn[,c(1,2,5,6,10)],
                  data.cbn.leu[,c(8,11,10,13,14)] %>%
                    filter(condition == 'SDmLeu')) %>%
  filter(base == 'SD', orf_name %in% c('FY4','FY4-met15del')) %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  # geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(fill = '#9E9E9E', size = 0.2, outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(labels = c('PlM_Glu' = '+ Met\n+ Ura\n+ Leu',
                              'MiM_Glu' = '- Met\n+ Ura\n+ Leu',
                              'MiU_Glu' = '+ Met\n- Ura\n+ Leu',
                              'MiM_Gal' = '- Met\n+ Ura\n+ Leu',
                              'MiU_Gal' = '+ Met\n- Ura\n+ Leu',
                              'MiM_Et' = '- Met\n+ Ura\n+ Leu',
                              'MiU_Et' = '+ Met\n- Ura\n+ Leu',
                              'SCmLeu' = '+ Met\n+ Ura\n- Leu',
                              'SDmLeu' = '+ Met\n+ Ura\n- Leu',
                              'SDpMET' = '+ Met\n+ Ura\n+ Leu',
                              'YPDA' = '+ Met\n+ Ura\n+ Leu')) +
  labs(x = 'Strain', y= 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_grid(orf_name ~ carbon, scales = 'free_x', space = 'free_x',
             labeller = labeller(orf_name = c('FY4' = 'FY4',
                                              'FY4-met15del' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = element_text(size = txt,  margin = margin(0.1,0,0.1,0, "mm")),
        strip.text.y = ggtext::element_markdown(size = txt, colour = 'black',
                                              margin = margin(0.1,0,0.1,0, "mm")))

figScbn1.g <- ggplot_gtable(ggplot_build(figScbn1))
stripr <- which(grepl('strip-r', figScbn1.g$layout$name))

figScbn1.g$grobs[[21]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
figScbn1.g$grobs[[21]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'

figScbn1.g$grobs[[22]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
figScbn1.g$grobs[[22]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'

# for (s in stripr) {
#   l <- figScbn1.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$label
#   a <- strain.labs.cbn$auxotrophy[strain.labs.cbn$orf_name == l]
#   
#   if (a == 'Methionine') {
#     figScbn1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
#   } else if (a == 'None') {
#     figScbn1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
#   } else if (a == 'Uracil') {
#     figScbn1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#E040FB'
#   } else {
#     figScbn1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF5722'
#   }
# }
grid.draw(figScbn1.g)

figScbn_cols <- readPNG('figures/final/final/Slide1.png')
figScbn_cols <- ggplot() + 
  background_image(figScbn_cols) +
  theme(plot.margin = margin(t=0, l=15, r=15, b=0, unit = "mm"),
        plot.background = element_blank())


figScbnSD <- annotate_figure(cowplot::plot_grid(figScbn_cols,
                                                figScbn1.g,
                                                legend.aux,
                                                labels = c('A','B',''),
                                                label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
                                                ncol = 1, rel_heights = c(7,7,1)),
                             top = text_grob('Synthetic Defined Media',
                                             face = "bold", family = "sans", size = titles))
ggsave(sprintf("%s/FigureS1.jpg",fig_path), figScbnSD,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)


##### CARBON SOURCE - SC
figScbn3 <- rbind(data.cbn[,c(1,2,5,6,10)],
                  data.cbn.leu[,c(8,11,10,13,14)] %>%
                    filter(condition == 'SCmLeu')) %>%
  filter(base == 'SC', orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  geom_boxplot(fill = '#9E9E9E', size = 0.2, outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(limits = c('BY4742','BY4741','FY4','FY4-met15del'),
                   labels = c('FY4' = 'FY4',
                              'FY4-met15del' = 'FY4-*met15Δ*',
                              'BY4742' = 'BY4742',
                              'BY4741' = 'BY4741')) +
  labs(x = 'Strain', y= 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_wrap(.~condition*carbon, scales = 'free_x',
             nrow = 1,
             labeller = labeller(condition = c('GLU' = 'SC - Met + Ura + Leu',
                                               'SCmLeu' = 'SC + Met + Ura - Leu',
                                               'GAL' = 'SC - Met + Ura + Leu'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(angle = 30, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = element_text(size = txt,  margin = margin(0.1,0,0.1,0, "mm")))



# figScbnSC <- annotate_figure(cowplot::plot_grid(figScbn3 + theme(legend.position = 'none'), legend.aux,
#                                               ncol = 1, rel_heights = c(7,1)),
#                            top = text_grob('Synthetic Complete Media',
#                                            face = "bold", family = "sans", size = titles))

ggsave(sprintf("%s/FigureS3.jpg",fig_path), 
       figScbn3,
       height = one.c, width = two.c, units = 'mm',
       dpi = 600)





##### MET3 MUTANT
figScbnMet3 <- rbind(data.cbn[,c(1,2,5,6,10)],
                     data.cbn.leu[,c(8,11,10,13,14)] %>%
                       filter(condition == 'SDmLeu')) %>%
  filter(base == 'SD', orf_name == 'FY4-met3del') %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  # geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(fill = '#9E9E9E', size = 0.2, outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(labels = c('PlM_Glu' = '+ Met\n+ Ura\n+ Leu',
                              'MiM_Glu' = '- Met\n+ Ura\n+ Leu',
                              'MiU_Glu' = '+ Met\n- Ura\n+ Leu',
                              'MiM_Gal' = '- Met\n+ Ura\n+ Leu',
                              'MiU_Gal' = '+ Met\n- Ura\n+ Leu',
                              'MiM_Et' = '- Met\n+ Ura\n+ Leu',
                              'MiU_Et' = '+ Met\n- Ura\n+ Leu',
                              'SCmLeu' = '+ Met\n+ Ura\n- Leu',
                              'SDmLeu' = '+ Met\n+ Ura\n- Leu',
                              'SDpMET' = '+ Met\n+ Ura\n+ Leu',
                              'YPDA' = '+ Met\n+ Ura\n+ Leu')) +
  labs(title = 'FY4-*met3Δ* in Synthetic Defined Media',
       x = 'Strain', y= 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_grid(.~carbon, scales = 'free_x', space = 'free_x',
             labeller = labeller(orf_name = c('FY4' = 'FY4',
                                              'FY4-met15del' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.title = ggtext::element_markdown(size = titles, hjust = 0.5, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,  margin = margin(0.1,0,0.1,0, "mm")))


ggsave(sprintf("%s/FigureS7.jpg",fig_path), 
       figScbnMet3,
       height = one.c, width = two.c, units = 'mm',
       dpi = 600)


##### RESPIRATION - BYs
data.res.gc.sum$methionine <- factor(data.res.gc.sum$methionine, levels = c('+Met','-Met'))

figSres1 <- merge(data.res.gc[str_detect(data.res.gc$orf_name, 'BY'),],
                  strain.labs.res, by = 'orf_name') %>%
  filter(expt_rep %in% c("1","2"), condition != 'SD+Met-Ura+Glu') %>% 
  ggplot(aes(x = Time, y = rel_cs)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               aes(group = orf_name, fill = met_aux), geom="ribbon", alpha = 0.4) +
  stat_summary(aes(group = orf_name, col = met_aux, linetype = pet), fun=mean, geom="line", lwd =0.7) +
  geom_text_repel(data = merge(data.res.gc[str_detect(data.res.gc$orf_name, 'BY') &
                                             data.res.gc$Time == max(data.res.gc$Time),],
                               strain.labs.res, by = 'orf_name') %>%
                    filter(expt_rep %in% c("1","2"), condition != 'SD+Met-Ura+Glu') %>%
                    group_by(orf_name, condition, Time, base, methionine, carbon, cysteine, labels, met_aux, pet) %>%
                    summarize(rel_cs = mean(rel_cs, na.rm = T), .groups = 'keep'),
                  aes(x = Time, y = rel_cs, label = labels),
                  parse = T, size = 2.5) +
  facet_grid(~carbon*methionine,
             labeller = labeller(methionine = c('+Met' = 'SD + Met',
                                                '-Met' = 'SD - Met'))) +
  scale_color_manual(name = 'Presumed Auxotrophy',
                     values = c('Prototroph' = '#FFC107',
                                'Presumed Auxotroph' = '#536DFE',
                                'Uracil-Leucine' = '#E040FB',
                                'Methionine-Uracil-Leucine' = '#FF5722'),
                     limits = c('Prototroph','Presumed Auxotroph'),
                     labels = c('Prototroph'='None',
                                'Presumed Auxotroph'='Methionine',
                                'Uracil-Leucine',
                                'Methionine-Uracil-Leucine')) +
  scale_fill_manual(name = 'Auxotrophy',
                    values = c('Prototroph' = '#FFC107',
                               'Presumed Auxotroph' = '#536DFE'),
                    limits = c('Prototroph','Presumed Auxotroph'),
                    guide = F) +
  scale_linetype_manual(name = 'Petite',
                        values = c('Yes' = 'dotdash',
                                   'No' = 'solid'),
                        limits = c('Yes','No')) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  labs(x = 'Time (hours)', y= 'Relative Colony Size') +
  coord_cartesian(xlim = c(0,350),
                  ylim = c(0,1.3)) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=1, byrow=TRUE, order = 1),
         linetype = guide_legend(nrow=1, byrow=TRUE, order = 2))

ggsave(sprintf("%s/Respiration_BYs.jpg",fig_path), figSres1,
       height = one.c, width = two.c, units = 'mm',
       dpi = 600)


# figSres2 <- readPNG('figures/final/final/RESP_STATS_SUPP.png')
# figSres2 <- ggplot() + 
#   background_image(figSres2) +
#   theme(plot.margin = margin(t=0, l=10, r=10, b=0, unit = "mm"),
#         plot.background = element_blank())
# 
# legend.pet <- g_legend(figSres1)
# figSres <- cowplot::plot_grid(figSres1 + theme(legend.position = 'none'),
#                               cowplot::plot_grid(legend.aux, legend.pet, nrow = 1, rel_widths = c(3,1)),
#                               figSres2,
#                               ncol = 1, rel_heights = c(6,1,6))
# ggsave(sprintf("%s/FigureS4.jpg",fig_path), figSres,
#        height = two.c, width = two.c, units = 'mm',
#        dpi = 600)


##### REPEATED PINNING
strain.labs.rp <- strain.labs.cbn
strain.labs.rp$orf_name[strain.labs.rp$orf_name == 'FY4-met15del'] <- 'FY4-met15D'
data.rp$aux <- factor(data.rp$aux, levels = c('Ura', 'Met'))
data.rp$carbon <- factor(data.rp$carbon, levels = c('Glu', 'Gal'))

data.rp$cum_hrs <- data.rp$hours
data.rp$cum_hrs[data.rp$pin == 2] <- data.rp$cum_hrs[data.rp$pin == 2] +
  max(data.rp$cum_hrs[data.rp$pin == 1], na.rm = T)
data.rp$cum_hrs[data.rp$pin == 3] <- data.rp$cum_hrs[data.rp$pin == 3] +
  max(data.rp$cum_hrs[data.rp$pin == 2], na.rm = T)
data.rp$cum_hrs[data.rp$pin == 4] <- data.rp$cum_hrs[data.rp$pin == 4] +
  max(data.rp$cum_hrs[data.rp$pin == 3], na.rm = T)
data.rp$cum_hrs[data.rp$pin == 5] <- data.rp$cum_hrs[data.rp$pin == 5] +
  max(data.rp$cum_hrs[data.rp$pin == 4], na.rm = T)
data.rp$cum_hrs[data.rp$pin == 6] <- data.rp$cum_hrs[data.rp$pin == 6] +
  max(data.rp$cum_hrs[data.rp$pin == 5], na.rm = T)
data.rp$cum_hrs[data.rp$pin == 7] <- data.rp$cum_hrs[data.rp$pin == 7] +
  max(data.rp$cum_hrs[data.rp$pin == 6], na.rm = T)

figSrp1 <- merge(data.rp, strain.labs.rp,
      by = 'orf_name') %>%
  filter(aux == 'Met', orf_name %in% c('BY4742')) %>%
  ggplot(aes(x = cum_hrs, y = average)) +
  stat_summary(aes(fill = auxotrophy, group = pin),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(aes(col = auxotrophy, group = pin),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Met', orf_name %in% c("BY4741")),
               aes(fill = auxotrophy, group = pin),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Met', orf_name %in% c("BY4741")),
               aes(col = auxotrophy, group = pin),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Met', orf_name %in% c("FY4")),
               aes(fill = auxotrophy, group = pin),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Met', orf_name %in% c("FY4")),
               aes(col = auxotrophy, group = pin),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Met', orf_name %in% c("FY4-met15D")),
               aes(fill = auxotrophy, group = pin),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Met', orf_name %in% c("FY4-met15D")),
               aes(col = auxotrophy, group = pin),
               fun=mean, geom="line", lwd =1) +
  geom_text_repel(data = merge(data.rp, strain.labs.rp,
                               by = 'orf_name') %>%
                    filter(aux == 'Met', orf_name %in% c('BY4741',"BY4742","FY4",'FY4-met15D'), 
                           hours == max_hrs, pin == 7) %>%
                    group_by(orf_name, labels, carbon, cum_hrs) %>%
                    summarize(average = mean(average, na.rm = T), .groups = 'keep'),
                  aes(x = cum_hrs, y = average + 200, label = labels),
                  parse = T, size = 2.5) +
  scale_x_continuous(breaks = seq(-1000,1000,100)) +
  scale_color_manual(name = '',
                     values = c('None' = '#FFC107',
                                'Methionine' = '#536DFE',
                                'Uracil' = '#E040FB',
                                'Both' = '#FF5722')) +
  scale_fill_manual(name = 'Auxotrophy',
                    values = c('None' = '#FFC107',
                               'Methionine' = '#536DFE',
                               'Uracil' = '#E040FB',
                               'Both' = '#FF5722'),
                    guide = F) +
  labs(x = 'Time (hours)', y = 'Colony Size (pixels)',
       title = 'SD-Met') +
  facet_grid(carbon~., labeller = labeller(carbon = c('Glu'='Glucose',
                                                      'Gal'='Galactose'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'none',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                                margin = margin(0.1,0,0.1,0, "mm")),
        strip.text.y = element_text(size = txt, margin = margin(0.1,0,0.1,0, "mm")))

figSrp2 <- readPNG('figures/final/final/FigureS_RP_MM.png')
figSrp2 <- ggplot() + 
  background_image(figSrp2) +
  theme(plot.margin = margin(t=0, l=10, r=10, b=0, unit = "mm"),
        plot.background = element_blank())

figSrp3 <- merge(data.rp, strain.labs.rp,
                 by = 'orf_name') %>%
  filter(aux == 'Ura', orf_name %in% c('BY4742')) %>%
  ggplot(aes(x = cum_hrs, y = average)) +
  stat_summary(aes(fill = auxotrophy, group = pin),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(aes(col = auxotrophy, group = pin),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Ura', orf_name %in% c("BY4741")),
               aes(fill = auxotrophy, group = pin),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Ura', orf_name %in% c("BY4741")),
               aes(col = auxotrophy, group = pin),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Ura', orf_name %in% c("FY4")),
               aes(fill = auxotrophy, group = pin),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Ura', orf_name %in% c("FY4")),
               aes(col = auxotrophy, group = pin),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Ura', orf_name %in% c("FY4-met15D")),
               aes(fill = auxotrophy, group = pin),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Ura', orf_name %in% c("FY4-met15D")),
               aes(col = auxotrophy, group = pin),
               fun=mean, geom="line", lwd =1) +
  geom_text_repel(data = merge(data.rp, strain.labs.rp,
                               by = 'orf_name') %>%
                    filter(aux == 'Ura', orf_name %in% c('BY4741',"BY4742","FY4",'FY4-met15D'), 
                           hours == max_hrs, pin == 7) %>%
                    group_by(orf_name, labels, carbon, cum_hrs) %>%
                    summarize(average = mean(average, na.rm = T), .groups = 'keep'),
                  aes(x = cum_hrs, y = average + 200, label = labels),
                  parse = T, size = 2.5) +
  scale_x_continuous(breaks = seq(-1000,1000,100)) +
  scale_color_manual(name = '',
                     values = c('None' = '#FFC107',
                                'Methionine' = '#536DFE',
                                'Uracil' = '#E040FB',
                                'Both' = '#FF5722')) +
  scale_fill_manual(name = 'Auxotrophy',
                    values = c('None' = '#FFC107',
                               'Methionine' = '#536DFE',
                               'Uracil' = '#E040FB',
                               'Both' = '#FF5722'),
                    guide = F) +
  labs(x = 'Time (hours)', y = 'Colony Size (pixels)',
       title = 'SD-Ura') +
  facet_grid(carbon~., labeller = labeller(carbon = c('Glu'='Glucose',
                                                      'Gal'='Galactose'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'none',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                                margin = margin(0.1,0,0.1,0, "mm")),
        strip.text.y = element_text(size = txt, margin = margin(0.1,0,0.1,0, "mm")))


figSrp4 <- readPNG('figures/final/final/FigureS_RP_MU.png')
figSrp4 <- ggplot() + 
  background_image(figSrp4) +
  theme(plot.margin = margin(t=0, l=10, r=10, b=0, unit = "mm"),
        plot.background = element_blank())

figSrpMM <- cowplot::plot_grid(legend.aux, figSrp1, figSrp2,
                             ncol = 1, rel_heights = c(0.3,1,0.9),
                             labels = c('','A','B'),
                             label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/FigureS5_1.jpg",fig_path), figSrpMM,
       height = 170, width = two.c, units = 'mm',
       dpi = 600)


figSrpMU <- cowplot::plot_grid(legend.aux, figSrp3, figSrp4,
                               ncol = 1, rel_heights = c(0.3,1,0.9),
                               labels = c('','A','B'),
                               label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/FigureS5_2.jpg",fig_path), figSrpMU,
       height = 170, width = two.c, units = 'mm',
       dpi = 600)


##### BISMUTH
data.bis$Strain <- as.character(data.bis$Strain)
data.bis.3 <- data.frame(rbind(c('FY4-met5D',1.0,'BiGGY','met5'),
                 c('FY4-met5D',1.0,'SD-Met-Cys+Bi','met5'),
                 c('FY4-met10D',1.0,'BiGGY','met10'),
                 c('FY4-met10D',1.0,'SD-Met-Cys+Bi','met10'),
                 c('FY4-yllD',4.0,'BiGGY','yll'),
                 c('FY4-yllD',4.0,'SD-Met-Cys+Bi','yll')))
colnames(data.bis.3) <- colnames(data.bis)
data.bis <- rbind(data.bis, data.bis.3)
data.bis$Strain <- factor(data.bis$Strain,
                          levels = c('FY4','BY4742','FY4-met12D','FY4-met2D','FY4-met6D','FY4-met13D',
                                     'FY4-cys4D','FY4-yllD','FY4-str3D','FY4-met3D','FY4-met5D','FY4-met10D',
                                     'FY4-met15D','BY4741'))
data.bis$HS <- as.numeric(data.bis$HS)

figSBi <- data.bis %>%
  filter(Condition == 'SD-Met-Cys+Bi') %>%
  ggplot(aes(x = Strain, y = HS)) +
  stat_summary(data = data.bis %>%
                 filter(Condition == 'SD-Met-Cys+Bi') %>%
                 group_by(Condition, Strain) %>%
                 summarize(HS = mean(HS, na.rm = T), .groups = 'keep'),
               aes(fill = round(HS)), col = 'white', alpha = 0.9, size = 1,
               fun = mean, geom = "bar") +
  # stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  # scale_x_discrete(labels = c('BiGGY' = 'BiGGY', 'SD-Met-Cys+Bi' = 'SD-Met+Bi'),
  #                  limits = c('BiGGY', 'SD-Met-Cys+Bi')) +
  scale_y_continuous(breaks = seq(0,10,1)) +
  scale_fill_gradient(low = "#D7CCC8", high = "#5D4037", guide = F) +
  labs(y = 'Relative Hydrogen Sulfide',
       x = 'Strains',
       title = 'SD-Met+Glu+Bi') +
  scale_x_discrete(limits = c('BY4742','BY4741','FY4','FY4-met15D',
                              'FY4-met2D','FY4-met3D','FY4-met5D','FY4-met6D',
                              'FY4-met10D','FY4-met12D','FY4-met13D','FY4-cys4D',
                              'FY4-str3D','FY4-yllD'),
                   labels = c('FY4'='FY4',
                              'FY4-met12D'='FY4-*met12Δ*',
                              'FY4-str3D'='FY4-*str3Δ*',
                              'FY4-met3D'='FY4-*met3Δ*',
                              'FY4-met15D'='FY4-*met15Δ*',
                              'FY4-met2D'='FY4-*met2Δ*',
                              'FY4-met6D'='FY4-*met6Δ*', 
                              'FY4-met13D'='FY4-*met13Δ*',
                              'FY4-cys4D'='FY4-*cys4Δ*',
                              'FY4-yllD'='FY4-*yll058wΔ*',
                              'FY4-met5D' = 'FY4-*met5Δ*',
                              'FY4-met10D' = 'FY4-*met10Δ*',
                              'BY4742'='BY4742',
                              'BY4741'='BY4741')) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'black',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(1,7))

# figSBi.g <- ggplot_gtable(ggplot_build(figSBi))
# stripr <- which(grepl('strip-t', figSBi.g$layout$name))
# 
# for (s in stripr) {
#   # figSBi.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
#   l <- figSBi.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
#   
#   if (length(l) != 0) {
#     if (l == 'FY4-') {
#       l2 <- figSBi.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$label
#       a <- strain.labs.bis$auxotrophy[str_detect(strain.labs.bis$labels, l2)]
#     } else {
#       a <- strain.labs.bis$auxotrophy[strain.labs.bis$orf_name == l]
#     }
#     
#     if (a == 'Presumed Auxotroph') {
#       figSBi.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
#     } else if (a == 'Prototroph') {
#       figSBi.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
#     } else {
#       figSBi.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#BDBDBD'
#     }
#   }
# }
# 
# ggsave(sprintf("%s/FigureS7.jpg",fig_path), figSBi.g,
#        height = one.5c, width = one.5c, units = 'mm',
#        dpi = 600)
figSBicp <- readPNG('figures/final/final/Bismuth.png')
figSBicp <- ggplot() + 
  background_image(figSBicp) +
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
        plot.background = element_blank())

figSBis <- plot_grid(figSBicp, figSBi,
                     ncol = 1, rel_heights = c(2.5,1),
                     labels = c('A','B'),
                     label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/FigureSBi.jpg",fig_path), figSBis,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)

##### MET5/10 MUTANTS
data.510$orf_name <- factor(data.510$orf_name, levels = c('FY4','BY4741_met10del','BY4741_met5del',
                                                          'met15del','met10del','met5del'))

figSmet510 <- merge(data.510, strain.labs.510, by = 'orf_name') %>%
  filter(hours == 160) %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  geom_boxplot(fill = '#616161', outlier.shape = NA, size = 0.2,
               position = position_dodge(width = 0.85)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(limits = c('YPDA', 
                              'SD+Met',
                              'SD-Met',
                              'SD-Ura')) +
  labs(x = 'Strain',
       y = 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_wrap(~orf_name, ncol = 3, 
             labeller = labeller(orf_name = c('FY4'='FY4',
                                              'met15del'='FY4-*met15Δ*',
                                              'met5del'='FY4-*met5Δ*',
                                              'met10del'='FY4-*met10Δ*',
                                              'BY4741_met5del'='BY4741-*met5Δ*',
                                              'BY4741_met10del'='BY4741-*met10Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(0.1,0,0.1,0, "mm")))

ggsave(sprintf("%s/Met510.jpg",fig_path), figSmet510,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)

##### NO SULFATE
strain.labs.ns <- strain.labs.cbn
strain.labs.ns$orf_name[strain.labs.ns$orf_name == 'FY4-met3del'] <- 'FY4_met3del'
strain.labs.ns$orf_name[strain.labs.ns$orf_name == 'FY4-met15del'] <- 'FY4_met15del'


figSns1 <- data.ns[!(data.ns$id == 'Agar_Difco_+sulfate_+Met' & data.ns$stage == 'S1'),] %>%
  filter(id %in% c('Agar_Difco_+sulfate_+Met',
                   'Agarose_Home_+sulfate_-Met',
                   'Agarose_Home_-sulfate_-Met'),
         orf_name %in% c('FY4', 'BY4742')) %>%
  group_by(stage, base, ynb_type, sulfate, methionine, orf_name, hours, cum_hrs, id) %>%
  summarize(average = median(average, na.rm = T), .groups = 'keep') %>%
  data.frame() %>%
  ggplot(aes(x = cum_hrs, y = average)) +
  geom_line(aes(col = id)) +
  # stat_summary(aes(fill = id, group = stage),
  #              fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  # stat_summary(aes(col = id, group = stage),
  #              fun=mean, geom="line", lwd =1) +
  scale_x_continuous(breaks = seq(0,1000,50)) +
  scale_color_manual(name = 'Condition',
                     breaks = c('Agarose_Home_-sulfate_-Met',
                                'Agarose_Home_+sulfate_-Met',
                                'Agar_Difco_+sulfate_+Met'),
                     values = c('Agar_Difco_+sulfate_+Met' = '#212121',
                                'Agarose_Home_+sulfate_-Met' = '#795548',
                                'Agarose_Home_-sulfate_-Met' = '#FF5722'),
                     labels = c('Agar_Difco_+sulfate_+Met' = 'SD + Met w/ Inorganic Sulfates',
                                'Agarose_Home_+sulfate_-Met' = 'SD - Met w/ Inorganic Sulfates',
                                'Agarose_Home_-sulfate_-Met' = 'SD - Met w/o Inorganic Sulfates'),
                     position = 'bottom') +
  scale_fill_manual(name = 'Condition',
                     breaks = c('Agarose_Home_-sulfate_-Met',
                                'Agarose_Home_+sulfate_-Met',
                                'Agar_Difco_+sulfate_+Met'),
                     values = c('Agar_Difco_+sulfate_+Met' = '#212121',
                                'Agarose_Home_+sulfate_-Met' = '#795548',
                                'Agarose_Home_-sulfate_-Met' = '#FF5722'),
                     labels = c('Agar_Difco_+sulfate_+Met' = 'SD + Met w/ Inorganic Sulfates',
                                'Agarose_Home_+sulfate_-Met' = 'SD - Met w/ Inorganic Sulfates',
                                'Agarose_Home_-sulfate_-Met' = 'SD - Met w/o Inorganic Sulfates'),
                     guide = F) +
  # scale_linetype_manual(name = 'Media Base',
  #                       values = c('Agar' = 'twodash',
  #                                  'Agarose' = 'solid')) +
  facet_wrap(.~orf_name, nrow = 2,
             labeller = labeller(orf_name = c('FY4' = 'FY4',
                                              'FY4_met15del' = 'FY4-*met15Δ*',
                                              'FY4_met3del' = 'FY4-*met3Δ*',
                                              'BY4742' = 'BY4742',
                                              'BY4741' = 'BY4741'))) +
  labs(x = 'Time (hours)',
       y = 'Colony Size (pixels)') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        plot.margin = margin(0, 5, 0, 0, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        legend.spacing.y = unit(0.5, "mm"),
        strip.text = ggtext::element_markdown(size = txt, color = 'black',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=3, byrow=TRUE, order = 1, override.aes=list(size = 3)),
         linetype = guide_legend(nrow=2, byrow=TRUE, order = 2))

figSns1.g <- ggplot_gtable(ggplot_build(figSns1 + theme(legend.position = 'none')))
stripr <- which(grepl('strip-t', figSns1.g$layout$name))

for (s in stripr) {
  l <- figSns1.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
  
  if (length(l) != 0) {
    if (l == 'FY4-') {
      l2 <- figSns1.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$label
      a <- sstrain.labs.ns$auxotrophy[strain.labs.ns$orf_name == l]
    } else {
      a <- strain.labs.ns$auxotrophy[strain.labs.ns$orf_name == l]
    }
    if (a == 'Methionine') {
      figSns1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
    } else if (a == 'None') {
      figSns1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
    } else if (a == 'Uracil') {
      figSns1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#E040FB'
    } else if (a == 'Both') {
      figSns1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF5722'
    } else {
      figSns1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF5722'
    }
  }
}
# grid.draw(figSns1.g)


figSns2 <- data.ns[!(data.ns$id == 'Agar_Difco_+sulfate_+Met' & data.ns$stage == 'S1'),] %>%
  filter(id %in% c('Agar_Difco_+sulfate_+Met',
                   'Agarose_Home_+sulfate_-Met',
                   'Agarose_Home_-sulfate_-Met'),
         orf_name %in% c('FY4_met15del','BY4741','FY4_met3del')) %>%
  group_by(stage, base, ynb_type, sulfate, methionine, orf_name, hours, cum_hrs, id) %>%
  summarize(average = median(average, na.rm = T), .groups = 'keep') %>%
  data.frame() %>%
  ggplot(aes(x = cum_hrs, y = average)) +
  geom_line(aes(col = id)) +
  # stat_summary(aes(fill = id, group = stage),
  #              fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  # stat_summary(aes(col = id, group = stage),
  #              fun=mean, geom="line", lwd =1) +
  scale_x_continuous(breaks = seq(0,1000,50)) +
  scale_color_manual(name = 'Condition',
                     breaks = c('Agarose_Home_-sulfate_-Met',
                                'Agarose_Home_+sulfate_-Met',
                                'Agar_Difco_+sulfate_+Met'),
                     values = c('Agar_Difco_+sulfate_+Met' = '#212121',
                                'Agarose_Home_+sulfate_-Met' = '#795548',
                                'Agarose_Home_-sulfate_-Met' = '#FF5722'),
                     labels = c('Agar_Difco_+sulfate_+Met' = 'SD + Met w/ Inorganic Sulfates',
                                'Agarose_Home_+sulfate_-Met' = 'SD - Met w/ Inorganic Sulfates',
                                'Agarose_Home_-sulfate_-Met' = 'SD - Met w/o Inorganic Sulfates'),
                     position = 'bottom') +
  scale_fill_manual(name = 'Condition',
                    breaks = c('Agarose_Home_-sulfate_-Met',
                               'Agarose_Home_+sulfate_-Met',
                               'Agar_Difco_+sulfate_+Met'),
                    values = c('Agar_Difco_+sulfate_+Met' = '#212121',
                               'Agarose_Home_+sulfate_-Met' = '#795548',
                               'Agarose_Home_-sulfate_-Met' = '#FF5722'),
                    labels = c('Agar_Difco_+sulfate_+Met' = 'SD + Met w/ Inorganic Sulfates',
                               'Agarose_Home_+sulfate_-Met' = 'SD - Met w/ Inorganic Sulfates',
                               'Agarose_Home_-sulfate_-Met' = 'SD - Met w/o Inorganic Sulfates'),
                    guide = F) +
  # scale_linetype_manual(name = 'Media Base',
  #                       values = c('Agar' = 'twodash',
  #                                  'Agarose' = 'solid')) +
  facet_wrap(.~orf_name, nrow = 3,
             labeller = labeller(orf_name = c('FY4' = 'FY4',
                                              'FY4_met15del' = 'FY4-*met15Δ*',
                                              'FY4_met3del' = 'FY4-*met3Δ*',
                                              'BY4742' = 'BY4742',
                                              'BY4741' = 'BY4741'))) +
  labs(x = 'Time (hours)',
       y = 'Colony Size (pixels)') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        plot.margin = margin(0, 5, 0, 0, "mm"),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        legend.spacing.y = unit(0.5, "mm"),
        strip.text = ggtext::element_markdown(size = txt, color = 'black',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=3, byrow=TRUE, order = 1, override.aes=list(size = 3)),
         linetype = guide_legend(nrow=2, byrow=TRUE, order = 2))

figSns2.g <- ggplot_gtable(ggplot_build(figSns2 + theme(legend.position = 'none')))
stripr <- which(grepl('strip-t', figSns2.g$layout$name))

for (s in stripr) {
  l <- figSns2.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
  
  if (length(l) != 0) {
    if (l == 'FY4-') {
      l2 <- figSns2.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$label
      a <- strain.labs.ns$auxotrophy[str_detect(strain.labs.ns$labels, l2)]
    } else {
      a <- strain.labs.ns$auxotrophy[str_detect(strain.labs.ns$labels, l)]
    }
    if (a == 'Methionine') {
      figSns2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
    } else if (a == 'None') {
      figSns2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
    } else if (a == 'Uracil') {
      figSns2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#E040FB'
    } else if (a == 'Both') {
      figSns2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF5722'
    } else {
      figSns2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF5722'
    }
  }
}
# grid.draw(figSns2.g)

legend.ns <- g_legend(figSns2)

figSns <- cowplot::plot_grid(cowplot::plot_grid(figSns1.g, figSns2.g,
                                                ncol = 1, rel_heights = c(2,3),
                                                labels = c('A','B'),
                                                label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                             cowplot::plot_grid(legend.aux, legend.ns, nrow = 1),
                             ncol = 1, rel_heights = c(5,0.3))
ggsave(sprintf("%s/FigureS8.jpg",fig_path), figSns,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)


##### PLASMID VALIDATION
data.pv %>%
  filter(hours %in% c(48,115), stage %in% c('WC','PS1'),
         deletion1 %in% c('met15', 'yll'), deletion2 %in% c('','yll'),
         plasmid_backbone == '2M', plasmid_orf %in% c('empty','yll'),
         orf_name %in% c('met15del_2M_empty',
                         'met15del_2M_yll')) %>%
  ggplot(aes(x = stage, y = average)) +
  geom_boxplot(aes(alpha = stage), fill = 'grey40', outlier.shape = NA, size = 0.3,
               position = position_dodge(width = 0.6)) +
  facet_wrap(.~orf_name, labeller = labeller(orf_name = c('met15del_2M_empty' = 'FY4-*met15Δ*<br />+ Empty Plasmid',
                                                          'met15del_2M_yll' = 'FY4-*met15Δ*<br />+ YLL058W Plasmid',
                                                          'met15del_ylldel_2M_empty' = 'FY4-*met15Δ yll058wΔ*<br />+ Empty Plasmid',
                                                          'met15del_ylldel_2M_yll'  = 'FY4-*met15Δ yll058wΔ*<br />+ YLL058W Plasmid',
                                                          'ylldel_2M_empty' = 'FY4-*yll058wΔ*<br />+ Empty Plasmid'))) +
  labs(y = 'Colony Size (pixels)') +
  scale_x_discrete(limits = c('WC','PS1'),
                   labels = c('WC' = 'YPDA',
                              'PS1' = 'SD-Met+Gal')) +
  scale_alpha_manual(name = 'Condition',
                     values = c('WC' = 0.4,
                                'PS1' = 1),
                     labels = c('WC' = 'YPDA',
                                'PS1' = 'SD-Met+Gal'),
                     guide = F) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")))


##### YLL FROM JAKES MUTANTS



temp1 <- merge(data.510, strain.labs.510, by = 'orf_name') %>%
  filter(hours == 160, condition %in% c('YPDA','SD-Met'), orf_name %in% c('met5del','met10del'))
temp2 <- merge(data.jm, strain.labs.jm, by = 'orf_name') %>%
  filter(orf_name %in% c('FY4','yll'), time == 't_final', attempt != 'pilot', condition %in% c('YPDA','SD-Met-Cys+Glu'))
temp2 <- temp2[,c(1:9,13,10,11,15,16)]
colnames(temp2) <- c(colnames(temp2)[1:9],'relative_fitness',colnames(temp2)[11:14])

data.jm.2 <- rbind(temp1, temp2)
data.jm.2$condition[data.jm.2$condition == 'SD-Met-Cys+Glu'] <- 'SD-Met'
data.jm.2$orf_name <- factor(data.jm.2$orf_name, levels = c('FY4','met5del','met10del','yll'))

data.bis.2 <- data.frame(orf_name = c('met10del','met5del','FY4','yll'),
                         H2S = c(1,1,4,4))


figSjmYLL <- data.jm.2 %>%
  filter(orf_name %in% c('yll','met5del','met10del')) %>%
  ggplot() +
  # geom_rect(data = data.bis.2, aes(fill = H2S),xmin = -Inf,xmax = Inf,
  #           ymin = -Inf,ymax = Inf) +
  # geom_hline(yintercept = seq(-1,2,0.2), col = '#FFECB3', size = 0.2, linetype = 'dashed') +
  # geom_vline(xintercept = seq(-1,5,1), col = '#FFECB3', size = 0.2, linetype = 'dashed') +
  # geom_boxplot(aes(x = condition, y = relative_fitness),
  #              fill = '#9E9E9E', col = 'white', outlier.shape = NA, size = 0.4,
  #              position = position_dodge(width = 0.85)) +
  geom_boxplot(aes(x = condition, y = relative_fitness),
               fill = '#9E9E9E', outlier.shape = NA, size = 0.2,
               position = position_dodge(width = 0.85)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(limits = c('YPDA', 'SD-Met'),
                   labels = c('YPDA' = 'YPDA',
                              'SD-Met' = 'SD-Met+Glu')) +
  scale_fill_gradient(name = 'Relative Hydrogen Sulfide',
                      low = "#D7CCC8", high = "#5D4037",
                      limits = c(1,7)) +
  labs(x = 'Strain',
       y = 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_wrap(~orf_name, nrow = 1, 
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
                                              'BY4741'='BY4741',
                                              'met5del'='FY4-*met5Δ*',
                                              'met10del'='FY4-*met10Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt), #angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        # legend.direction = 'vertical',
        # legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(0.1,0,0.1,0, "mm")))

figSjmYLLcp <- readPNG('figures/final/final/BiGGY_Supp.png')
figSjmYLLcp <- ggplot() + 
  background_image(figSjmYLLcp) +
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
        plot.background = element_blank())

figSjmYLLcol <- data.bis %>%
  filter(Condition == 'BiGGY', orf_name %in% c('yll','met5','met10')) %>%
  ggplot(aes(x = Strain, y = HS)) +
  stat_summary(data = data.bis %>%
                 filter(Condition == 'BiGGY', orf_name %in% c('yll','met5','met10')) %>%
                 group_by(Condition, Strain) %>%
                 summarize(HS = mean(HS, na.rm = T), .groups = 'keep'),
               aes(fill = round(HS)), col = 'white', alpha = 0.9, size = 1,
               fun = mean, geom = "bar") +
  # stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  # scale_x_discrete(labels = c('BiGGY' = 'BiGGY', 'SD-Met-Cys+Bi' = 'SD-Met+Bi'),
  #                  limits = c('BiGGY', 'SD-Met-Cys+Bi')) +
  scale_y_continuous(breaks = seq(0,10,1)) +
  scale_fill_gradient(low = "#D7CCC8", high = "#5D4037", guide = F) +
  labs(y = 'Relative Hydrogen Sulfide',
       x = 'Strains',
       title = 'BiGGY Media') +
  scale_x_discrete(limits = c('FY4-met5D','FY4-met10D','FY4-yllD'),
                   labels = c('FY4'='FY4',
                              'FY4-met12D'='FY4-*met12Δ*',
                              'FY4-str3D'='FY4-*str3Δ*',
                              'FY4-met3D'='FY4-*met3Δ*',
                              'FY4-met15D'='FY4-*met15Δ*',
                              'FY4-met2D'='FY4-*met2Δ*',
                              'FY4-met6D'='FY4-*met6Δ*',
                              'FY4-met13D'='FY4-*met13Δ*',
                              'FY4-cys4D'='FY4-*cys4Δ*',
                              'FY4-yllD'='FY4-*yll058wΔ*',
                              'BY4742'='BY4742',
                              'BY4741'='BY4741',
                              'FY4-met5D' = 'FY4-*met5Δ*',
                              'FY4-met10D' = 'FY4-*met10Δ*')) +
  theme_linedraw() +
  theme(plot.title = element_blank(),#element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'black',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(1,7))

figSYLL <- plot_grid(figSjmYLL, figSjmYLLcp, figSjmYLLcol,
                     ncol = 1, rel_heights = c(1,1,1),
                     labels = c('A','B','C'),
                     label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')

ggsave(sprintf("%s/FigureSYLLBiGGY.jpg",fig_path), figSYLL,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)

# figSjmYLL.g <- ggplot_gtable(ggplot_build(figSjmYLL))
# stripr <- which(grepl('strip-t', figSjmYLL.g$layout$name))
# 
# for (s in stripr) {
#   # fig3b.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
#   l <- figSjmYLL.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
#   
#   if (length(l) != 0) {
#     if (l == 'FY4-') {
#       l2 <- figSjmYLL.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$label
#       a <- strain.labs.jm$auxotrophy[str_detect(strain.labs.jm$labels, l2)]
#     } else {
#       a <- strain.labs.jm$auxotrophy[strain.labs.jm$orf_name == l]
#     }
#     
#     if (a == 'Presumed Auxotroph') {
#       figSjmYLL.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
#     } else if (a == 'Prototroph') {
#       figSjmYLL.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
#     } else {
#       figSjmYLL.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107' #Unkown = '#BDBDBD'
#     }
#   }
# }

ggsave(sprintf("%s/FigureS9.jpg",fig_path), figSjmYLL,
       height = 70, width = two.c, units = 'mm',
       dpi = 600)


##### MET5/10 JAKE'S MUTANT
data.510$orf_name <- factor(data.510$orf_name, levels = c('FY4','met15del','met10del','met5del',
                                                          'BY4741_met10del','BY4741_met5del'))

figSmet510 <- merge(data.510, strain.labs.510, by = 'orf_name') %>%
  filter(hours == 160, condition %in% c('SD+Met','SD-Met','SD-Ura'),
         orf_name %in% c('FY4','met15del','BY4741_met5del','BY4741_met10del')) %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  geom_boxplot(fill = '#9E9E9E', outlier.shape = NA, size = 0.2,
               position = position_dodge(width = 0.85)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(limits = c('SD+Met',
                              'SD-Met',
                              'SD-Ura'),
                   labels = c('SD+Met' = '+ Met\n+ Ura\n+ Leu',
                              'SD-Met' = '- Met\n+ Ura\n+ Leu',
                              'SD-Ura' = '+ Met\n- Ura\n+ Leu')) +
  labs(title = 'Glucose',
       x = 'Strain',
       y = 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_wrap(~orf_name, nrow = 1, 
             labeller = labeller(orf_name = c('FY4'='FY4',
                                              'met15del'='FY4-*met15Δ*',
                                              'met5del'='FY4-*met5Δ*',
                                              'met10del'='FY4-*met10Δ*',
                                              'BY4741_met5del'='BY4741-*met5Δ*',
                                              'BY4741_met10del'='BY4741-*met10Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/FigureS10.jpg",fig_path), figSmet510,
       height = 70, width = two.c, units = 'mm',
       dpi = 600)


