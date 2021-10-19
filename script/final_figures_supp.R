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
figScbn1 <- data.cbn[str_detect(data.cbn$orf_name, 'FY'),] %>%
  filter(carbon == 'Glucose', base == 'SD') %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  # geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(fill = '#9E9E9E', size = 0.2, outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(labels = c('PlM_Glu' = '+ Met\n+ Ura',
                              'MiM_Glu' = '- Met\n+ Ura',
                              'MiU_Glu' = '+ Met\n- Ura',
                              'MiM_Gal' = '- Met\n+ Ura',
                              'MiU_Gal' = '+ Met\n- Ura',
                              'MiM_Et' = '- Met\n+ Ura',
                              'MiU_Et' = '+ Met\n- Ura')) +
  labs(x = 'Strain', y= 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_grid(orf_name ~ carbon, scales = 'free_x') +
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
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")),
        strip.text.y = element_text(color = 'black'))

figScbn1.g <- ggplot_gtable(ggplot_build(figScbn1))
stripr <- which(grepl('strip-r', figScbn1.g$layout$name))

for (s in stripr) {
  l <- figScbn1.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$label
  a <- strain.labs.cbn$auxotrophy[strain.labs.cbn$orf_name == l]
  
  if (a == 'Methionine') {
    figScbn1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
  } else if (a == 'None') {
    figScbn1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
  } else if (a == 'Uracil') {
    figScbn1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#E040FB'
  } else {
    figScbn1.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF5722'
  }
}
grid.draw(figScbn1.g)


figScbn2 <- data.cbn[str_detect(data.cbn$orf_name, 'FY'),]  %>%
  filter(carbon != 'Glucose', base == 'SD') %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  # geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(fill = '#9E9E9E', size = 0.2, outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(labels = c('PlM_Glu' = '+ Met\n+ Ura',
                              'MiM_Glu' = '- Met\n+ Ura',
                              'MiU_Glu' = '+ Met\n- Ura',
                              'MiM_Gal' = '- Met\n+ Ura',
                              'MiU_Gal' = '+ Met\n- Ura',
                              'MiM_Et' = '- Met\n+ Ura',
                              'MiU_Et' = '+ Met\n- Ura')) +
  labs(x = 'Strain', y= 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_grid(orf_name ~ carbon, scales = 'free_x') +
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
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")),
        strip.text.y = element_text(color = 'black'))

figScbn2.g <- ggplot_gtable(ggplot_build(figScbn2))
stripr <- which(grepl('strip-r', figScbn2.g$layout$name))

for (s in stripr) {
  l <- figScbn2.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$label
  a <- strain.labs.cbn$auxotrophy[strain.labs.cbn$orf_name == l]
  
  if (a == 'Methionine') {
    figScbn2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
  } else if (a == 'None') {
    figScbn2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
  } else if (a == 'Uracil') {
    figScbn2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#E040FB'
  } else {
    figScbn2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF5722'
  }
}
grid.draw(figScbn2.g)


figScbn_cols <- readPNG('figures/final/final/Figure1C_Supp.png')
figScbn_cols <- ggplot() + 
  background_image(figScbn_cols) +
  theme(plot.margin = margin(t=0, l=30, r=30, b=0, unit = "mm"),
        plot.background = element_blank())


figScbnSD <- annotate_figure(cowplot::plot_grid(figScbn_cols,
                                                legend.aux,
                                                cowplot::plot_grid(figScbn1.g, figScbn2.g,
                                                                   rel_widths = c(3,4), labels = c('B','C'), nrow = 1,
                                                                   label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                                                labels = c('A','',''),
                                                label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
                                                ncol = 1, rel_heights = c(5,1,7)),
                             top = text_grob('Synthetic Defined Media',
                                             face = "bold", family = "sans", size = titles))
ggsave(sprintf("%s/FigureS2.jpg",fig_path), figScbnSD,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)


##### CARBON SOURCE - SC
figScbn3 <- data.cbn %>%
  filter(base == 'SC', orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  # geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(aes(fill = auxotrophy), size = 0.2, outlier.shape = NA,
               position = position_dodge(width = 0.6)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(limits = c('BY4742','BY4741',
                              'FY4','FY4-met15del'),
                   labels = c( 'FY4' = 'FY4',
                               'FY4-met15del' = 'FY4-*met15Δ*',
                               'BY4742' = 'BY4742',
                               'BY4741' = 'BY4741')) +
  scale_fill_manual(name = 'Methionine-Uracil\nAuxotrophy',
                    values = c('None' = '#FFC107',
                               'Methionine' = '#536DFE',
                               'Uracil' = '#E040FB',
                               'Both' = '#FF5722'),
                    position = 'bottom') +
  labs(x = 'Strain', y= 'Relative Colony Size') +
  facet_grid(.~carbon) +
  coord_cartesian(ylim = c(0,1.3)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(fill = guide_legend(nrow=2, byrow=TRUE, order = 1))


figScbnSC <- annotate_figure(cowplot::plot_grid(legend.aux, figScbn3 + theme(legend.position = 'none'),
                                              ncol = 1, rel_heights = c(1,6)),
                           top = text_grob('Synthetic Complete Media',
                                           face = "bold", family = "sans", size = titles))

ggsave(sprintf("%s/FigureS3.jpg",fig_path), figScbnSC,
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
  facet_grid(~carbon*methionine, labeller = labeller(methionine = c('+Met' = 'SD + Met',
                                                                   '-Met' = 'SD - Met'))) +
  scale_color_manual(name = 'Methionine\nAuxotrophy',
                     values = c('Prototroph' = '#FFC107',
                                'Presumed Auxotroph' = '#536DFE'),
                     limits = c('Prototroph','Presumed Auxotroph'),
                     guide = F) +
  scale_fill_manual(name = 'Methionine\nAuxotrophy',
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
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(linetype = guide_legend(nrow=1, byrow=TRUE, order = 1))


figSres2 <- readPNG('figures/final/final/RESP_STATS_SUPP.png')
figSres2 <- ggplot() + 
  background_image(figSres2) +
  theme(plot.margin = margin(t=0, l=10, r=10, b=0, unit = "mm"),
        plot.background = element_blank())

legend.pet <- g_legend(figSres1)
figSres <- cowplot::plot_grid(figSres1 + theme(legend.position = 'none'),
                              cowplot::plot_grid(legend.aux, legend.pet, nrow = 1, rel_widths = c(3,1)),
                              figSres2,
                              ncol = 1, rel_heights = c(6,1,6))
ggsave(sprintf("%s/FigureS4.jpg",fig_path), figSres,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)


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
data.bis$Strain <- factor(data.bis$Strain,
                          levels = c('FY4','BY4742','FY4-met12D','FY4-met2D','FY4-met6D','FY4-met13D',
                                     'FY4-cys4D','FY4-yllD','FY4-str3D','FY4-met3D','FY4-met15D','BY4741'))
figSBi <- data.bis %>%
  filter(Condition == 'SD-Met-Cys+Bi') %>%
  ggplot(aes(x = 1, y = HS)) +
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
  facet_wrap(.~Strain, ncol = 3,
             labeller = labeller(Strain = c('FY4'='FY4',
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
                                            'BY4741'='BY4741'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'black',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(1,7),
                  xlim = c(0.2,1.8))

figSBi.g <- ggplot_gtable(ggplot_build(figSBi))
stripr <- which(grepl('strip-t', figSBi.g$layout$name))

for (s in stripr) {
  # figSBi.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
  l <- figSBi.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
  
  if (length(l) != 0) {
    if (l == 'FY4-') {
      l2 <- figSBi.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$label
      a <- strain.labs.bis$auxotrophy[str_detect(strain.labs.bis$labels, l2)]
    } else {
      a <- strain.labs.bis$auxotrophy[strain.labs.bis$orf_name == l]
    }
    
    if (a == 'Presumed Auxotroph') {
      figSBi.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
    } else if (a == 'Prototroph') {
      figSBi.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
    } else {
      figSBi.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#BDBDBD'
    }
  }
}

ggsave(sprintf("%s/FigureS7.jpg",fig_path), figSBi.g,
       height = one.5c, width = one.5c, units = 'mm',
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
      a <- strain.labs.ns$auxotrophy[str_detect(strain.labs.ns$labels, l2)]
    }
    if (a == 'Methionine') {
      figSns2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
    } else if (a == 'None') {
      figSns2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
    } else if (a == 'Uracil') {
      figSns2.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#E040FB'
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
