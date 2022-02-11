##### FINAL FIGURES FOR METHIONINE PAPER
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 09/14/2021 

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
library(cowplot)
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
# load(file = 'figures/final/data.RData')

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
  # scale_fill_manual(name = 'Presumed\nAuxotrophy',
  #                   limits = c('None',
  #                              'Methionine',
  #                              'Uracil-Leucine',
  #                              'Methionine-Uracil-Leucine'),
  #                   values = c('None' = '#FFC107',
  #                              'Methionine' = '#536DFE',
  #                              'Uracil-Leucine' = '#E040FB',
  #                              'Methionine-Uracil-Leucine' = '#FF5722',
  #                              'Unknown' = '#BDBDBD'),
  #                   position = 'bottom') +
  scale_fill_manual(name = 'Presumed\nAuxotroph',
                  limits = c('None',
                             'Methionine'),
                  values = c('None' = '#FFC107',
                             'Methionine' = '#536DFE'),
                  labels = c('None' = 'No',
                             'Methionine' = 'Yes'),
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
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(fill = guide_legend(nrow=2, byrow=F, order = 1))

legend.aux <- g_legend(legend.aux)

ggsave(sprintf("%s/ColorLegend.jpg",fig_path), legend.aux,
       height = one.c, width = one.c, units = 'mm',
       dpi = 600)

##### FIGURE 1. INTRODUCTION TO METHIONINE PATHWAY AND EXPT PROTOCOL
fig1a <- readPNG('figures/final/final/Figure1A.png')
fig1a <- ggplot() + 
  background_image(fig1a) +
  theme(plot.margin = margin(t=0, l=20, r=20, b=0, unit = "mm"),
        plot.background = element_blank())

fig1b <- readPNG('figures/final/final/ExptPipeline.png')
fig1b <- ggplot() + 
  background_image(fig1b) +
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
        plot.background = element_blank())

fig1c <- readPNG('figures/final/final/Figure1C.png')
fig1c <- ggplot() + 
  background_image(fig1c) +
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
        plot.background = element_blank())

fig1 <- cowplot::plot_grid(fig1a, cowplot::plot_grid(fig1b, fig1c, nrow = 1, rel_widths = c(5,2),
                                                     labels = c('B','C'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                           labels = c('A',''), ncol = 1, rel_heights = c(1.2,1),
                           label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure1.jpg",fig_path), fig1,
       height = 170, width = two.c, units = 'mm',
       dpi = 600)

##### FIGURE 2. CARBON SOURCE & RESPIRATION
data.cbn$methionine <- factor(data.cbn$methionine, levels = c('+Met +Ura', '+Met -Ura', '-Met +Ura'))
data.cbn$condition <- factor(data.cbn$condition, levels = c('PlM_Glu', 'MiM_Glu', 'MiU_Glu',
                                                            'MiM_Gal', 'MiU_Gal',
                                                            'MiM_Et', 'MiU_Et',
                                                            'GLU', 'GAL'))
fig2a <- rbind(data.cbn[,c(1,2,5,6,10)],
               data.cbn.leu[,c(8,11,10,13,14)] %>%
                 filter(condition == 'SDmLeu')) %>%
  filter(carbon == 'Glucose', base == 'SD', orf_name %in% c('BY4741','BY4742')) %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  # geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(fill = '#9E9E9E', size = 0.2, outlier.shape = NA) +
  # stat_compare_means(method = 't.test',
  #                    comparisons = list(c('MiM_Glu','MiU_Glu'), c('MiM_Glu','SDmLeu'))) +
  scale_y_continuous(breaks = seq(-1.2,2,0.3)) +
  scale_x_discrete(limits = c('PlM_Glu', 'MiM_Glu', 'MiU_Glu', 'SDmLeu'),
                   labels = c('PlM_Glu' = '+ Met\n+ Ura\n+ Leu',
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
  coord_cartesian(ylim = c(0,1.6)) +
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
                                              margin = margin(1,0,0.1,0, "mm")),
        strip.background.y = element_blank(),
        strip.text.y = element_blank())

# fig2a.g <- ggplot_gtable(ggplot_build(fig2a))
# stripr <- which(grepl('strip-r', fig2a.g$layout$name))
# 
# fig2a.g$grobs[[11]]$grobs[[1]]$children[[2]]$children[[1]]$label
# fig2a.g$grobs[[11]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
# fig2a.g$grobs[[11]]$grobs[[1]]$children[[1]]$gp$fill <- '#E040FB'
# 
# fig2a.g$grobs[[12]]$grobs[[1]]$children[[2]]$children[[1]]$label
# fig2a.g$grobs[[12]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
# fig2a.g$grobs[[12]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF5722'
# 
# ggsave(sprintf("%s/CarbonSource_SD_Glu.jpg",fig_path), fig2a.g,
#        height = one.c, width = one.c, units = 'mm',
#        dpi = 600)

# grid.draw(fig2a.g)

fig2b <- data.cbn %>%
  filter(carbon != 'Glucose', base == 'SD', orf_name %in% c('BY4741','BY4742')) %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  # geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(fill = '#9E9E9E', size = 0.2, outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-1.2,2,0.3)) +
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
  coord_cartesian(ylim = c(0,1.6)) +
  facet_grid(orf_name ~ carbon, scales = 'free_x') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                                margin = margin(1,0,0.1,0, "mm")),
        strip.text.y = ggtext::element_markdown(size = txt, face = 'bold',
                                              margin = margin(0.1,1,0.1,0, "mm")))

fig2b.g <- ggplot_gtable(ggplot_build(fig2b))
stripr <- which(grepl('strip-r', fig2b.g$layout$name))

fig2b.g$grobs[[16]]$grobs[[1]]$children[[2]]$children[[1]]$label
fig2b.g$grobs[[16]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
fig2b.g$grobs[[16]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF0000'#'#E040FB'

fig2b.g$grobs[[17]]$grobs[[1]]$children[[2]]$children[[1]]$label
fig2b.g$grobs[[17]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
fig2b.g$grobs[[17]]$grobs[[1]]$children[[1]]$gp$fill <- '#990000'#'#FF5722'

fig2b.leg <- g_legend(data.frame(labels = c('Uracil-Leucine','Methionine-Uracil-Leucine')) %>%
  ggplot(aes(x = 1, y = labels, col = labels)) +
  geom_point(shape = 15) +
  scale_color_manual(name = 'Presumed Auxotrophy',
                     values = c('Uracil-Leucine' = '#FF0000',
                                'Methionine-Uracil-Leucine' = '#990000')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title =  ggtext::element_markdown(size = titles),
        legend.text =  ggtext::element_markdown(size = txt),
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(4, "mm"),
        legend.direction = 'horizontal',
        legend.position = 'bottom') +
  guides(col = guide_legend(nrow=1, override.aes=list(size = 3))))
  

# ggsave(sprintf("%s/CarbonSource_SD_GalEth.jpg",fig_path), fig2b.g,
#        height = one.c, width = one.c, units = 'mm',
#        dpi = 600)

# grid.draw(fig2b.g)
strain.labs.res$labels2 <- as.character(strain.labs.res$labels)
# strain.labs.res$labels2[str_detect(strain.labs.res$labels2, '(rho)')] <-
#   paste0(strain.labs.res$labels2[str_detect(strain.labs.res$labels2, '(rho)')], '')

fig2c <- merge(data.res.gc[str_detect(data.res.gc$orf_name, 'FY'),],
              strain.labs.res, by = 'orf_name') %>%
  filter(expt_rep %in% c("1","2"), condition != 'SD+Met-Ura+Glu') %>% 
  ggplot(aes(x = Time, y = rel_cs)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               aes(group = orf_name, fill = met_aux), geom="ribbon", alpha = 0.4) +
  stat_summary(aes(group = orf_name, col = met_aux, linetype = pet), fun=mean, geom="line", lwd =0.7) +
  geom_text_repel(data = merge(data.res.gc[str_detect(data.res.gc$orf_name, 'FY') &
                                             data.res.gc$Time == max(data.res.gc$Time),],
                               strain.labs.res, by = 'orf_name') %>%
                    filter(expt_rep %in% c("1","2"), condition != 'SD+Met-Ura+Glu') %>%
                    group_by(orf_name, condition, Time, base, methionine, carbon, cysteine, labels, labels2, met_aux, pet) %>%
                    summarize(rel_cs = mean(rel_cs, na.rm = T), .groups = 'keep'),
                  aes(x = Time, y = rel_cs, label = labels2),
                  parse = T, size = 2, min.segment.length = 10) +
  facet_grid(~carbon*methionine,
             labeller = labeller(methionine = c('+Met' = 'SD + Met',
                                                '-Met' = 'SD - Met'))) +
  scale_color_manual(name = 'Presumed Auxotroph',
                     values = c('Prototroph' = '#FFC107',
                                'Presumed Auxotroph' = '#536DFE',
                                'Uracil-Leucine' = '#E040FB',
                                'Methionine-Uracil-Leucine' = '#FF5722'),
                     limits = c('Prototroph','Presumed Auxotroph'),
                     labels = c('Prototroph'='No',
                                'Presumed Auxotroph'='Yes',
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
  guides(color = guide_legend(nrow=1, byrow=TRUE, order = 1,
                              override.aes = list(size = 3)),
         linetype = guide_legend(nrow=1, byrow=TRUE, order = 2))

ggsave(sprintf("%s/Respiration_FYs.jpg",fig_path), fig2c,
       height = one.c, width = two.c, units = 'mm',
       dpi = 600)

# legend.pet <- g_legend(fig2c)
# fig2 <- cowplot::plot_grid(legend.aux,
#                            cowplot::plot_grid(fig2a.g, fig2b.g,
#                                               rel_widths = c(3,4), labels = c('A','B'), nrow = 1,
#                                               label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
#                            fig2c + theme(legend.position="bottom", legend.direction = 'horizontal'),
#                            labels = c('','','C'), ncol = 1, rel_heights = c(0.1,1,1),
#                            label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure2.jpg",fig_path), fig2,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)


##### FIGURE 3. JAKES MUTANT, REPEATED PINNING & NO SULFATES
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

data.rp %>%
  group_by(carbon, pin) %>%
  summarize(max_hrs = max(hours, na.rm = T)) %>%
  data.frame()

data.rp.rf <- merge(data.rp %>%
        filter(aux == 'Met', orf_name %in% c('FY4'), carbon == 'Glu', hours == max_hrs) %>%
        group_by(pin, orf_name, hours, cum_hrs) %>%
        summarize(average = median(average, na.rm = T), .groups = 'keep'),
      
      data.rp %>%
        filter(aux == 'Met', orf_name %in% c('FY4-met15D'), carbon == 'Glu', hours == max_hrs) %>%
        group_by(pin, orf_name, hours, cum_hrs) %>%
        summarize(average = median(average, na.rm = T), .groups = 'keep'),
      by = c('pin','hours','cum_hrs'), suffixes = c('_ref',''))
data.rp.rf$relative_fitness <- data.rp.rf$average/data.rp.rf$average_ref

fig3a <- merge(data.rp, strain.labs.rp,
               by = 'orf_name') %>%
  filter(aux == 'Met', orf_name %in% c('FY4'), carbon == 'Glu') %>%
  ggplot(aes(x = cum_hrs, y = average)) +
  stat_summary(aes(fill = auxotrophy, group = pin),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(aes(col = auxotrophy, group = pin),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Met', orf_name %in% c("FY4-met15D"), carbon == 'Glu'),
               aes(fill = auxotrophy, group = pin),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = merge(data.rp, strain.labs.rp,
                            by = 'orf_name') %>%
                 filter(aux == 'Met', orf_name %in% c("FY4-met15D"), carbon == 'Glu'),
               aes(col = auxotrophy, group = pin),
               fun=mean, geom="line", lwd =1) +
  geom_text_repel(data = merge(data.rp, strain.labs.rp,
                               by = 'orf_name') %>%
                    filter(aux == 'Met', orf_name %in% c('FY4',"FY4-met15D"), carbon == 'Glu', hours == max_hrs, pin == 7) %>%
                    group_by(orf_name, labels, carbon, cum_hrs) %>%
                    summarize(average = mean(average, na.rm = T), .groups = 'keep'),
                  aes(x = cum_hrs, y = average + 300, label = labels),
                  parse = T, size = 2.5) +
  geom_text(data = data.rp.rf,
                  aes(x = cum_hrs + 10, y = average, 
                      label = sprintf('%0.2f',relative_fitness)),
                  size = 2.5) +
  scale_x_continuous(breaks = seq(-1000,1000,100)) +
  scale_color_manual(name = 'Auxotrophy',
                     values = c('Methionine' = '#536DFE',
                                'None' = '#FFC107')) +
  scale_fill_manual(name = 'Auxotrophy',
                     values = c('Methionine' = '#536DFE',
                                'None' = '#FFC107')) +
  labs(x = 'Time (hours)', y = 'Colony Size (pixels)') +
  facet_grid(.~carbon, labeller = labeller(carbon = c('Glu'='SD-Met+Glu'))) +
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
                                  margin = margin(1,0,0.1,0, "mm")),
        strip.text.y = element_text(size = txt, margin = margin(0.1,0,0.1,0, "mm")))

  
# fig3a.g <- ggplot_gtable(ggplot_build(fig3a))
# stripr <- which(grepl('strip-t', fig3a.g$layout$name))
# 
# for (s in stripr) {
#   # fig3a.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
#   l <- fig3a.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
#   
#   if (l == 'FY4-') {
#     l2 <- fig3a.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$label
#     a <- strain.labs.rp$auxotrophy[str_detect(strain.labs.rp$labels, l2)]
#   } else {
#     a <- strain.labs.rp$auxotrophy[strain.labs.rp$orf_name == l]
#   }
#   
#   if (a == 'Methionine') {
#     fig3a.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
#   } else if (a == 'None') {
#     fig3a.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
#   } else if (a == 'Uracil') {
#     fig3a.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#E040FB'
#   } else {
#     fig3a.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF5722'
#   }
# }
# # grid.draw(fig3a.g)

data.bis$orf_name[data.bis$Strain == 'FY4'] <- 'FY4'
data.bis$orf_name[data.bis$Strain == 'FY4-met15D'] <- 'met15'
data.bis$orf_name[data.bis$Strain == 'FY4-met3D'] <- 'met3'
data.bis$orf_name[data.bis$Strain == 'FY4-met13D'] <- 'met13'
data.bis$orf_name[data.bis$Strain == 'FY4-met12D'] <- 'met12'
data.bis$orf_name[data.bis$Strain == 'FY4-met6D'] <- 'met6'
data.bis$orf_name[data.bis$Strain == 'FY4-met2D'] <- 'met2'
data.bis$orf_name[data.bis$Strain == 'FY4-str3D'] <- 'str3'
data.bis$orf_name[data.bis$Strain == 'FY4-cys4D'] <- 'cys4'
data.bis$orf_name[data.bis$Strain == 'FY4-met15D'] <- 'met15'
data.bis$orf_name[data.bis$Strain == 'BY4741'] <- 'BY4741'
data.bis$orf_name[data.bis$Strain == 'BY4742'] <- 'BY4742'



data.jm$orf_name <- factor(data.jm$orf_name,
                           levels = c('FY4','BY4742','met12','met2','met6','met13',
                                      'cys4','yll','str3','met3','met15','BY4741'))
fig3b <- merge(merge(data.jm, strain.labs.jm, by = 'orf_name'),
      data.bis %>% filter(Condition == 'BiGGY') %>%
        group_by(orf_name) %>%
        summarize(HS = median(HS, na.rm = T), .groups = 'keep') %>%
        data.frame(), by = 'orf_name') %>%
  filter(orf_name != 'yll', time == 't_final', attempt != 'pilot', condition %in% c('YPDA','SD-Met-Cys+Glu')) %>%
  ggplot(aes(x = condition, y = relative_cs, alpha = condition)) +
  # geom_rect(aes(fill = HS),xmin = -Inf,xmax = Inf,
  #           ymin = -Inf,ymax = Inf) +
  # geom_hline(yintercept = seq(-1,2,0.2), col = '#FFECB3', size = 0.2, linetype = 'dashed') +
  # geom_vline(xintercept = seq(-1,5,1), col = '#FFECB3', size = 0.2, linetype = 'dashed') +
  # geom_boxplot(fill = '#9E9E9E', col = 'white', outlier.shape = NA, size = 0.4,
  #              position = position_dodge(width = 0.85)) +
  geom_boxplot(fill = '#9E9E9E', outlier.shape = NA, size = 0.2,
               position = position_dodge(width = 0.85)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(limits = c('YPDA', 'SD-Met-Cys+Glu'),
                   labels = c('YPDA' = 'YPDA',
                              'SD-Met-Cys+Glu' = 'SD-Met+Glu')) +
  scale_alpha_manual(name = 'Condition',
                     values = c('YPDA' = 1,
                                'SD-Met-Cys+Glu' = 1),
                     labels = c('YPDA' = 'YPDA',
                                'SD-Met-Cys+Glu' = 'SD-Met-Cys+Ura+Glu'),
                     guide = F) +
  scale_fill_gradient(name = 'Relative Hydrogen Sulfide',
                      low = "#D7CCC8", high = "#5D4037") +
  # scale_fill_manual(name = 'Methionine\nAuxotrophy',
  #                   values = c('Prototroph' = '#FFC107',
  #                              'Presumed Auxotroph' = '#536DFE',
  #                              'Unknown' = '#616161'),
  #                   position = 'bottom') +
  labs(x = 'Strain',
       y = 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_wrap(~orf_name, nrow = 2, 
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
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        # legend.direction = 'vertical',
        # legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'black',
                                              margin = margin(0.1,0,0.1,0, "mm")))

fig3b.g <- ggplot_gtable(ggplot_build(fig3b))
stripr <- which(grepl('strip-t', fig3b.g$layout$name))

for (s in stripr) {
  # fig3b.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
  l <- fig3b.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
  
  if (length(l) != 0) {
    if (l == 'FY4-') {
      l2 <- fig3b.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$label
      a <- strain.labs.jm$auxotrophy[str_detect(strain.labs.jm$labels, l2)]
    } else {
      a <- strain.labs.jm$auxotrophy[strain.labs.jm$orf_name == l]
    }
    
    if (a == 'Presumed Auxotroph') {
      fig3b.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
    } else if (a == 'Prototroph') {
      fig3b.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
    } else {
      fig3b.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107' #Unkown = '#BDBDBD'
    }
  }
}
# grid.draw(fig3b.g)
# ggsave(sprintf("%s/JakesMutants.jpg",fig_path), fig3b.g,
#        height = one.5c, width = one.5c, units = 'mm',
#        dpi = 600)


data.bis$Strain <- factor(data.bis$Strain,
                           levels = c('FY4','BY4742','FY4-met12D','FY4-met2D','FY4-met6D','FY4-met13D',
                                      'FY4-cys4D','FY4-yllD','FY4-str3D','FY4-met3D','FY4-met15D','BY4741'))
fig3c <- data.bis %>%
  filter(Condition == 'BiGGY') %>%
  ggplot(aes(x = Strain, y = HS)) +
  stat_summary(data = data.bis %>%
                 filter(Condition == 'BiGGY') %>%
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
  scale_x_discrete(labels = c('FY4'='FY4',
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
                                   'BY4741'='BY4741')) +
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

# fig3c.g <- ggplot_gtable(ggplot_build(fig3c))
# stripr <- which(grepl('strip-t', fig3c.g$layout$name))
# 
# for (s in stripr) {
#   # fig3c.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
#   l <- fig3c.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
# 
#   if (length(l) != 0) {
#     if (l == 'FY4-') {
#       l2 <- fig3c.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$label
#       a <- strain.labs.bis$auxotrophy[str_detect(strain.labs.bis$labels, l2)]
#     } else {
#       a <- strain.labs.bis$auxotrophy[strain.labs.bis$orf_name == l]
#     }
# 
#     if (a == 'Presumed Auxotroph') {
#       fig3c.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
#     } else if (a == 'Prototroph') {
#       fig3c.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
#     } else {
#       fig3c.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#BDBDBD'
#     }
#   }
# }
# grid.draw(fig3c.g)

fig3cp <- readPNG('figures/final/final/BiGGY_Main.png')
fig3cp <- ggplot() + 
  background_image(fig3cp) +
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
        plot.background = element_blank())

strain.labs.ns <- strain.labs.cbn
strain.labs.ns$orf_name[strain.labs.ns$orf_name == 'FY4-met3del'] <- 'FY4_met3del'
strain.labs.ns$orf_name[strain.labs.ns$orf_name == 'FY4-met15del'] <- 'FY4_met15del'
data.ns$orf_name <- factor(data.ns$orf_name, levels = c('FY4','FY4_met15del','FY4_met3del',
                                                        'BY4742','BY4741'))
fig3d <- merge(data.ns[data.ns$stage %in% c('S1','S2') & data.ns$hours %in% c(165,149) & data.ns$average != 0 &
                str_detect(data.ns$orf_name, 'FY'),],
               strain.labs.ns, by = 'orf_name') %>%
  filter(base == 'Agarose', ynb_type == 'Home') %>%
  ggplot(aes(x = stage, y = average)) +
  geom_boxplot(size = 0.2, outlier.shape = NA, fill = '#9E9E9E',
               position = position_dodge(width = 0.85)) +
  labs(y = 'Colony Size (pixels)') +
  scale_x_discrete(labels = c( 'S1' = 'Pin #1',
                               'S2' = 'Pin #2')) +
  facet_grid(sulfate~orf_name,
             labeller = labeller(orf_name = c('FY4' = 'FY4',
                                              'FY4_met15del' = 'FY4-*met15Δ*',
                                              'FY4_met3del' = 'FY4-*met3Δ*',
                                              'BY4742' = 'BY4742',
                                              'BY4741' = 'BY4741'),
                                 sulfate = c('+sulfate' = 'SD - Met + Glu w/<br /> Inorganic Sulfates',
                                             '-sulfate' = 'SD - Met + Glu w/o<br /> Inorganic Sulfates'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.box = 'vertical',
        # legend.direction = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        legend.spacing.y = unit(0.5, "mm"),
        strip.text = ggtext::element_markdown(size = txt, color = 'black',
                                              margin = margin(0.1,0,0.1,0, "mm")),
        strip.text.y = ggtext::element_markdown(color = 'white'))


fig3d.g <- ggplot_gtable(ggplot_build(fig3d + theme(legend.position = 'None')))
stripr <- which(grepl('strip-t', fig3d.g$layout$name))

for (s in stripr) {
  # fig3d.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
  l <- fig3d.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
  
  if (str_detect(l, 'FY4-')) {
    a <- 'Methionine'
  } else {
    a <- strain.labs.ns$auxotrophy[strain.labs.ns$labels == l]
  }
  
  if (a == 'Methionine') {
    fig3d.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
  } else if (a == 'None') {
    fig3d.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
  } else if (a == 'Uracil') {
    fig3d.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#E040FB'
  } else {
    fig3d.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF5722'
  }
}
# grid.draw(fig3d.g)

# fig3 <- cowplot::plot_grid(legend.aux,
#                            fig3a.g,
#                            cowplot::plot_grid(fig3b.g, fig3c.g,
#                                               ncol = 2, rel_widths = c(2,1.5), align = 'h',
#                                               labels = c('B','C'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
#                            cowplot::plot_grid(fig3d.g, legend.sul, nrow = 1, rel_widths = c(5,1)),
#                            ncol = 1, rel_heights = c(0.2,1,2,1),
#                            labels = c('','A','','D'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
# ggsave(sprintf("%s/Figure3.jpg",fig_path), fig3,
#        height = 220, width = two.c, units = 'mm',
#        dpi = 600)

ggsave(sprintf("%s/NoSulfate.jpg",fig_path), fig3d.g,
       height = one.c, width = two.c, units = 'mm',
       dpi = 600)


###### SCREEN RESULTS
fig4a <- merge(data.del.diff %>% filter(orf_name %in% strain.labs.del$orf_name),
      strain.labs.del, by = 'orf_name') %>%
  ggplot(aes(x = fitness_MM, y = fitness_PM)) +
  geom_point(data = data.del.diff, size = 1, col = '#9E9E9E') +
  # geom_point(size = 2) +
  geom_line(data = data.del.diff.dist, aes(x = x , y = y.ul),
            linetype = 'dashed', size = 0.5) +
  geom_line(data = data.del.diff.dist, aes(x = x , y = y.ll),
            linetype = 'dashed', size = 0.5) +
  geom_point(aes(x = fitness_MM, y = fitness_PM), shape = 1, size = 2) +
  geom_point(data = data.del.diff %>% filter(orf_name %in%
                                               strain.labs.del$orf_name[strain.labs.del$standard_name %in% 
                                                                          c('YLL058W','MET12','MET5','MET10')]), 
                                             size = 1, col = '#FFC107') +
  geom_text_repel(aes(x = fitness_MM, y = fitness_PM, label = standard_name), size = 2,
                  min.segment.length = unit(0, 'lines'), seed = 10,
                  force = 2, max.overlaps = 30) +
  scale_x_continuous(trans = 'pseudo_log') +
  labs(x = 'Relative Fitness in SD-Met+Gal (log)',
       y = 'Relative Fitness in SD+Met+Gal') +
  coord_cartesian(xlim = c(0, 15),
                  ylim = c(0, 1.5)) +
  facet_wrap(.~stage, ncol = 3,
             labeller = labeller(stage = c('Pre-Screen #1' = 'Pin #1',
                                                    'Pre-Screen #2' = 'Pin #2',
                                                    'Final Screen' = 'Pin #3'))) +
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
  guides(color = guide_legend(nrow=3, byrow=TRUE, order = 1))

head(data.del.diff)
data.del.diff %>%
  group_by(stage) %>%
  summarize(max_hrs = max(hours_MM, na.rm = T), .groups = 'keep') %>%
  data.frame()
# ggsave(sprintf("%s/DeletionScreen.jpg",fig_path), fig4a,
#        height = 80, width = two.c, units = 'mm',
#        dpi = 600)



fig4c <- data.pv2 %>%
  filter(orf_name %in% c('Plasmid_1','Plasmid_5','Plasmid_7','Plasmid_11'),
         arm == 'PV_FY_MM', stage != 'PS2') %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  geom_boxplot(outlier.shape = NA, fill = '#9E9E9E', size = 0.3) +
  scale_x_discrete(limits = c('Plasmid_1','Plasmid_5','Plasmid_7','Plasmid_11'),
                   labels = c('Plasmid_1' = 'FY4<br />w/ Empty Plasmid',
                              'Plasmid_5' = 'FY4-*met15Δ*<br />w/ Empty Plasmid',
                              'Plasmid_7' = 'FY4-*met15Δ yll058wΔ*<br />w/ Empty Plasmid',
                              'Plasmid_11' = 'FY4-*met15Δ*<br />w/ *YLL058W* Plasmid')) +
  labs(title = 'SD - Met + Gal + G418',
       y = 'Relative Colony Size') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        axis.text.x = ggtext::element_markdown(size = txt-1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%s/YLL_PlasmidValidation.jpg",fig_path), fig4c,
#        height = one.c, width = one.5c, units = 'mm',
#        dpi = 600)

fig4c.2 <- data.pv2 %>%
  filter(orf_name %in% c('Plasmid_2','Plasmid_6','Plasmid_12'),
         arm == 'PV_FY_MM', stage != 'PS2') %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  geom_boxplot(outlier.shape = NA, fill = '#9E9E9E', size = 0.3) +
  scale_x_discrete(limits = c('Plasmid_2','Plasmid_6','Plasmid_12'),
                   labels = c('Plasmid_2' = 'FY4<br />w/ Empty Plasmid',
                              'Plasmid_6' = 'FY4-*met15Δ*<br />w/ Empty Plasmid',
                              'Plasmid_12' = 'FY4-*met15Δ*<br />w/ *MET15* Plasmid')) +
  labs(title = 'SD - Met + Gal + Hyg',
       y = 'Relative Colony Size') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        axis.text.x = ggtext::element_markdown(size = txt-1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%s/MET15_PlasmidValidation.jpg",fig_path), fig4c.2,
#        height = one.c, width = one.5c, units = 'mm',
#        dpi = 600)


fig4c.3 <- data.pv %>%
  filter(hours %in% c(48,115), stage == 'PS1',
         deletion1 %in% c('met15', 'met12'), deletion2 %in% c('','met12'),
         plasmid_backbone == '2M', plasmid_orf %in% c('empty','met12')) %>%
  ggplot(aes(x = orf_name, y = average)) +
  geom_boxplot(fill = '#9E9E9E', outlier.shape = NA, size = 0.3) +
  labs(title = 'SD - Met + Gal + G418',
       y = 'Colony Size (pixels)') +
  scale_x_discrete(limits = c('met15del_2M_empty',
                              'met15del_2M_met12',
                              'met15del_met12del_2M_empty',
                              'met15del_met12del_2M_met12'),
                   labels = c('met15del_2M_empty' = 'FY4-*met15Δ*<br />w/ Empty Plasmid',
                              'met15del_2M_met12' = 'FY4-*met15Δ*<br />w/ *MET12* Plasmid',
                              'met15del_met12del_2M_empty' = 'FY4-*met15Δ met12Δ*<br />w/ Empty Plasmid',
                              'met15del_met12del_2M_met12'  = 'FY4-*met15Δ met12Δ*<br />w/ *MET12* Plasmid',
                              'met12del_2M_empty' = 'FY4-*met12Δ*<br />w/ Empty Plasmid')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt-1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,6000))

# ggsave(sprintf("%s/MET12_PlasmidValidation.jpg",fig_path), fig4c.3,
#        height = one.c, width = two.c, units = 'mm',
#        dpi = 600)


fig5d <- data.bioc %>%
  # filter(outlier == FALSE) %>%
  ggplot(aes(x = `Time (min)`, y = uM)) +
  stat_summary(aes(group = Sample, fill = Sample), fun.data=mean_se, fun.args = list(mult=1), geom="ribbon",
               alpha = 0.4) +
  stat_summary(aes(group = Sample, col = Sample), fun=mean, geom="line", lwd = 0.7) +
  # stat_summary(aes(group = Sample, col = Sample), fun.data = mean_se, geom = "errorbar", lwd = 0.7) +
  stat_summary(aes(group = Sample), fun=mean, geom="point", size =2) +
  stat_summary(aes(group = Sample, col = Sample), fun=mean, geom="point", size = 0.5) +
  labs(y = 'Homocysteine (log2(μM))') +
  scale_color_manual(name = 'Sample',
                     values = c('Met15' = "#7C4DFF",
                                'Yll058w' = "#FF5722",
                                'None' = "#757575"),
                     guide = F) +
  scale_fill_manual(name = 'Sample',
                     values = c('Met15' = "#7C4DFF",
                                'Yll058w' = "#FF5722",
                                'None' = "#757575")) +
  scale_y_continuous(trans = 'log2') +
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
  coord_cartesian(ylim = c(5,125)) +
  guides(fill = guide_legend(override.aes=list(shape = 15, alpha = 1)))

ggsave(sprintf("%s/Biochemistry.jpg",fig_path), fig5d,
       height = one.c, width = one.c, units = 'mm',
       dpi = 600)

fig5e <- data.mdl %>%
  filter(label %in% c('C. B + Hypothesized YLL058W reaction',
                        'D. C - All MET15 reactions')) %>%
  ggplot(aes(x = log(Yll058w.Met15,10), y = Growth)) +
  geom_line(aes(col = label), lwd = 1) +
  geom_point(size = 2) +
  geom_point(aes(col = label), size = 0.5) +
  scale_x_continuous(trans = 'reverse',
                     breaks = seq(10,-10,-1),
                     labels = 10^seq(10,-10,-1)) +
  scale_color_manual(name = 'Model with',
                     values = c('C. B + Hypothesized YLL058W reaction' = '#333333',
                                'D. C - All MET15 reactions' = '#999999'),
                     labels = c('Hypothesized *YLL058W* reaction<br/>w/ all *MET15* reactions',
                                'Hypothesized *YLL058W* reaction<br/>w/o all *MET15* reactions')) +
  labs(x = 'Yll058w/Met15 Kcat Ratio',
       y = 'Simulated Biomass Flux\n("Growth")') +
  # facet_zoom(ylim = c(0.3225,0.3275), zoom.size = .5) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow = 2, byrow=F, order = 1,
                             override.aes = list(shape = 15, size = 3, alpha = 1)))

ggsave(sprintf("%s/Simulation.jpg",fig_path), fig5e,
       height = one.c, width = one.c, units = 'mm',
       dpi = 600)

# fig4 <- cowplot::plot_grid(cowplot::plot_grid(fig4a,
#                                               cowplot::plot_grid(fig4c, fig4b,
#                                                                  nrow = 1, rel_widths = c(1,1),
#                                                                  labels = c('B','C'),
#                                                                  label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
#                                               nrow = 2, rel_heights = c(1,1),
#                                               labels = c('A',''),
#                                               label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
#                            cowplot::plot_grid(fig4d, fig4e,
#                                               labels = c('D','E'), ncol = 2,
#                                               rel_widths = c(1,1),
#                                               label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
#                                               align = 'hv', axis = 'tb'),
#                            ncol = 1, rel_heights = c(2,1))
# ggsave(sprintf("%s/Figure4.jpg",fig_path), fig4,
#        height = 250, width = two.c, units = 'mm',
#        dpi = 600)

##### EVO ANALYSIS
fig5a <- readPNG('figures/final/final/YLL_Locus.png')
fig5a <- ggplot() + 
  background_image(fig5a) +
  theme(plot.margin = margin(t=0, l=5, r=5, b=0, unit = "mm"),
        plot.background = element_blank())

# fig5b <- readPNG('figures/final/final/YLL_MDS.png')
# fig5b <- ggplot() + 
#   background_image(fig5b) +
#   theme(plot.margin = margin(t=10, l=0, r=0, b=10, unit = "mm"),
#         plot.background = element_blank())

fig5b <- readPNG('figures/final/final/tree2.png')
fig5b <- ggplot() + 
  background_image(fig5b) +
  theme(plot.margin = margin(t=0, l=5, r=5, b=0, unit = "mm"),
        plot.background = element_rect(fill = 'white'))

fig5c <- readPNG('figures/final/final/active site panel v1.png')
fig5c <- ggplot() + 
  background_image(fig5c) +
  theme(plot.margin = margin(t=5, l=0, r=0, b=5, unit = "mm"),
        plot.background = element_blank())
