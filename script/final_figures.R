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

##### FIGURES MADE OUTSIDE R
fig.cbn.a <- readPNG('figures/final/fig1a.png')
fig.cbn.a <- ggplot() + 
  background_image(fig.cbn.a) +
  theme(plot.margin = margin(t=0, l=20, r=20, b=0, unit = "mm"),
        plot.background = element_blank())

fig.cbn.b <- readPNG('figures/final/final/Figure1B.png')
fig.cbn.b <- ggplot() + 
  background_image(fig.cbn.b) +
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
        plot.background = element_blank())

# fig.yll.e <- readPNG('figures/branden/Fig6D.png')
# fig.yll.e <- ggplot() + 
#   background_image(fig.yll.e) +
#   theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
#         plot.background = element_blank())

##### CARBON SOURCE EXPERIMENT
data.cbn$methionine <- factor(data.cbn$methionine, levels = c('+Met +Ura', '+Met -Ura', '-Met +Ura'))
fig.cbn.c <- data.cbn %>%
  filter(methionine != '+Met +Ura', base == 'SD', orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  # geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(aes(fill = auxotrophy, alpha = methionine), size = 0.4, outlier.shape = NA,
               position = position_dodge(width = 0.6)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(labels = c( 'FY4' = 'FY4',
                               'FY4-met15del' = 'FY4-*met15Δ*',
                               'BY4742' = 'BY4742',
                               'BY4741' = 'BY4741')) +
  scale_fill_manual(name = 'Methionine-Uracil\nAuxotrophy',
                    values = c('None' = '#FFC107',
                               'Methionine' = '#536DFE',
                               'Uracil' = '#E040FB',
                               'Both' = '#FF5722'),
                    position = 'bottom') +
  scale_alpha_manual(name = 'Condition',
                     values = c('-Met +Ura' = 1,
                                '+Met -Ura' = 0.4)) +
  labs(x = 'Strain', y= 'Relative Colony Size') +
  facet_grid(.~carbon) +
  coord_cartesian(ylim = c(0,1.3)) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 40, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(fill = guide_legend(nrow=2, byrow=TRUE, order = 1),
         alpha = guide_legend(nrow=2, byrow=TRUE, order = 2,
                              override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(1,0.4)),
                                                colour=NA)))

fig.cbn <- cowplot::plot_grid(fig.cbn.a, fig.cbn.b, fig.cbn.c,
                             labels = c('A','B','C'), ncol = 1, rel_heights = c(1,1,1),
                             label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure_CBN.jpg",fig_path), fig.cbn,
       height = 250, width = two.c, units = 'mm',
       dpi = 600)

##### JAKES MUTANTS Pinning
fig.jm.a <- merge(data.jm, strain.labs.jm, by = 'orf_name') %>%
  filter(orf_name != 'yll', time == 't_final', attempt != 'pilot', condition %in% c('YPDA','SD-Met-Cys+Glu')) %>%
  ggplot(aes(x = orf_name, y = relative_cs,fill = auxotrophy, alpha = condition)) +
  geom_boxplot(outlier.shape = NA, size = 0.4,
               position = position_dodge(width = 0.6)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(labels = c('FY4'='FY4',
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
                              'BY4741'='BY4741'),
                   limits = c('FY4','met12','met3','met15','met2','met6', 
                              'met13','cys4','str3','BY4742','BY4741')) +
  scale_alpha_manual(name = 'Condition',
                        values = c('YPDA' = 0.4,
                                   'SD-Met-Cys+Glu' = 1)) +
  scale_fill_manual(name = 'Methionine-Uracil\nAuxotrophy',
                    values = c('Unknown' = '#FFC107',
                               'Presumed Auxotroph' = '#536DFE',
                               'Prototroph' = '#E040FB'),
                    position = 'bottom') +
  labs(x = 'Strain',
       y = 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 40, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(fill = guide_legend(nrow=2, byrow=TRUE, order = 1),
         alpha = guide_legend(nrow=2, byrow=TRUE, order = 2,
                              override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.4,1)),
                                                colour=NA)))

##### JAKES MUTANT Bismuth Biggy
fig.jm.b <- data.bis %>%
  # filter(Strain != 'FY4-cys4D') %>%
  ggplot(aes(x = Strain, y = HS)) +
  stat_summary(data = data.bis %>%
                 # filter(Strain != 'FY4-cys4D') %>% 
                 group_by(Condition, Strain) %>%
                 summarize(HS = mean(HS, na.rm = T), .groups = 'keep'),
               aes(fill = round(HS)), col = 'white', alpha = 0.9, size = 1,
               fun = mean, geom = "bar") +
  # stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_x_discrete(labels = strain.labs.bis$labels[c(1,2,4:9,3,10,11)],
                   limits = strain.labs.bis$orf_name[c(1,2,4:9,3,10,11)]) +
  scale_y_continuous(breaks = seq(0,10,1)) +
  scale_fill_gradient(low = "#D7CCC8", high = "#5D4037", guide = F) +
  labs(y = 'Relative Hydrogen Sulfide',
       x = 'Strains') +
  facet_wrap(.~Condition, ncol = 2) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(angle = 40, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(1,7))

fig.jm <- cowplot::plot_grid(fig.jm.a, fig.jm.b,
                             labels = c('A','B'), ncol = 1, rel_heights = c(1.3,1),
                             label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure_JM.jpg",fig_path), fig.jm,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)

##### RESPIRATION
fig.res.a <- merge(data.res.gc[str_detect(data.res.gc$orf_name, 'FY'),],
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
                    group_by(orf_name, condition, Time, base, methionine, carbon, cysteine, labels, met_aux, pet) %>%
                    summarize(rel_cs = mean(rel_cs, na.rm = T), .groups = 'keep'),
                  aes(x = Time, y = rel_cs, label = labels),
                  parse = T, size = 2.5) +
  facet_grid(carbon~methionine) +
  scale_color_manual(name = 'Methionine\nAuxotrophy',
                     values = c('Prototroph' = '#FFC107',
                                'Presumed Auxotroph' = '#536DFE'),
                     limits = c('Prototroph','Presumed Auxotroph'),
                     guide = F) +
  scale_fill_manual(name = 'Methionine\nAuxotrophy',
                    values = c('Prototroph' = '#FFC107',
                               'Presumed Auxotroph' = '#536DFE'),
                    limits = c('Prototroph','Presumed Auxotroph'),
                    position = 'bottom') +
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
  guides(fill = guide_legend(nrow=2, byrow=TRUE, order = 1),
         linetype = guide_legend(nrow=2, byrow=TRUE, order = 2))


data.res.gc.sum$methionine <- factor(data.res.gc.sum$methionine, levels = c('+Met','-Met'))
fig.res.b <- merge(data.res.gc.sum[str_detect(data.res.gc.sum$orf_name, 'FY'),],
      strain.labs.res, by = 'orf_name') %>%
  filter(expt_rep %in% c("1","2"), condition != 'SD+Met-Ura+Glu') %>% 
  ggplot(aes(x = orf_name, y = rel_auc)) +
  geom_boxplot(aes(fill = met_aux, alpha = methionine), size = 0.2, outlier.shape = NA,
               position = position_dodge(width = 0.6)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(labels = c( 'FY4' = 'FY4',
                               'FY4_pet' = 'FY4 (&rho;-)' ,
                               'FY4-met15del' = 'FY4-*met15Δ*', 
                               'FY4-met15del_pet' = 'FY4-*met15Δ* (&rho;-)', 
                               'BY4742' = 'BY4742', 
                               'BY4742_pet' = 'BY4742 (&rho;-)', 
                               'BY4741' = 'BY4741', 
                               'BY4741_pet' = 'BY4741 (&rho;-)')) +
  scale_fill_manual(name = 'Methionine\nAuxotrophy',
                    values = c('Prototroph' = '#FFC107',
                               'Presumed Auxotroph' = '#536DFE'),
                    limits = c('Prototroph','Presumed Auxotroph'),
                    guide = F) +
  scale_alpha_manual(name = 'Condition',
                     values = c('-Met' = 1,
                                '+Met' = 0.4)) +
  labs(x = 'Strain', y= 'Relative Area Under the Curve') +
  facet_grid(carbon~.) +
  coord_cartesian(ylim = c(0,1.3)) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 40, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(alpha = guide_legend(nrow=2, byrow=TRUE, order = 1,
                              override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(1,0.4)),
                                                colour=NA)))

fig.res.legend.a <- g_legend(fig.res.a)
fig.res.legend.b <- g_legend(fig.res.b)

fig.res <- cowplot::plot_grid(fig.res.a + theme(legend.position="none"),
                              fig.res.b + theme(legend.position="none"),
                              fig.res.legend.a, fig.res.legend.b,
                              labels = c('A','B','',''), ncol = 2, nrow = 2, 
                              rel_widths = c(2,1), rel_heights = c(10,1),
                              label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
                              align = 'hv', axis = 'tb')

ggsave(sprintf("%s/Figure_RES.jpg",fig_path), fig.res,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)


###### YLL RESULTS
fig.yll.a <- merge(data.del.diff %>% filter(orf_name %in% strain.labs.del$orf_name),
      strain.labs.del, by = 'orf_name') %>%
  filter(stage == 'Final Screen') %>%
  ggplot(aes(x = fitness_MM, y = fitness_PM)) +
  geom_point(aes(col = phenotype), size = 2) +
  geom_line(data = data.del.diff.dist %>% filter(stage == 'Final Screen'), aes(x = x , y = y.ul),
            linetype = 'dashed', size = 0.5) +
  geom_line(data = data.del.diff.dist %>% filter(stage == 'Final Screen'), aes(x = x , y = y.ll),
            linetype = 'dashed', size = 0.5) +
  geom_point(aes(x = fitness_MM, y = fitness_PM), shape = 1, size = 2) +
  geom_text_repel(aes(x = fitness_MM, y = fitness_PM, label = standard_name), size = 2,
                  force = 2, max.overlaps = 30) +
  scale_x_continuous(trans = 'pseudo_log') +
  labs(x = 'Relative Fitness in -Met -Cys',
       y = 'Relative Fitness in +Met -Cys') +
  scale_color_manual(name = 'Phenotype',
                     breaks = c('Beneficial','Neutral','Deleterious','Dead'),
                     values = c('Beneficial' = '#FF9800',
                                'Neutral'= '#E64A19',
                                'Deleterious' = '#795548',
                                'Dead' = 'black'),
                     na.translate = F) +
  coord_cartesian(xlim = c(0, 5),
                  ylim = c(0, 1.5)) +
  # facet_wrap(.~stage) +
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

fig.yll.b <- merge(data.jm, strain.labs.jm, by = 'orf_name') %>%
  filter(orf_name %in% c('FY4','met15','yll'), time == 't_final', attempt != 'pilot', condition %in% c('YPDA','SD-Met-Cys+Glu')) %>%
  ggplot(aes(x = orf_name, y = relative_cs,fill = auxotrophy, alpha = condition)) +
  geom_boxplot(outlier.shape = NA, size = 0.4,
               position = position_dodge(width = 0.6)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(labels = c('FY4'='FY4',
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
                              'BY4741'='BY4741'),
                   limits = c('FY4','met15','yll')) +
  scale_alpha_manual(name = 'Condition',
                     values = c('YPDA' = 0.4,
                                'SD-Met-Cys+Glu' = 1)) +
  scale_fill_manual(name = 'Methionine\nAuxotrophy',
                    values = c('Unknown' = '#FFC107',
                               'Presumed Auxotroph' = '#536DFE',
                               'Prototroph' = '#E040FB'),
                    position = 'bottom') +
  labs(x = 'Strain',
       y = 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        legend.box="vertical", legend.margin=margin(),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(fill = guide_legend(nrow=2, byrow=TRUE, order = 1),
         alpha = guide_legend(nrow=2, byrow=TRUE, order = 2,
                              override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.4,1)),
                                                colour=NA)))

fig.yll.c <- data.pv %>%
  filter(hours %in% c(48,115), stage %in% c('WC','PS1'),
         deletion1 %in% c('met15', 'yll'), deletion2 %in% c('','yll'),
         plasmid_backbone == '2M', plasmid_orf %in% c('empty','yll')) %>%
  ggplot(aes(x = orf_name, y = average)) +
  geom_boxplot(aes(alpha = stage), fill = 'grey40', outlier.shape = NA, size = 0.3,
               position = position_dodge(width = 0.6)) +
  # facet_grid(stage*plasmid_backbone ~ deletion1*deletion2*plasmid_orf,
  #            labeller = labeller(stage = c('WC' = 'YPDA', 'PS1' = 'SD-Met-Cys+Gal',
  #                                          'PS2' = 'SD-Met-Cys+Gal', 'FS' = 'SD-Met-Cys+Gal'))) +
  labs(y = 'Colony Size (pixels)') +
  scale_x_discrete(limits = c('met15del_2M_empty',
                              'met15del_2M_yll',
                              'met15del_ylldel_2M_empty',
                              'met15del_ylldel_2M_yll',
                              'ylldel_2M_empty'),
                   labels = c('met15del_2M_empty' = 'FY4-*met15Δ*<br />+ Empty Plasmid',
                              'met15del_2M_yll' = 'FY4-*met15Δ*<br />+ YLL058W Plasmid',
                              'met15del_ylldel_2M_empty' = 'FY4-*met15Δ yll058wΔ*<br />+ Empty Plasmid',
                              'met15del_ylldel_2M_yll'  = 'FY4-*met15Δ yll058wΔ*<br />+ YLL058W Plasmid',
                              'ylldel_2M_empty' = 'FY4-*yll058wΔ*<br />+ Empty Plasmid')) +
  # scale_y_continuous(trans = 'log10') +
  scale_alpha_manual(name = 'Condition',
                     values = c('WC' = 0.4,
                                'PS1' = 1),
                     labels = c('WC' = 'YPDA',
                                'PS1' = 'SD-Met-Cys+Gal')) +
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
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(alpha = guide_legend(nrow=1, byrow=TRUE, order = 1,
                              override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(0.4,1)),
                                                colour=NA)))

fig.yll.d <- data.bioc %>%
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


fig.yll.e <- data.mdl %>%
  ggplot(aes(x = log(Yll058w.Met15,10), y = Growth)) +
  geom_line(aes(col = label), lwd = 1) +
  geom_point(size = 2) +
  geom_point(aes(col = label), size = 0.5) +
  scale_x_continuous(trans = 'reverse',
                     breaks = seq(10,-10,-1),
                     labels = 10^seq(10,-10,-1)) +
  scale_color_manual(name = 'Model',
                     values = c('A. Default' = '#212121',
                                'B. A - All YLL058W reactions' = '#5D4037',
                                'C. B + Hypothesized YLL058W reaction' = '#FF5722',
                                'D. C - All MET15 reactions' = '#FFC107')) +
  labs(x = 'Yll058w/Met15 Kcat Ratio',
       y = 'Biomass Flux ("Growth")') +
  # facet_zoom(ylim = c(0.3225,0.3275), zoom.size = .5) +
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
  guides(color = guide_legend(nrow = 2, byrow=F, order = 1,
                             override.aes = list(shape = 15, size = 3, alpha = 1)))


fig.yll.legend.a <- g_legend(fig.yll.a)
fig.yll.legend.b <- g_legend(fig.yll.b)


fig.yll <- cowplot::plot_grid(cowplot::plot_grid(cowplot::plot_grid(fig.yll.a + theme(legend.position="none"),
                                                                    fig.yll.b + theme(legend.position="none"),
                                                                    fig.yll.legend.a, fig.yll.legend.b,
                                                                    labels = c('A','B','',''), ncol = 2, nrow = 2, 
                                                                    rel_widths = c(1,1), rel_heights = c(3,1),
                                                                    label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
                                                                    align = 'hv', axis = 'tb'),
                                                 fig.yll.c, labels = c('','C'),nrow = 2, 
                                                 rel_heights = c(1.4,1),
                                                 label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                              cowplot::plot_grid(fig.yll.d, fig.yll.e,
                                                 labels = c('D','E'), ncol = 2,
                                                 rel_widths = c(1,1),
                                                 label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
                                                 align = 'hv', axis = 'tb'),
                              ncol = 1, rel_heights = c(2,1))

ggsave(sprintf("%s/Figure_YLL.jpg",fig_path), fig.yll,
       height = 250, width = two.c, units = 'mm',
       dpi = 600)


##### NO SULFATE
strain.labs.ns <- strain.labs.cbn
strain.labs.ns$orf_name[strain.labs.ns$orf_name == 'FY4-met3del'] <- 'FY4_met3del'
strain.labs.ns$orf_name[strain.labs.ns$orf_name == 'FY4-met15del'] <- 'FY4_met15del'

fig3a <- merge(data.ns[data.ns$stage == 'S1' & data.ns$hours == 165 & data.ns$average != 0,],
      strain.labs.ns, by = 'orf_name') %>%
  ggplot(aes(x = orf_name, y = average, fill = met_aux, alpha = ynb_type)) +
  geom_boxplot(aes(linetype = base), size = 0.2, outlier.shape = NA,
               position = position_dodge(width = 0.6)) +
  labs(y = 'Colony Size (pixels)') +
  scale_x_discrete(labels = c( 'FY4' = 'FY4',
                               'FY4_met15del' = 'FY4-*met15Δ*',
                               'FY4_met3del' = 'FY4-*met3Δ*',
                               'BY4742' = 'BY4742',
                               'BY4741' = 'BY4741')) +
  scale_alpha_manual(name = 'YNB Type',
                     values = c('Difco' = 0.4,
                                'Home' = 1),
                     labels = c('Difco' = 'Commercial',
                                'Home' = 'Homemade')) +
  scale_fill_manual(name = 'Methionine\nAuxotrophy',
                    values = c('Presumed Auxotroph' = '#536DFE',
                               'Prototroph' = '#E040FB'),
                    labels = c('Presumed Auxotroph' = 'Presumed\nAuxotroph',
                               'Prototroph' = 'Prototroph'),
                    position = 'bottom') +
  scale_linetype_manual(name = 'Media Base',
                        values = c('Agar' = 'twodash',
                                  'Agarose' = 'solid')) +
  facet_wrap(methionine~sulfate, nrow = 3,
             labeller = labeller(methionine = c('+Met' = '+ Methionine',
                                                '-Met' = '- Methionine'),
                                 sulfate = c('+sulfate' = '+ Inorganic Sulfates',
                                             '-sulfate' = '- Inorganic Sulfates'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 40, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.box = 'vertical',
        # legend.direction = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        legend.spacing.y = unit(0.5, "mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(fill = guide_legend(nrow=2, byrow=TRUE, order = 1),
         alpha = guide_legend(nrow=2, byrow=TRUE, order = 3,
                              override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(1,0.4)),
                                                colour=NA)),
         linetype = guide_legend(nrow=2, byrow=TRUE, order = 2))
  

fig3b <- data.ns[!(data.ns$id == 'Agar_Difco_+sulfate_+Met' & data.ns$stage == 'S1'),] %>%
  filter(id %in% c('Agar_Difco_+sulfate_+Met',
                   'Agarose_Home_+sulfate_-Met',
                   'Agarose_Home_-sulfate_-Met')) %>%
  group_by(stage, base, ynb_type, sulfate, methionine, orf_name, hours, cum_hrs, id) %>%
  summarize(average = median(average, na.rm = T), .groups = 'keep') %>%
  data.frame() %>%
  ggplot(aes(x = cum_hrs, y = average)) +
  geom_line(aes(linetype = base, col = id)) +
  scale_x_continuous(breaks = seq(0,1000,50)) +
  scale_color_manual(name = 'Condition',
                     breaks = c('Agarose_Home_-sulfate_-Met',
                                'Agarose_Home_+sulfate_-Met',
                                'Agar_Difco_+sulfate_+Met'),
                    values = c('Agar_Difco_+sulfate_+Met' = '#212121',
                               'Agarose_Home_+sulfate_-Met' = '#795548',
                               'Agarose_Home_-sulfate_-Met' = '#FF5722'),
                    labels = c('Agar_Difco_+sulfate_+Met' = 'Commercial YNB + Inorganic Sulfates - Methionine',
                               'Agarose_Home_+sulfate_-Met' = 'Homemade YNB + Inorganic Sulfates - Methionine',
                               'Agarose_Home_-sulfate_-Met' = 'Homemade YNB - Inorganic Sulfates - Methionine'),
                    position = 'bottom') +
  scale_linetype_manual(name = 'Media Base',
                        values = c('Agar' = 'twodash',
                                   'Agarose' = 'solid')) +
  facet_wrap(.~orf_name, nrow = 5,
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
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        legend.spacing.y = unit(0.5, "mm"),
        strip.text = ggtext::element_markdown(size = txt, 
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=3, byrow=TRUE, order = 1, override.aes=list(size = 3)),
         linetype = guide_legend(nrow=2, byrow=TRUE, order = 2))


fig3a.legend <- g_legend(fig3a)
fig3b.legend <- g_legend(fig3b)

fig3 <- cowplot::plot_grid(fig3a + theme(legend.position="none"),
                   fig3b + theme(legend.position="none"),
                   fig3a.legend, fig3b.legend,
                   labels = c('A','B','',''), ncol = 2, nrow = 2, 
                   rel_widths = c(1,2), rel_heights = c(4,1),
                   label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
                   align = 'hv', axis = 'tb')

ggsave(sprintf("%s/Figure_NS.jpg",fig_path), fig3,
       height = 250, width = two.c, units = 'mm',
       dpi = 600)

