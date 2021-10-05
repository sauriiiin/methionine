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

##### CARBON SOURCE
fig.cbn <- data.cbn %>%
  filter(methionine == '+Met +Ura', base == 'SD', orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  # geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(aes(fill = auxotrophy), size = 0.2, outlier.shape = NA,
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
  # scale_alpha_manual(name = 'Condition',
  #                    values = c('-Met +Ura' = 1,
  #                               '+Met -Ura' = 0.4)) +
  labs(x = 'Strain', y= 'Relative Colony Size',
       title = 'Synthetic Defined Media') +
  facet_grid(.~methionine*carbon) +
  coord_cartesian(ylim = c(0,1.3)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        axis.title.x = element_blank(),
        axis.title = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(fill = guide_legend(nrow=2, byrow=TRUE, order = 1),
         alpha = guide_legend(nrow=2, byrow=TRUE, order = 2,
                              override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(1,0.4)),
                                                colour=NA)))

ggsave(sprintf("%s/FigureS_CBN.jpg",fig_path), fig.cbn,
       height = one.c, width = one.c, units = 'mm',
       dpi = 600)


fig.S.rp <- merge(data.rp[str_detect(data.rp$orf_name, 'BY'),], strain.labs.rp,
                  by = 'orf_name') %>%
  filter(hours == max_hrs, orf_name != 'BOR') %>%
  ggplot(aes(x = orf_name, y = fitness)) +
  geom_boxplot(aes(fill = auxotrophy, alpha = aux), size = 0.2, outlier.shape = NA,
               position = position_dodge(width = 0.6)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(limits = c('BY4742','BY4741'),
                   labels = c( 'FY4' = 'FY4',
                               'FY4-met15D' = 'FY4-*met15Δ*',
                               'BY4742' = 'BY4742',
                               'BY4741' = 'BY4741')) +
  scale_fill_manual(name = 'Methionine-Uracil\nAuxotrophy',
                    values = c('None' = '#FFC107',
                               'Methionine' = '#536DFE',
                               'Uracil' = '#E040FB',
                               'Both' = '#FF5722'),
                    position = 'bottom',
                    drop = F) +
  scale_alpha_manual(name = 'Condition',
                     values = c('Met' = 1,
                                'Ura' = 0.4),
                     labels = c('Met' = '-Met -Cys +Ura',
                                'Ura' = '+Met -Cys -Ura')) +
  # geom_boxplot(outlier.shape = NA) +
  labs(x = 'Strain', y = 'Relative Colony Size') +
  facet_grid(carbon~pin,
             labeller = labeller(pin = c('1' = 'Pin #1','2' = 'Pin #2',
                                         '3' = 'Pin #3','4' = 'Pin #4',
                                         '5' = 'Pin #5','6' = 'Pin #6',
                                         '7' = 'Pin #7'),
                                 carbon = c('Glu' = 'Glucose',
                                            'Gal' = 'Galactose'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(angle = 30, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(fill = guide_legend(nrow=2, byrow=TRUE, order = 1),
         alpha = guide_legend(nrow=2, byrow=TRUE, order = 2,
                              override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(1,0.4)),
                                                colour=NA)))

ggsave(sprintf("%s/FigureS_RP.jpg",fig_path), fig.S.rp,
       height = one.c, width = two.c, units = 'mm',
       dpi = 600)


##### RESPIRATION
data.res.gc.sum$methionine <- factor(data.res.gc.sum$methionine, levels = c('+Met','-Met'))

fig.res.a <- merge(data.res.gc[str_detect(data.res.gc$orf_name, 'BY'),],
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

fig.res.b <- merge(data.res.gc.sum, strain.labs.res, by = 'orf_name') %>%
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
  guides(alpha = guide_legend(nrow=2, byrow=TRUE, order = 1,
                              override.aes=list(fill=hcl(c(15,195),100,0,alpha=c(1,0.4)),
                                                colour=NA)))

fig.res.a.legend <- g_legend(fig.res.a)
fig.res.b.legend <- g_legend(fig.res.b)

fig.res <- cowplot::plot_grid(fig.res.a + theme(legend.position="none"),
                              fig.res.b + theme(legend.position="none"),
                              cowplot::plot_grid(fig.res.a.legend, fig.res.b.legend,
                                                 rel_widths = c(1.5,1)),
                              labels = c('A','B',''), ncol = 1, nrow = 3, 
                              rel_heights = c(1.5,1,0.3),
                              label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
                              align = 'hv', axis = 'tb')
ggsave(sprintf("%s/FigureS_RES.jpg",fig_path), fig.res,
       height = 250, width = two.c, units = 'mm',
       dpi = 600)
