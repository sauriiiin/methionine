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

fig.cbn.sc <- data.cbn %>%
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
  labs(x = 'Strain', y= 'Relative Colony Size',
       title = 'Synthetic Complete Media') +
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


fig.rp <- merge(data.rp[str_detect(data.rp$orf_name, 'BY'),], strain.labs.rp,
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


##### NO SULFATE
fig.ns <- data.ns[!(data.ns$id == 'Agar_Difco_+sulfate_+Met' & data.ns$stage == 'S1'),] %>%
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
                     labels = c('Agar_Difco_+sulfate_+Met' = 'Commercial YNB + Inorganic Sulfates + Methionine - Cysteine',
                                'Agarose_Home_+sulfate_-Met' = 'Homemade YNB + Inorganic Sulfates - Methionine - Cysteine',
                                'Agarose_Home_-sulfate_-Met' = 'Homemade YNB - Inorganic Sulfates - Methionine - Cysteine'),
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
        plot.margin = margin(0, 5, 0, 0, "mm"),
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


# ggsave(sprintf("%s/Figure3B.jpg",fig_path), fig3b,
#        height = 150, width = two.c*1.5, units = 'mm',
#        dpi = 600)


##### YLL JAKES MUTANT
fig.jm.yll <- merge(data.jm, strain.labs.jm, by = 'orf_name') %>%
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
                    values = c('Prototroph' = '#FFC107',
                               'Presumed Auxotroph' = '#536DFE',
                               'Unknown' = '#616161'),
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

##### PLASMID VALIDATION
fig.pv <- data.pv %>%
  filter(hours %in% c(48,115), stage %in% c('WC','PS1'),
         deletion1 %in% c('met15', 'met12'), deletion2 %in% c('','met12'),
         plasmid_backbone == '2M', plasmid_orf %in% c('empty','met12')) %>%
  ggplot(aes(x = orf_name, y = average)) +
  geom_boxplot(aes(alpha = stage), fill = 'grey40', outlier.shape = NA, size = 0.3,
               position = position_dodge(width = 0.6)) +
  # facet_grid(stage*plasmid_backbone ~ deletion1*deletion2*plasmid_orf,
  #            labeller = labeller(stage = c('WC' = 'YPDA', 'PS1' = 'SD-Met-Cys+Gal',
  #                                          'PS2' = 'SD-Met-Cys+Gal', 'FS' = 'SD-Met-Cys+Gal'))) +
  labs(y = 'Colony Size (pixels)') +
  scale_x_discrete(limits = c('met15del_2M_empty',
                              'met15del_2M_met12',
                              'met15del_met12del_2M_empty',
                              'met15del_met12del_2M_met12',
                              'met12del_2M_empty'),
                   labels = c('met15del_2M_empty' = 'FY4-*met15Δ*<br />+ Empty Plasmid',
                              'met15del_2M_met12' = 'FY4-*met15Δ*<br />+ MET12 Plasmid',
                              'met15del_met12del_2M_empty' = 'FY4-*met15Δ met12Δ*<br />+ Empty Plasmid',
                              'met15del_met12del_2M_met12'  = 'FY4-*met15Δ met12Δ*<br />+ MET12 Plasmid',
                              'met12del_2M_empty' = 'FY4-*met12Δ*<br />+ Empty Plasmid')) +
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
ggsave(sprintf("%s/FigureS_PV.jpg",fig_path), fig.pv,
       height = 100, width = two.c, units = 'mm',
       dpi = 600)
