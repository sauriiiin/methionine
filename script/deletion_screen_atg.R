

goi <- c('YGL180W','YNL242W','YNR007C','YNL223W','YPL149W','YPL120W','YHR171W','YBL078C','YDL149W','YLL042C','YBR217W','YPR185W','YBR128C','YMR159C','YFR021W')
# c('YLR423C','YPL166W','YDR022C') - starvation-induced autophagy
# c('YPR049C','YOL082W','YDL113C','YPL100W','YLR431C','YJL036W','YJL178C') - CVT pathway
# c('YGL180W','YNL242W','YNR007C','YNL223W','YPL149W','YPL120W','YHR171W','YBL078C','YDL149W','YLL042C','YBR217W','YPR185W','YBR128C','YMR159C','YFR021W')  - core machinery of membrane formation

data.cont.lim

data.del.diff2 <- merge(data.del.diff, data.cont.lim, by = c('stage','arm_MM','arm_PM'))
data.del.diff2$norm_MM <- data.del.diff2$fitness_MM * data.del.diff2$cs_m_MM
data.del.diff2$norm_PM <- data.del.diff2$fitness_PM * data.del.diff2$cs_m_PM

data.del.diff.dist2 <- merge(data.del.diff.dist, data.cont.lim, by = c('stage'))
data.del.diff.dist2$x <- data.del.diff.dist2$x * data.del.diff.dist2$cs_m_MM
data.del.diff.dist2$y.ll <- data.del.diff.dist2$y.ll * data.del.diff.dist2$cs_m_PM
data.del.diff.dist2$y.ul <- data.del.diff.dist2$y.ul * data.del.diff.dist2$cs_m_PM

data.del.diff %>% filter(orf_name %in% goi)

merge(data.del.diff2 %>% filter(orf_name %in% strain.labs.del$orf_name),
      strain.labs.del, by = 'orf_name') %>%
  filter(orf_name %in% goi, stage == 'Final Screen') %>%
  ggplot(aes(x = norm_MM, y = norm_PM)) +
  geom_point(data = data.del.diff2 %>% filter(fitness_diff <= 10, stage == 'Final Screen'), size = 1, col = '#9E9E9E') +
  # geom_point(size = 2) +
  geom_line(data = data.del.diff.dist2 %>% filter(stage == 'Final Screen'),
            aes(x = x , y = y.ul),
            linetype = 'dashed', size = 0.5) +
  geom_line(data = data.del.diff.dist2 %>% filter(stage == 'Final Screen'),
            aes(x = x , y = y.ll),
            linetype = 'dashed', size = 0.5) +
  geom_point(aes(x = norm_MM, y = norm_PM), shape = 1, size = 2) +
  geom_point(data = data.del.diff2 %>% filter(orf_name %in%
                                                strain.labs.del$orf_name[strain.labs.del$standard_name %in% 
                                                                           c('YLL058W','MET12','MET5','MET10')],
                                              stage == 'Final Screen'),
             size = 1, col = '#FFC107') +
  geom_text_repel(aes(x = norm_MM, y = norm_PM, label = orf_name), size = 2,
                  min.segment.length = unit(0, 'lines'), seed = 10,
                  force = 2, max.overlaps = 30) +
  # scale_x_continuous(trans = 'pseudo_log') +
  labs(x = 'Relative Fitness in SD-Met+Gal (log)',
       y = 'Relative Fitness in SD+Met+Gal') +
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
  guides(color = guide_legend(nrow=3, byrow=TRUE, order = 1)) +
  coord_cartesian(xlim = c(0, 3200),
                  ylim = c(0, 1600))




##### DELETION SCREEN FIGURE
pmcs <- 1
data.cont.lim

merge(data.del.diff %>% filter(orf_name %in% strain.labs.del$orf_name, stage == 'Final Screen'),
      strain.labs.del, by = 'orf_name') 

fig.del1 <- merge(data.del.diff %>% filter(orf_name %in% strain.labs.del$orf_name, stage == 'Final Screen'),
      strain.labs.del, by = 'orf_name') %>%
  ggplot(aes(x = fitness_MM * 283/pmcs, y = fitness_PM * 1240/pmcs)) +
  geom_point(data = data.del.diff %>% filter(stage == 'Final Screen'), size = 1, col = '#9E9E9E') +
  # geom_point(size = 2) +
  geom_line(data = data.del.diff.dist %>% filter(stage == 'Final Screen'), 
            aes(x = x * 283/pmcs, y = y.ul * 1240/pmcs),
            linetype = 'dashed', size = 0.5) +
  geom_line(data = data.del.diff.dist %>% filter(stage == 'Final Screen'), 
            aes(x = x * 283/pmcs, y = y.ll * 1240/pmcs),
            linetype = 'dashed', size = 0.5) +
  geom_abline(linetype = 'dashed', size = 0.5, col = 'red') +
  geom_point(aes(x = fitness_MM * 283/pmcs, y = fitness_PM * 1240/pmcs), shape = 1, size = 2) +
  geom_point(data = data.del.diff %>% filter(orf_name %in%
                                               strain.labs.del$orf_name[strain.labs.del$standard_name %in% 
                                                                          c('YLL058W','MET12','MET5','MET10')],
                                             stage == 'Final Screen'), 
             size = 1, col = '#FFC107') +
  geom_text_repel(aes(x = fitness_MM * 283/pmcs, y = fitness_PM * 1240/pmcs, label = standard_name), size = 2,
                  min.segment.length = unit(0, 'lines'), seed = 10,
                  force = 2, max.overlaps = 30) +
  scale_x_continuous(breaks = seq(0,5000/pmcs,500/pmcs),
                     minor_breaks = seq(0,5000/pmcs,250/pmcs)) +
  scale_y_continuous(breaks = seq(0,5000/pmcs,500/pmcs),
                     minor_breaks = seq(0,5000/pmcs,250/pmcs)) +
  # scale_x_continuous(trans = 'pseudo_log') +
  labs(x = 'Normalized Colony Size\nin SD-Met+Gal',
       y = 'Normalized Colony Size\nin SD+Met+Gal') +
  coord_cartesian(xlim = c(0, 1800/pmcs),
                  ylim = c(0, 1800/pmcs)) +
  # facet_wrap(.~stage, ncol = 3,
  #            labeller = labeller(stage = c('Pre-Screen #1' = 'Pin #1',
  #                                          'Pre-Screen #2' = 'Pin #2',
  #                                          'Final Screen' = 'Pin #3'))) +
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


fig.del2 <- merge(data.del.diff %>% filter(orf_name %in% strain.labs.del$orf_name, stage == 'Final Screen'),
      strain.labs.del, by = 'orf_name') %>%
  ggplot(aes(x = fitness_MM * 283/pmcs, y = fitness_PM * 1240/pmcs)) +
  geom_point(data = data.del.diff %>% filter(stage == 'Final Screen'), size = 1, col = '#9E9E9E') +
  # geom_point(size = 2) +
  geom_line(data = data.del.diff.dist %>% filter(stage == 'Final Screen'), 
            aes(x = x * 283/pmcs, y = y.ul * 1240/pmcs),
            linetype = 'dashed', size = 0.5) +
  geom_line(data = data.del.diff.dist %>% filter(stage == 'Final Screen'), 
            aes(x = x * 283/pmcs, y = y.ll * 1240/pmcs),
            linetype = 'dashed', size = 0.5) +
  geom_point(aes(x = fitness_MM * 283/pmcs, y = fitness_PM * 1240/pmcs), shape = 1, size = 2) +
  geom_point(data = data.del.diff %>% filter(orf_name %in%
                                               strain.labs.del$orf_name[strain.labs.del$standard_name %in% 
                                                                          c('YLL058W','MET12','MET5','MET10')],
                                             stage == 'Final Screen'), 
             size = 1, col = '#FFC107') +
  scale_x_continuous(breaks = seq(0,5000/pmcs,2000/pmcs),
                     minor_breaks = seq(0,5000/pmcs,250/pmcs)) +
  labs(x = '\n',
       y = 'Normalized Colony Size\nin SD+Met+Gal') +
  coord_cartesian(xlim = c(1800/pmcs,4100/pmcs),
                  ylim = c(0, 1800/pmcs)) +
  # facet_wrap(.~stage, ncol = 3,
  #            labeller = labeller(stage = c('Pre-Screen #1' = 'Pin #1',
  #                                          'Pre-Screen #2' = 'Pin #2',
  #                                          'Final Screen' = 'Pin #3'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=3, byrow=TRUE, order = 1))


figDSP <- readPNG('figures/final/final/DeletionScreenPipeline.png')
figDSP <- ggplot() + 
  background_image(figDSP) +
  theme(plot.margin = margin(t=10, l=1, r=1, b=10, unit = "mm"),
        plot.background = element_blank()) 


fig6 <- plot_grid(figDSP, 
          plot_grid(fig.del1, fig.del2,
          nrow = 1, rel_widths = c(1,0.2)),
          ncol = 1, rel_heights = c(0.8,1),
          labels = c('A','B'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure6.jpg",fig_path), fig6,
       height = two.c, width = one.5c, units = 'mm',
       dpi = 600)



