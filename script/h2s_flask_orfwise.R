
fig6b.2 <- merge(data.pred.h2s %>%
                 melt(id.vars = 'Time', variable.name = 'Flask', value.name = 'OD'), strain.labs.h2s, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Time == 3075, Expected.Initial.OD == 0.1, Replicate != 'Rep_1')  %>%
  ggplot(aes(x = FeEDTA, y = OD)) +
  stat_summary(col = 'black', alpha = 0.9, size = 0.3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  stat_compare_means(method = 't.test', label = 'p.format',
                     label.x = 1.5, label.y = 8,
                     hjust = 0.5, size = 1.5) +
  scale_x_discrete(labels = c('FY4' = 'FY4',
                              'met15del' = 'FY4-*met15Δ*')) +
  labs(x = 'Biological Replicate',
       y = 'OD<sub>600</sub> at Saturation') +
  facet_grid(.~Strain, labeller = labeller(FeEDTA = c('Absent' = 'SD - Met + Glu<br />w/o H<sub>2</sub>S Chelator',
                                                      'Present' = 'SD - Met + Glu<br />w/ H<sub>2</sub>S Chelator'),
                                           Strain = c('FY4' = 'FY4',
                                                      'met15del' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        # axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1, vjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        # legend.direction = 'vertical',
        # legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt-2,
                                                margin = margin(1,0,0.1,0, "mm")),
        strip.text.y = ggtext::element_markdown(size = txt-2,
                                                margin = margin(0.5,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,8.1))


fig6c.2 <- data.h2s.re[data.h2s.re$Time..Minutes. == 3075 &
                       !(data.h2s.re$FeEDTA == 'Absent' & data.h2s.re$Replicate == 'R2') &
                       !(data.h2s.re$FeEDTA == 'Present' & data.h2s.re$Replicate == 'R3'),-2] %>%
  group_by(FeEDTA, ORF, Replicate) %>%
  summarize(OD = median(OD, na.rm = T), .groups = 'keep')  %>%
  ggplot(aes(x = FeEDTA, y = OD)) +
  stat_summary(col = 'black', alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  stat_compare_means(method = 't.test', label = 'p.format',
                     label.x = 1.5, label.y = 8,
                     hjust = 0.5, size = 1.5) +
  scale_x_discrete(labels = c('FY4' = 'FY4',
                              'met15D' = 'FY4-*met15Δ*')) +
  labs(x = 'Biological Replicate',
       y = 'OD<sub>600</sub> at Saturation') +
  facet_grid(.~ORF, labeller = labeller(FeEDTA = c('Absent' = 'SD - Met + Glu<br />(H<sub>2</sub>S Chelator<br />Previously Absent)',
                                                      'Present' = 'SD - Met + Glu<br />(H<sub>2</sub>S Chelator<br />Previously Present)'),
                                           ORF = c('FY4' = 'FY4',
                                                      'met15D' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        # axis.text.x = ggtext::element_markdown(size = txt, angle = 30, hjust = 1, vjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        # legend.direction = 'vertical',
        # legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt-2,
                                                margin = margin(1,0,0.1,0, "mm")),
        strip.text.y = ggtext::element_markdown(size = txt-2,
                                                margin = margin(0.1,0,0.1,0, "mm")))  +
  coord_cartesian(ylim = c(0,8.1))


# fig6 <- cowplot::plot_grid(fig6a, cowplot::plot_grid(fig6b, fig6c, nrow = 1, rel_widths = c(1,1),
#                                              labels = c('B','C'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
#                    ncol = 1, rel_heights = c(1,1),
#                    labels = c('A',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
# ggsave(sprintf("%s/Figure6.jpg",fig_path), fig6,
#        height = two.c, width = two.c, units = 'mm',
#        dpi = 600)


fig6a <- readPNG('figures/final/final/H2S_Flask_Schema2.png')
fig6a <- ggplot() + 
  background_image(fig6a) +
  theme(plot.margin = margin(t=0, l=40, r=40, b=0, unit = "mm"),
        plot.background = element_blank())

fig6d <- readPNG('figures/final/final/H2S_Flask_Spots2.png')
fig6d <- ggplot() + 
  background_image(fig6d) +
  theme(plot.margin = margin(t=0, l=5, r=5, b=0, unit = "mm"),
        plot.background = element_blank())

fig6.2 <- cowplot::plot_grid(fig6a, cowplot::plot_grid(cowplot::plot_grid(fig6b.2, fig6c.2, ncol = 1, rel_heights = c(1,1),
                                                                        labels = c('B','C'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                                                     fig6d, nrow = 1, rel_widths = c(1,2.6), 
                                                     labels = c('','D'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                           ncol = 1, rel_heights = c(1,1.5),
                           labels = c('A',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure6_.jpg",fig_path), fig6.2,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)
