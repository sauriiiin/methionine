data.bioc %>%
  group_by(`Time (min)`, Sample) %>%
  summarize(uM = mean(uM, na.rm = T), .groups = 'keep') %>%
  filter(`Time (min)` == 15)

fig5d <- data.bioc %>%
  # filter(outlier == FALSE) %>%
  ggplot(aes(x = `Time (min)`, y = uM)) +
  stat_summary(aes(group = Sample, fill = Sample), fun.data=mean_se, fun.args = list(mult=1), geom="ribbon",
               alpha = 0.4) +
  stat_summary(aes(group = Sample, col = Sample), fun=mean, geom="line", lwd = 0.7) +
  geom_segment(aes(x = 15.3, y = 8, xend = 15.3, yend = 119), size = 0.5) +
  annotate("text", x = 15.6, y = 30, 
           label = "p == 4.46~X~10^-5",
           angle = 270, size = 2,
           parse = TRUE) +
  geom_segment(aes(x = 16, y = 8, xend = 16, yend = 12), size = 0.5) +
  annotate("text", x = 16.3, y = 9.8, 
           label = "p == 0.00018",
           angle = 270, size = 2,
           parse = TRUE) +
  stat_summary(aes(group = Sample), fun=mean, geom="point", size =2) +
  stat_summary(aes(group = Sample, col = Sample), fun=mean, geom="point", size = 0.5) +
  labs(x = 'Time (minutes)',
       y = 'Homocysteine (μM)') +
  scale_color_manual(name = 'Sample',
                     values = c('Met15' = "#330099", #6600CC
                                'Yll058w' = "#CC33FF",
                                'None' = "#CC99CC"),
                     guide = 'none') +
  scale_fill_manual(name = 'Sample',
                    values = c('Met15' = "#330099",
                               'Yll058w' = "#CC33FF",
                               'None' = "#CC99CC")) +
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
  coord_cartesian(xlim = c(0,16),
                  ylim = c(5,125)) +
  guides(fill = guide_legend(override.aes=list(shape = 15, alpha = 1)))

fig5a <- readPNG('figures/final/final/YLL_Locus.png')
fig5a <- ggplot() + 
  background_image(fig5a) +
  theme(plot.margin = margin(t=5, l=1, r=1, b=1, unit = "mm"),
        plot.background = element_blank())

# fig5b <- readPNG('figures/final/final/tree3.png')
# fig5b <- ggplot() + 
#   background_image(fig5b) +
#   theme(plot.margin = margin(t=0, l=5, r=5, b=0, unit = "mm"),
#         plot.background = element_rect(fill = 'white'))

fig5c <- readPNG('figures/final/final/active site panel v1.png')
fig5c <- ggplot() + 
  background_image(fig5c) +
  theme(plot.margin = margin(t=0, l=20, r=20, b=0, unit = "mm"),
        plot.background = element_blank())

figPV2PS1_2 <- readPNG('figures/final/final/ColonyPics_PV2.png')
figPV2PS1_2 <- ggplot() + 
  background_image(figPV2PS1_2) +
  theme(plot.margin = margin(t=0, l=16, r=3, b=0, unit = "mm"),
        plot.background = element_blank()) 


fig4 <- plot_grid(plot_grid(plot_grid(fig5a, fig5c, 
                                      ncol = 1, rel_heights = c(0.8,1),
                                      labels = c('A','C'),
                                      label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                            fig5b, 
                            nrow = 1, rel_widths = c(1.5,1),
                            labels = c('','B'),
                            label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  figPV2PS1_2,
                  fig4c,
                  fig5d,
                  ncol = 1, rel_heights = c(2.2,0.6,0.9,1),
                  labels = c('','D','','E'),
                  label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')

ggsave(sprintf("%s/Figure4_.jpg",fig_path), fig4,
       height = 230, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 600)



