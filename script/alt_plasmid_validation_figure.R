
unique(data.pv$deletion1)


data.pv2 %>%
  filter(orf_name %in% c('Plasmid_3'),
         arm == 'PV_FY_MM', stage != 'PS2') %>%
  group_by(stage) %>%
  summarize(cs = median(average, na.rm = T), .groups = 'keep')

data.pv2$yll_norm <- data.pv2$average
data.pv2$yll_norm[data.pv2$arm == 'PV_FY_MM' & data.pv2$stage == 'PS1'] <-
  data.pv2$yll_norm[data.pv2$arm == 'PV_FY_MM' & data.pv2$stage == 'PS1']/6638
data.pv2$yll_norm[data.pv2$arm == 'PV_FY_MM' & data.pv2$stage == 'PS1_2'] <-
  data.pv2$yll_norm[data.pv2$arm == 'PV_FY_MM' & data.pv2$stage == 'PS1_2']/2950

fig4c <- data.pv2 %>%
  filter(orf_name %in% c('Plasmid_3','Plasmid_5','Plasmid_7','Plasmid_11','Plasmid_13'),
         arm == 'PV_FY_MM', stage != 'PS2') %>%
  ggplot(aes(x = orf_name, y = yll_norm)) +
  geom_boxplot(outlier.shape = NA, fill = '#9E9E9E', size = 0.3) +
  scale_x_discrete(limits = c('Plasmid_3','Plasmid_5','Plasmid_11','Plasmid_7','Plasmid_13'),
                   labels = c('Plasmid_3' = 'FY4-*yll058wΔ*<br />w/ Empty Plasmid',
                              'Plasmid_5' = 'FY4-*met15Δ*<br />w/ Empty Plasmid',
                              'Plasmid_7' = 'FY4-*met15Δ yll058wΔ*<br />w/ Empty Plasmid',
                              'Plasmid_13' = 'FY4-*met15Δ yll058wΔ*<br />w/ *YLL058W* Plasmid',
                              'Plasmid_11' = 'FY4-*met15Δ*<br />w/ *YLL058W* Plasmid')) +
  labs(title = 'SD - Met + Gal + G418',
       y = 'Relative Colony Size') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        axis.text.x = ggtext::element_markdown(size = txt-1, angle = 30, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")))



unique(data.pv$orf_name)

data.pv %>%
  filter(hours %in% c(115), stage == 'PS1',
         orf_name %in% c('met12del_2M_empty')) %>%
  summarize(cs = median(average, na.rm = T))

data.pv$met12_norm <- data.pv$average
data.pv$met12_norm[data.pv$hours == 115 & data.pv$stage == 'PS1'] <-
  data.pv$met12_norm[data.pv$hours == 115 & data.pv$stage == 'PS1']/10623

fig4b <- data.pv %>%
  filter(hours %in% c(115), stage == 'PS1',
         orf_name %in% c('met12del_2M_empty',
                         'met15del_2M_empty',
                         'met15del_2M_met12',
                         'met15del_met12del_2M_empty',
                         'met15del_met12del_2M_met12')) %>%
  ggplot(aes(x = orf_name, y = met12_norm)) +
  geom_boxplot(fill = '#9E9E9E', outlier.shape = NA, size = 0.3) +
  labs(title = 'SD - Met + Gal + G418',
       y = 'Colony Size (pixels)') +
  scale_x_discrete(limits = c('met12del_2M_empty',
                              'met15del_2M_empty',
                              'met15del_2M_met12',
                              'met15del_met12del_2M_empty',
                              'met15del_met12del_2M_met12'),
                   labels = c('met12del_2M_empty' = 'FY4-*met12Δ*<br />w/ Empty Plasmid',
                              'met15del_2M_empty' = 'FY4-*met15Δ*<br />w/ Empty Plasmid',
                              'met15del_2M_met12' = 'FY4-*met15Δ*<br />w/ *MET12* Plasmid',
                              'met15del_met12del_2M_empty' = 'FY4-*met15Δ met12Δ*<br />w/ Empty Plasmid',
                              'met15del_met12del_2M_met12'  = 'FY4-*met15Δ met12Δ*<br />w/ *MET12* Plasmid',
                              'met12del_2M_empty' = 'FY4-*met12Δ*<br />w/ Empty Plasmid')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt-1, angle = 30, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))


figDSP <- readPNG('figures/final/final/DeletionScreenPipeline.png')
figDSP <- ggplot() + 
  background_image(figDSP) +
  theme(plot.margin = margin(t=0, l=25, r=25, b=0, unit = "mm"),
        plot.background = element_blank())


fig4 <- plot_grid(figDSP, fig4a + facet_wrap(.~stage, ncol = 3,
                                             labeller = labeller(stage = c('Pre-Screen #1' = 'Pin #1',
                                                                           'Pre-Screen #2' = 'Pin #2',
                                                                           'Final Screen' = 'Pin #3'))),
                  plot_grid(fig4b, fig4c, ncol = 2, rel_widths = c(1,1),
                            labels = c('C','D'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  ncol = 1, rel_heights = c(1,1,1),
                  labels = c('A','B',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure4_.jpg",fig_path), fig4,
       height = 220, width = two.c, units = 'mm',
       dpi = 600)
