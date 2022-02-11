

temp5 <- merge(data.h2sr %>%
                 melt(id.vars = c('Time', 'Attempt'), variable.name = 'Flask', value.name = 'OD'), 
               strain.labs.h2sr, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Time == 3075, Expected.Initial.OD == 0.1, Strain %in% c('BY4742','BY4741'),
         FeEDTA == 'Absent') %>%
  group_by(Attempt, Strain, Replicate) %>%
  summarise(OD = OD, .groups = 'keep') %>%
  data.frame()
temp5$Strain <- factor(temp5$Strain, levels = c('BY4742','BY4741'))

fig2BYOD <- temp5 %>%
  filter(Attempt == 'Original') %>%
  ggplot(aes(x = Strain, y = OD)) +
  stat_summary(col = 'black', alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  # scale_x_discrete(limits = c('BY4742','BY4741')) +
  stat_compare_means(method = 't.test', label = 'p.format',
                     label.x = 1.5, label.y = 7.9,
                     hjust = 0.5, size = 2.5) +
  labs(x = '',
       y = 'OD<sub>600</sub> at Saturation') +
  # facet_wrap(.~Attempt, ncol = 1) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'none',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                                margin = margin(1,0,0.1,0, "mm")))  +
  coord_cartesian(ylim = c(0,8.2))

fig.rp <- merge(data.rp, strain.labs.rp,
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
                  aes(x = cum_hrs, y = average - 400, label = labels), position = 'identity',
                  parse = T, size = 2.5) +
  geom_text(data = data.rp.rf,
            aes(x = cum_hrs, y = average + 300, 
                label = sprintf('%0.2f',relative_fitness)),
            size = 2.5) +
  scale_x_continuous(breaks = seq(-1000,1000,100)) +
  scale_color_manual(name = 'Presumed Auxotroph',
                     values = c('Methionine' = '#536DFE',
                                'None' = '#FFC107'),
                     labels = c('None' = 'No',
                                'Methionine' = 'Yes')) +
  scale_fill_manual(name = 'Presumed Auxotroph',
                    values = c('Methionine' = '#536DFE',
                               'None' = '#FFC107'),
                    labels = c('None' = 'No',
                               'Methionine' = 'Yes'),
                    guide = "none") +
  labs(x = 'Time (hours)', y = 'Colony Size (pixels)') +
  facet_grid(.~carbon, labeller = labeller(carbon = c('Glu'='SD-Met+Glu'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'none',
        # legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.1,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                                margin = margin(1,0,0.1,0, "mm")),
        strip.text.y = element_text(size = txt, margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(col = guide_legend(nrow=1, byrow=F, order = 1,
                            override.aes = list(size = 3))) +
  coord_cartesian(xlim = c(0,570))

fig.rp.expt <- merge(data.rp, strain.labs.rp,
      by = 'orf_name') %>%
  filter(aux == 'Met', orf_name %in% c('FY4','FY4-met15D'), carbon == 'Glu') %>%
  group_by(orf_name, auxotrophy, pin) %>%
  summarise(med_hrs = median(cum_hrs, na.rm = T),.groups = 'keep') %>%
  data.frame() %>%
  ggplot(aes(x = med_hrs, y = auxotrophy)) +
  geom_point(aes(col = auxotrophy), 
             shape = 15, size = 5) +
  geom_point(aes(x = med_hrs - 5, col = auxotrophy), 
             shape = 15, size = 5) +
  geom_point(aes(x = med_hrs + 5, col = auxotrophy), 
             shape = 15, size = 5) +
  geom_point(size = .5) +
  geom_text(data = merge(data.rp, strain.labs.rp,
                               by = 'orf_name') %>%
                    filter(aux == 'Met', orf_name %in% c('FY4',"FY4-met15D"), carbon == 'Glu', hours == max_hrs, pin == 7) %>%
                    group_by(orf_name, auxotrophy, labels, carbon, cum_hrs) %>%
                    summarize(average = mean(average, na.rm = T), .groups = 'keep'),
                  aes(x = cum_hrs, y = auxotrophy, label = labels), position = 'identity',
                  parse = T, size = 2.5) +
  scale_y_discrete(breaks = c('None',
                              'Methionine'),
                   limits = c('Methionine',
                              'None')) +
  scale_color_manual(name = 'Presumed Auxotroph',
                     values = c('Methionine' = '#536DFE',
                                'None' = '#FFC107'),
                     labels = c('None' = 'No',
                                'Methionine' = 'Yes'),
                     guide = 'none') + 
  theme_minimal() +
  theme(plot.title = element_blank(),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        panel.spacing = grid::unit(c(0,0,0,0), "mm"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0,570))


fig1b <- readPNG('figures/final/final/ExptPipe.png')
fig1b <- ggplot() + 
  background_image(fig1b) +
  theme(plot.margin = margin(t=1, l=17, r=17, b=0, unit = "mm"),
        plot.background = element_blank())

figBYue <- readPNG('figures/final/final/Fig2A.png')
figBYue <- ggplot() + 
  background_image(figBYue) +
  theme(plot.margin = margin(t=0, l=7, r=7, b=0, unit = "mm"),
        plot.background = element_blank())

figBYe <- readPNG('figures/final/final/Fig2D.png')
figBYe <- ggplot() + 
  background_image(figBYe) +
  theme(plot.margin = margin(t=0, l=8, r=8, b=0, unit = "mm"),
        plot.background = element_blank())

figBYschema <- readPNG('figures/final/final/H2SFlask/Fig 2B.png')
figBYschema <- ggplot() + 
  background_image(figBYschema) +
  theme(plot.margin = margin(t=0, l=12, r=12, b=0, unit = "mm"),
        plot.background = element_blank())

# figRPschema <- readPNG('figures/final/final/RepeatedPinningSchema.png')
# figRPschema <- ggplot() + 
#   background_image(figRPschema) +
#   theme(plot.margin = margin(t=1, l=10, r=15, b=0, unit = "mm"),
#         plot.background = element_blank())
# 
# figLegMUL <- readPNG('figures/final/final/MetUraLeu.png')
# figLegMUL <- ggplot() + 
#   background_image(figLegMUL) +
#   theme(plot.margin = margin(t=0, l=23, r=23, b=0, unit = "mm"),
#         plot.background = element_blank())

figLegYN <- readPNG('figures/final/final/YesNo.png')
figLegYN <- ggplot() + 
  background_image(figLegYN) +
  theme(plot.margin = margin(t=0, l=70, r=70, b=0, unit = "mm"),
        plot.background = element_blank())

fig2 <- plot_grid(fig1b,
                  plot_grid(figBYue, plot_grid(plot_grid(fig2a, fig2b.g, nrow = 1, rel_widths = c(1,1)),
                                               fig2b.leg, ncol = 1, rel_heights = c(10,0.8)),
                            nrow = 1, 
                            rel_widths = c(1,2),
                            labels = c('B','C'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  plot_grid(figBYschema, fig2BYOD, figBYe,
                            nrow = 1, 
                            rel_widths = c(1,1,1),
                            labels = c('D','E','F'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  plot_grid(fig.rp.expt, fig.rp,
                            ncol = 1, align = 'hv', axis = 'lr',
                            rel_heights = c(0.3,1)),
                  figLegYN,
                  ncol = 1, rel_heights = c(0.8,1,0.8,1.03,0.07),
                  labels = c('A','','','G',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure2.jpg",fig_path), fig2,
       height = 230, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 600)








