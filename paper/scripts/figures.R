
# source('/home/sbp29/R/Projects/methionine/paper/scripts/data.R')
source('/home/sbp29/R/Projects/methionine/paper/scripts/initialize.R')
load('/home/sbp29/R/Projects/methionine/paper/data/FigureData.RData')
outsidepanels <- 'paper/figures/PanelsMadeOutsideR/'

##### PANEL BACKGROUND COLOR
pinningpanelbg <- '#FFF6DB'
flaskpanelbg <- '#E1E4E8' ##D6E4F4

##### FIGURE 1: The unexpected growth of met15Δ cells on organosulfur-free media is both stable and dependent on cell propagation method.
fig1a <- readPNG(sprintf('%sFig1A.png',outsidepanels))
fig1a <- ggplot() + 
  background_image(fig1a) +
  theme(plot.margin = margin(t=1, l=10, r=10, b=0, unit = "mm"),
        plot.background = element_blank())

fig1b <- readPNG(sprintf('%sFig1B.png',outsidepanels))
fig1b <- ggplot() + 
  background_image(fig1b) +
  theme(plot.margin = margin(t=0, l=5, r=5, b=0, unit = "mm"),
        plot.background = element_blank())

fig1c.glu <- rbind(data.cbn[,c(1,2,5,6,10)],
               data.cbn.leu[,c(8,11,10,13,14)] %>%
                 filter(condition == 'SDmLeu')) %>%
  filter(carbon == 'Glucose', base == 'SD', orf_name %in% c('BY4741','BY4742')) %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  geom_boxplot(fill = '#9E9E9E', size = 0.2, outlier.shape = NA) +
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
  labs(x = '', y= 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.6)) +
  facet_grid(orf_name ~ carbon, scales = 'free_x') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        plot.margin = margin(2,5,0,5),
        panel.background = element_rect(fill = pinningpanelbg),
        # axis.title.x = element_blank(),
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

fig1c.galeth <- data.cbn %>%
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
  labs(x = '', y= 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.6)) +
  facet_grid(orf_name ~ carbon, scales = 'free_x') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        plot.margin = margin(2,5,0,5),
        panel.background = element_rect(fill = pinningpanelbg),
        # axis.title.x = element_blank(),
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

fig1c.galeth <- ggplot_gtable(ggplot_build(fig1c.galeth))
stripr <- which(grepl('strip-r', fig1c.galeth$layout$name))

fig1c.galeth$grobs[[16]]$grobs[[1]]$children[[2]]$children[[1]]$label
fig1c.galeth$grobs[[16]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
fig1c.galeth$grobs[[16]]$grobs[[1]]$children[[1]]$gp$fill <- '#FF0000'#'#E040FB'

fig1c.galeth$grobs[[17]]$grobs[[1]]$children[[2]]$children[[1]]$label
fig1c.galeth$grobs[[17]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
fig1c.galeth$grobs[[17]]$grobs[[1]]$children[[1]]$gp$fill <- '#990000'#'#FF5722'

fig1c.leg <- g_legend(data.frame(labels = c('Uracil-Leucine','Methionine-Uracil-Leucine')) %>%
                        ggplot(aes(x = 1, y = labels, col = labels)) +
                        geom_point(shape = 15) +
                        scale_color_manual(name = 'Presumed Auxotrophy',
                                           values = c('Uracil-Leucine' = '#FF0000',
                                                      'Methionine-Uracil-Leucine' = '#990000')) +
                        theme_linedraw() +
                        theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
                              axis.title = element_text(size = titles),
                              axis.text = element_text(size = txt),
                              legend.margin = margin(0,0,0,0),
                              legend.title =  ggtext::element_markdown(size = titles),
                              legend.text =  ggtext::element_markdown(size = txt),
                              legend.spacing = unit(0.1,"mm"),
                              legend.key.size = unit(4, "mm"),
                              legend.direction = 'horizontal',
                              legend.position = 'bottom') +
                        guides(col = guide_legend(nrow=1, override.aes=list(size = 3))))

fig1d <- readPNG(sprintf('%sFig1D.png',outsidepanels))
fig1d <- ggplot() + 
  background_image(fig1d) +
  theme(plot.margin = margin(t=0, l=12, r=12, b=0, unit = "mm"),
        plot.background = element_blank())

fig1e <- temp5 %>%
  filter(Attempt == 'Original') %>%
  ggplot(aes(x = Strain, y = OD)) +
  stat_summary(col = 'black', fill = '#D6E4F4',
               alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  # scale_x_discrete(limits = c('BY4742','BY4741')) +
  stat_compare_means(method = 't.test', label = 'p.format',
                     label.x = 1.5, label.y = 6.5,
                     hjust = 0.5, size = 2.5) +
  labs(x = '',
       y = 'OD<sub>600</sub> at Saturation',
       title = 'SD-Met-Cys+Glu') +
  # facet_wrap(.~Attempt, ncol = 1) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5),
        panel.background = element_rect(fill = flaskpanelbg),
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
  coord_cartesian(ylim = c(0,7))

fig1f <- readPNG(sprintf('%sFig1F.png',outsidepanels))
fig1f <- ggplot() + 
  background_image(fig1f) +
  theme(plot.margin = margin(t=0, l=8, r=8, b=0, unit = "mm"),
        plot.background = element_blank())

fig1g <-merge(data.rp, strain.labs.rp,
              by = 'orf_name') %>%
  filter(aux == 'Met', orf_name %in% c('FY4'), carbon == 'Glu') %>%
  ggplot(aes(x = cum_hrs, y = average)) +
  stat_summary(aes(fill = auxotrophy, group = pin),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.4) +
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
  facet_grid(.~carbon, labeller = labeller(carbon = c('Glu'='SD-Met-Cys+Glu'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        panel.background = element_rect(fill = pinningpanelbg),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.1,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                                margin = margin(1,0,0.1,0, "mm")),
        strip.text.y = element_text(size = txt, margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(col = guide_legend(nrow=1, byrow=F, order = 1,
                            override.aes = list(size = 3))) +
  coord_cartesian(xlim = c(0,570))

fig1g.expt <- merge(data.rp, strain.labs.rp,
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
            aes(x = cum_hrs-5, y = auxotrophy, label = labels), position = 'identity',
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

fig1g.leg <- g_legend(data.frame(labels = c('Yes','No')) %>%
                        ggplot(aes(x = 1, y = labels, col = labels)) +
                        geom_point(shape = 15) +
                        scale_color_manual(name = 'Presumed Auxotroph',
                                           values = c('Yes' = '#536DFE',
                                                      'No' = '#FFC107'),
                                           breaks = c('Yes','No')) +
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

fig1 <- plot_grid(fig1a,
                  plot_grid(NULL, plot_grid(plot_grid(fig1c.glu, fig1c.galeth, nrow = 1, rel_widths = c(1,1)),
                                               fig1c.leg, ncol = 1, rel_heights = c(10,0.7)),
                            nrow = 1, 
                            rel_widths = c(1,2),
                            labels = c('B','C'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  plot_grid(NULL, fig1e, NULL,
                            nrow = 1, 
                            rel_widths = c(1,1,1),
                            labels = c('D','E','F'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  plot_grid(fig1g.expt, fig1g,
                            ncol = 1, align = 'hv', axis = 'lr',
                            rel_heights = c(0.3,1)),
                  ncol = 1, rel_heights = c(0.8,1,0.8,1.03),
                  labels = c('A','','','G'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure1.jpg",fig_path), fig1,
       height = 230, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)

##### FIGURE 2: Growth of met15Δ colonies in organosulfur-free media is specific and dependent on the utilization of inorganic sulfates.
# fig2a <- readPNG(sprintf('%sFig2A.png',outsidepanels))
# fig2a <- ggplot() + 
#   background_image(fig2a) +
#   theme(plot.margin = margin(t=0, l=35, r=35, b=0, unit = "mm"),
#         plot.background = element_blank())

fig2a <- readPNG(sprintf('%sFig2A.png',outsidepanels))
fig2a <- ggplot() + 
  background_image(fig2a) +
  theme(plot.margin = margin(t=1, l=25, r=25, b=0, unit = "mm"),
        plot.background = element_blank())

fig2b <- data.jm.2[!(data.jm.2$orf_name %in% c('yll','BY4742','BY4741')),] %>%
  ggplot() +
  geom_boxplot(aes(x = condition, y = relative_fitness),
               fill = '#9E9E9E', outlier.shape = NA, size = 0.2,
               position = position_dodge(width = 0.85)) +
  scale_y_continuous(breaks = seq(-1,2,0.4)) +
  scale_x_discrete(limits = c('YPDA', 'SD-Met'),
                   labels = c('YPDA' = 'YPDA',
                              'SD-Met' = 'SD-Met-Cys+Glu')) +
  scale_fill_gradient(name = 'Relative Hydrogen Sulfide',
                      low = "#D7CCC8", high = "#5D4037",
                      limits = c(1,7)) +
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
                                              'BY4741'='BY4741',
                                              'met5del'='FY4-*met5Δ*',
                                              'met10del'='FY4-*met10Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        panel.background = element_rect(fill = pinningpanelbg),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = titles),
        axis.text.x = element_text(size = txt, angle = 15, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        # legend.direction = 'vertical',
        # legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'black',face = 'bold',
                                              margin = margin(0.5,0,0.1,0, "mm")))

fig2b <- ggplot_gtable(ggplot_build(fig2b))
stripr <- which(grepl('strip-t', fig2b$layout$name))

for (s in stripr) {
  # fig3b.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
  l <- fig2b$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
  
  if (length(l) != 0) {
    if (l == 'FY4-') {
      l2 <- fig2b$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$label
      a <- strain.labs.jm$auxotrophy[str_detect(strain.labs.jm$labels, l2)]
    } else {
      a <- strain.labs.jm$auxotrophy[strain.labs.jm$orf_name == l]
    }
    
    if (a == 'Presumed Auxotroph') {
      fig2b$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
    } else if (a == 'Prototroph') {
      fig2b$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
      # fig2b$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$gp$col <- 'black'
    } else {
      fig2b$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107' #Unkown = '#BDBDBD'
    }
  }
}

fig2c <- readPNG(sprintf('%sFig2C.png',outsidepanels))
fig2c <- ggplot() + 
  background_image(fig2c) +
  theme(plot.margin = margin(t=0, l=10, r=8, b=0, unit = "mm"),
        plot.background = element_blank())

fig2d <- data.ns2 %>%
  ggplot(aes(x = cum_hrs, y = average)) +
  # stat_summary(data = data.ns2 %>% filter(stage == 'S1'),
  #              aes(fill = id),
  #              fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S1'),
               aes(col = id),
               fun=mean, geom="line", lwd = 1) +
  # stat_summary(data = data.ns2 %>% filter(stage == 'S2'),
  #              aes(fill = id),
  #              fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S2'),
               aes(col = id),
               fun=mean, geom="line", lwd = 1) +
  # stat_summary(data = data.ns2 %>% filter(stage == 'S3'),
  #              aes(fill = id),
  #              fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S3'),
               aes(col = id),
               fun=mean, geom="line", lwd = 1) +
  # stat_summary(data = data.ns2 %>% filter(stage == 'S4'),
  #              aes(fill = id),
  #              fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S4'),
               aes(col = id),
               fun=mean, geom="line", lwd = 1) +
  # stat_summary(data = data.ns2 %>% filter(stage == 'S5'),
  #              aes(fill = id),
  #              fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S5'),
               aes(col = id),
               fun=mean, geom="line", lwd = 1) +
  # stat_summary(data = data.ns2 %>% filter(stage == 'S6'),
  #              aes(fill = id),
  #              fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S6'),
               aes(col = id),
               fun=mean, geom="line", lwd = 1) +
  scale_fill_manual(name = 'Condition',
                    breaks = c('Agarose_Home_-sulfate_-Met',
                               'Agarose_Home_+sulfate_-Met',
                               'Agar_Difco_+sulfate_+Met'),
                    values = c('Agar_Difco_+sulfate_+Met' = '#9E9E9E',
                               'Agarose_Home_+sulfate_-Met' = '#212121',
                               'Agarose_Home_-sulfate_-Met' = '#FF5722'),
                    labels = c('Agar_Difco_+sulfate_+Met' = 'SD+Met-Cys+Glu w/\nInorganic Sulfates',
                               'Agarose_Home_+sulfate_-Met' = 'SD-Met-Cys+Glu w/\nInorganic Sulfates',
                               'Agarose_Home_-sulfate_-Met' = 'SD-Met-Cys+Glu w/o\nInorganic Sulfates'),
                    guide = F) +
  scale_color_manual(name = 'Condition',
                     breaks = c('Agarose_Home_-sulfate_-Met',
                                'Agarose_Home_+sulfate_-Met',
                                'Agar_Difco_+sulfate_+Met'),
                     values = c('Agar_Difco_+sulfate_+Met' = '#9E9E9E',
                                'Agarose_Home_+sulfate_-Met' = '#212121',
                                'Agarose_Home_-sulfate_-Met' = '#FF5722'),
                     labels = c('Agar_Difco_+sulfate_+Met' = 'SD+Met-Cys+Glu w/ Inorganic Sulfates',
                                'Agarose_Home_+sulfate_-Met' = 'SD-Met-Cys+Glu w/ Inorganic Sulfates',
                                'Agarose_Home_-sulfate_-Met' = 'SD-Met-Cys+Glu w/o Inorganic Sulfates')) +
  # scale_x_discrete(breaks = c('FY4','FY4_met15del','FY4_met3del'),
  #                  labels = c('FY4' = 'FY4',
  #                             'FY4_met15del' = 'FY4-*met15Δ*',
  #                             'FY4_met3del' = 'FY4-*met3Δ*')) +
  scale_x_continuous(breaks = seq(-100,1000,100)) +
  scale_y_continuous(breaks = seq(-4500,13500,4500)) +
  facet_wrap(.~orf_name, ncol = 1,
             labeller = labeller(stage = c('S1' = 'Pin #1',
                                           'S2' = 'Pin #2',
                                           'S3' = 'Pin #3',
                                           'S4' = 'Pin #4',
                                           'S5' = 'Pin #5',
                                           'S6' = 'Pin #6'),
                                 orf_name = c('FY4' = 'FY4',
                                              'FY4_met15del' = 'FY4-*met15Δ*',
                                              'FY4_met3del' = 'FY4-*met3Δ*'))) +
  labs(x = 'Time (hours)',
       y = 'Colony Size (pixels)') +
  theme_bw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        panel.background = element_rect(fill = pinningpanelbg, color = 'transparent'),
        # panel.grid = element_line(color = '#FFFFFF'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.margin = margin(2,20,2,-5),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.2,"mm"),
        legend.spacing.y = unit(1, "mm"),
        panel.spacing.x = unit(0, "lines"),
        strip.text = ggtext::element_markdown(size = txt, 
                                              margin = margin(1,0,0,0, "mm")),
        strip.text.y = ggtext::element_markdown(size = txt,
                                                margin = margin(0,1,0,0, "mm"))) +
  guides(color = guide_legend(nrow=1, byrow=FALSE, order = 1, 
                              override.aes = list(size = 3))) +
  coord_cartesian(xlim = c(0,910))

fig2d.expt <- merge(data.ns2 %>%
  filter(id %in% c('Agarose_Home_-sulfate_-Met',
                   'Agarose_Home_+sulfate_-Met',
                   'Agar_Difco_+sulfate_+Met')) %>%
  group_by(stage, id) %>%
  summarize(stage_hrs = median(cum_hrs, na.rm = T), .groups = 'keep') %>%
  data.frame(),
  data.ns2 %>%
    filter(id %in% c('Agarose_Home_-sulfate_-Met',
                     'Agarose_Home_+sulfate_-Met',
                     'Agar_Difco_+sulfate_+Met')) %>%
    group_by(stage) %>%
    summarize(med_hrs = median(cum_hrs, na.rm = T), .groups = 'keep') %>%
    data.frame(), by = 'stage') %>%
  ggplot(aes(x = med_hrs, y = id)) +
  geom_point(aes(col = id), 
             shape = 15, size = 5) +
  geom_point(aes(x = med_hrs - 8, col = id), 
             shape = 15, size = 5) +
  geom_point(aes(x = med_hrs + 8, col = id), 
             shape = 15, size = 5) +
  geom_point(size = .5) +
  scale_y_discrete(breaks = c('Agarose_Home_-sulfate_-Met',
                              'Agarose_Home_+sulfate_-Met',
                              'Agar_Difco_+sulfate_+Met'),
                   limits = c('Agar_Difco_+sulfate_+Met',
                              'Agarose_Home_-sulfate_-Met',
                              'Agarose_Home_+sulfate_-Met')) +
  scale_color_manual(name = 'Condition',
                     breaks = c('Agarose_Home_-sulfate_-Met',
                                'Agarose_Home_+sulfate_-Met',
                                'Agar_Difco_+sulfate_+Met'),
                     values = c('Agar_Difco_+sulfate_+Met' = '#9E9E9E',
                                'Agarose_Home_+sulfate_-Met' = '#212121',
                                'Agarose_Home_-sulfate_-Met' = '#FF5722'),
                     labels = c('Agar_Difco_+sulfate_+Met' = 'SD+Met+Glu w/\nInorganic Sulfates',
                                'Agarose_Home_+sulfate_-Met' = 'SD-Met+Glu w/\nInorganic Sulfates',
                                'Agarose_Home_-sulfate_-Met' = 'SD-Met+Glu w/o\nInorganic Sulfates'),
                     guide = 'none') +
  theme_minimal() +
  theme(plot.title = element_blank(),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        panel.spacing = grid::unit(c(0,0,0,0), "mm"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0,910))

fig2 <- plot_grid(fig2a, fig2b, NULL,
                  plot_grid(fig2d.expt, fig2d,
                            ncol = 1, align = 'hv', axis = 'lr',
                            rel_heights = c(0.275,1)) ,
                  ncol = 1, rel_heights = c(1.2,1,0.8,1.6),
                  labels = c('A','B','C','D'),
                  label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure2.jpg",fig_path), fig2,
       height = 230, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


##### FIGURE 3: Yll058w functions as an inefficient homocysteine synthase.
fig3a <- readPNG(sprintf('%sFig3A.png',outsidepanels))
fig3a <- ggplot() + 
  background_image(fig3a) +
  theme(plot.margin = margin(t=5, l=1, r=1, b=1, unit = "mm"),
        plot.background = element_blank())

fig3b.tree <- ggtree(tree1)+
  geom_tiplab(size = 2.8)+
  xlim(-8, 36-2)+
  annotate("rect",xmin=29.4-2.2,xmax=30.6-1.8,ymin=(1:30)-.5,ymax=(1:30)+.3,fill=rev(clus1colors2),col='transparent') +
  annotate("rect",xmin=31.4-2.2,xmax=32.6-1.8,ymin=(1:30)-.5,ymax=(1:30)+.3,fill=rev(clus2colors2),col='transparent') +
  annotate("rect",xmin=33.4-2.2,xmax=34.6-1.8,ymin=(1:30)-.5,ymax=(1:30)+.3,fill=rev(clus3colors2),col='transparent') +
  annotate("point",x=30-2,y=(1:30)-.1,col=rev(clus1colors),size=1) +
  annotate("point",x=32-2,y=(1:30)-.1,col=rev(clus2colors),size=1) +
  annotate("point",x=34-2,y=(1:30)-.1,col=rev(clus3colors),size=1)

fig3b.leg <- g_legend(data.frame(levels = c('a','b','c','d','e')) %>%
  ggplot(aes(x = 1, y = levels, col = levels, shape = levels)) +
  geom_point() +
  scale_color_manual(name = 'Classes of YLL058W<br/>homologs',
                     values = c('a' = '#99CCFF',
                                'b' = '#CC33FF',
                                'c' = '#006666',
                                'd' = '#212121',
                                'e' = '#9E9E9E'),
                     labels = c('a' = 'Ancestral class',
                                'b' = 'YLL058W class',
                                'c' = '*STR2*/YML082W class',
                                'd' = 'Yes',
                                'e' = 'No'),
                     limits = c('a','b','c')) +
  scale_shape_manual(name = 'Present in sulfur-<br/>related cluster',
                     values = c('a' = 16,
                                'b' = 16,
                                'c' = 16,
                                'd' = 15,
                                'e' = 15),
                     labels = c('a' = 'Ancestral class',
                                'b' = '*YLL058W* class',
                                'c' = '*STR2/YML082W* class',
                                'd' = 'Yes',
                                'e' = 'No'),
                     limits = c('d','e')) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title =  ggtext::element_markdown(size = titles),
        legend.text =  ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.spacing = unit(0.1,"mm"),
        legend.key.size = unit(4, "mm")) +
  guides(col = guide_legend(order = 1, nrow=3, override.aes=list(size = 2)),
         shape = guide_legend(nrow=3, override.aes=list(size = 2, 
                                                        color = c('#212121','#9E9E9E')))))

fig3b <- plot_grid(fig3b.tree, fig3b.leg,
                   ncol = 1, rel_heights = c(5,1))

fig3c <- readPNG(sprintf('%sFig3C.png',outsidepanels))
fig3c <- ggplot() + 
  background_image(fig3c) +
  theme(plot.margin = margin(t=0, l=6, r=6, b=2.5, unit = "mm"),
        plot.background = element_blank())

fig3c.leg <- g_legend(data.frame(labels = c('Met15','Yll058w','MetY')) %>%
                        ggplot(aes(x = 1, y = labels, col = labels)) +
                        geom_point(shape = 15) +
                        scale_color_manual(name = '',
                                           values = c('Met15' = "#330099",
                                                      'Yll058w' = "#CC33FF",
                                                      'MetY' = 'grey40')) +
                        theme_linedraw() +
                        theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
                              axis.title = element_text(size = titles),
                              axis.text = element_text(size = txt),
                              legend.title =  element_blank(),
                              legend.text =  ggtext::element_markdown(size = txt),
                              legend.spacing = unit(0.1,"mm"),
                              legend.key.size = unit(4, "mm"),
                              legend.direction = 'horizontal',
                              legend.position = 'bottom') +
                        guides(col = guide_legend(nrow=1, override.aes=list(size = 3))))


fig3d <- data.pv3 %>%
  filter(hours == 120, id %in% c('MM')) %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  geom_boxplot(outlier.shape = NA, fill = '#9E9E9E', size = 0.3) +
  scale_x_discrete(limits = c(
    'FY4_empty',
    # 'FY4_ylldel_empty',
    'FY4_met15del_empty',
    'FY4_met15ylldel_empty',
    # 'FY4_ylldel_yll',
    'FY4_met15del_yll',
    'FY4_met15ylldel_yll',
    'FY4_met15del_yll_k376a',
    'FY4_met15ylldel_yll_k376a'),
    labels = c('FY4_empty' = 'FY4<br />w/ Empty<br />Plasmid',
               'FY4_ylldel_empty' = 'FY4-*yll058wΔ*<br />w/ Empty<br />Plasmid',
               'FY4_met15del_empty' = 'FY4-*met15Δ*<br />w/ Empty<br />Plasmid',
               'FY4_met15ylldel_empty' = 'FY4-*met15Δ yll058wΔ*<br />w/ Empty<br />Plasmid',
               'FY4_ylldel_yll' = 'FY4-*yll058wΔ*<br />w/ YLL058W<br />Plasmid',
               'FY4_met15del_yll' = 'FY4-*met15Δ*<br />w/ YLL058W<br />Plasmid',
               'FY4_met15ylldel_yll' = 'FY4-*met15Δ yll058wΔ*<br />w/ YLL058W<br />Plasmid',
               'FY4_met15del_yll_k376a' = 'FY4-*met15Δ*<br />w/ YLL058W<br />(K376A) Plasmid',
               'FY4_met15ylldel_yll_k376a' = 'FY4-*met15Δ yll058wΔ*<br />w/ YLL058W<br />(K376A) Plasmid'),
    position = "top") +
  labs(y = 'Relative Colony Size in\nSD-Met-Cys+Gal+G418') +
  theme_linedraw() +
  theme(plot.margin = margin(0,3,10,0),
        plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        panel.background = element_rect(fill = pinningpanelbg),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.top = ggtext::element_markdown(size = txt-1),
        axis.ticks.x.top = element_line(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.2))

fig3e <- oah_data %>%
  ggplot(aes(x = OAH, y = Rate)) +
  stat_summary(aes(group = Sample, fill = Sample), fun.data=mean_se, fun.args = list(mult=1), geom="ribbon",
               alpha = 0.4) +
  stat_summary(aes(group = Sample, col = Sample, linetype = Sample), fun=mean, geom="line", lwd = 0.7) +
  # stat_summary(aes(group = Sample, shape = Sample), fun=mean, geom="point", size =2) +
  stat_summary(aes(group = Sample, col = Sample, shape = Sample), fun=mean, geom="point", size = 2) +
  annotate("text", x = 10.1, y = 0.25, 
           label = "Met15:V[max]==list(89.210,K[m]==8.033)",
           size = 2, parse = T, hjust = 0) +
  annotate("text", x = 10.1, y = 0.25/2, 
           label = "Yll058w:V[max]==list(1.624,K[m]==4.292)",
           size = 2, parse = T, hjust = 0) +
  annotate("text", x = 10.1, y = 0.25/4, 
           label = "Yll058w (K376A):V[max]==list(0.003,K[m]<0.001)",
           size = 2, parse = T, hjust = 0) +
  labs(x = 'OAH (mM)',
       y = 'Rate (μM/min)') +
  scale_color_manual(name = 'Sample',
                     values = c('Met15' = "#330099", #6600CC
                                'Yll058w' = "#CC33FF",
                                'Yll058w (K376A)' = "#CC33FF"),
                     guide = 'none') +
  scale_fill_manual(name = 'Sample',
                    values = c('Met15' = "#330099",
                               'Yll058w' = "#CC33FF",
                               'Yll058w (K376A)' = "#CC33FF"),
                    guide = 'none') +
  scale_linetype_manual(name = 'Sample',
                        values = c('Met15' = 1,
                                   'Yll058w' = 1,
                                   'Yll058w (K376A)' = 6)) +
  scale_shape_manual(name = 'Sample',
                     values = c('Met15' = 16,
                                'Yll058w' = 16,
                                'Yll058w (K376A)' = 15)) +
  scale_y_continuous(trans = 'log2', breaks = c(0.016,0.25,4,64)) +
  theme_linedraw() +
  theme(plot.margin = margin(0,7,5,2),
        plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_text(margin = margin(-2,0,0,0)),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        legend.margin = margin(0,0,-5,0),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(linetype = guide_legend(keywidth = 1.5,
                                 override.aes = list(color = c("#330099", "#CC33FF", "#CC33FF"))))


fig3 <- plot_grid(plot_grid(plot_grid(fig3a, fig3c, fig3c.leg, 
                                      ncol = 1, rel_heights = c(0.8,1,0.05),
                                      labels = c('A','C',''),
                                      label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                            fig3b, 
                            nrow = 1, rel_widths = c(1.5,1),
                            labels = c('','B'),
                            label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  NULL,
                  plot_grid(fig3d, fig3e,
                            ncol = 1, rel_heights = c(0.9,1),
                            labels = c('','E'),
                            label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  ncol = 1, rel_heights = c(2.2,0.6,1.9),
                  labels = c('','D',''),
                  label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')

ggsave(sprintf("%s/Figure3.jpg",fig_path), fig3,
       height = 230, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


##### FIGURE 4: Chelation of H2S facilitates the growth of met15Δ cells in the absence of exogenous organosulfurs.
# fig4a <- readPNG(sprintf('%sFig4A.png',outsidepanels))
# fig4a <- ggplot() + 
#   background_image(fig4a) +
#   theme(plot.margin = margin(t=0, l=6, r=6, b=0, unit = "mm"),
#         plot.background = element_blank())

fig4b.ori <- rbind(temp1, temp2) %>%
  filter(escapee == 'No') %>%
  ggplot(aes(x = FeEDTA, y = OD, fill = FeEDTA)) +
  stat_summary(col = 'black',
               alpha = 0.9, size = 0.3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = 't.test', label = 'p.format',
                     paired = F, method.args = list(var.equal = FALSE),
                     label.x = 1.5, label.y = 7.9,
                     hjust = 0.5, size = 2.5) +
  scale_fill_manual(guide = 'none',
                    values = c('Absent' = '#D6E4F4',
                               'Present' = '#71B3AF')) +
  scale_color_manual(guide = 'none',
                    values = c('Absent' = '#D6E4F4',
                               'Present' = '#71B3AF')) +
  labs(x = 'H<sub>2</sub>S Chelator',
       y = 'OD<sub>600</sub> at Saturation<br/>of Original Cultures') +
  facet_grid(.~Strain, labeller = labeller(Strain = c('FY4' = 'FY4',
                                                      'met15del' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.margin = margin(1,1,1,1),
        plot.title = element_blank(),
        panel.background = element_rect(fill = flaskpanelbg, color = 'transparent'),
        # axis.title.x = ggtext::element_markdown(size = txt),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                                margin = margin(1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,9))

fig4b.re <- rbind(temp3, temp4) %>%
  filter(escapee == 'No') %>%
  ggplot(aes(x = FeEDTA, y = OD, fill = FeEDTA)) +
  stat_summary(col = 'black', alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  stat_compare_means(method = 't.test', label = 'p.format',
                     paired = F, method.args = list(var.equal = FALSE),
                     label.x = 1.5, label.y = 7.9,
                     hjust = 0.5, size = 2.5) +
  scale_x_discrete(labels = c('Absent' = 'Absent\nPreviously',
                              'Present' = 'Present\nPreviously')) +
  scale_fill_manual(guide = 'none',
                    values = c('Absent' = '#D6E4F4',
                               'Present' = '#D6E4F4')) +
  labs(x = 'H<sub>2</sub>S Chelator',
       y = 'OD<sub>600</sub> at Saturation<br/>of Repassaged Cultures') +
  facet_grid(.~ORF, labeller = labeller(ORF = c('FY4' = 'FY4',
                                                'met15D' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.margin = margin(1,1,1,1),
        plot.title = element_blank(),
        panel.background = element_rect(fill = flaskpanelbg, color = 'transparent'),
        axis.title.x = ggtext::element_markdown(size = txt),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                                margin = margin(1,0,0.1,0, "mm")))  +
  coord_cartesian(ylim = c(0,9))

# fig4c <- readPNG(sprintf('%sFig4C.png',outsidepanels))
# fig4c <- ggplot() + 
#   background_image(fig4c) +
#   theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
#         plot.background = element_blank())

fig4d <- data.chelator.dose %>%
  filter(petite == 'Yes') %>%
  ggplot(aes(x = Media_ID, y = OD,
             fill = FeEDTA)) +
  stat_summary(col = 'black', alpha = 0.9, size = 0.3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = 't.test', label = 'p.format',
                     paired = F, method.args = list(var.equal = FALSE),
                     label.x = 1.5, label.y = 6,
                     hjust = 0.5, size = 2.5) +
  scale_x_discrete(labels = c('SD-Met+Glu' = 'Absent',
                              'SD-Met+Glu+0.048M FeEDTA' = 'Present')) +
  scale_fill_gradient(guide = 'none',
                    low = '#D6E4F4',
                    high = '#71B3AF') +
  facet_wrap(.~labels, labeller = label_parsed,
             ncol = 1) +
  labs(x = 'H<sub>2</sub>S Chelator',
       y = 'OD<sub>600</sub> at Saturation') +
  theme_linedraw() +
  theme(plot.margin = margin(5,1,5,5),
        plot.title = element_blank(),
        panel.background = element_rect(fill = flaskpanelbg, color = 'transparent'),
        axis.title.x = ggtext::element_markdown(size = txt),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text.x = element_text(size = txt),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = element_text(size = txt,
                                    margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,6.5))


fig4 <- plot_grid(plot_grid(NULL, plot_grid(fig4b.ori, fig4b.re, ncol = 1, rel_heights = c(1,1.2),
                                             labels = c('B',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                            nrow = 1, rel_widths = c(1.25,1), 
                            labels = c('A',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  plot_grid(NULL, fig4d, nrow = 1, rel_widths = c(3,1),
                            labels = c('C','D'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  ncol = 1, rel_heights = c(1,1.2),
                  labels = c('',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure4.jpg",fig_path), fig4,
       height = 163, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)

# 6.85/two.c * 67

##### FIGURE S1: Growth of met15Δ cells in an organosulfur-deficient medium at earlier incubation times. 

##### FIGURE S2: met15Δ-associated growth phenotype segregates 2:2 in a heterozygous cross. 

##### FIGURE S3: Robust growth of met15Δ cells in media lacking organosulfurs
figS3a <- rbind(data.cbn[,c(1,2,5,6,10)],
                  data.cbn.leu[,c(8,11,10,13,14)] %>%
                    filter(condition == 'SCmLeu')) %>%
  filter(base == 'SC', orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  geom_boxplot(fill = '#9E9E9E', size = 0.2, outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(limits = c('BY4742','BY4741','FY4','FY4-met15del'),
                   labels = c('FY4' = 'FY4',
                              'FY4-met15del' = 'FY4-*met15Δ*',
                              'BY4742' = 'BY4742',
                              'BY4741' = 'BY4741')) +
  labs(x = 'Strain', y= 'Relative Colony Size\nat Saturation') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_wrap(.~condition, scales = 'free_x',
             nrow = 1,
             labeller = labeller(condition = c('GLU' = 'SC-Met-Cys+Glu',
                                               'SCmLeu' = 'SC+Met-Cys-Leu+Glu',
                                               'GAL' = 'SC-Met-Cys+Gal'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, hjust = 0.5, face = 'bold'),
        panel.background = element_rect(fill = pinningpanelbg, color = 'transparent'),
        axis.title.x = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(), #angle = 30, vjust = 1, hjust = 1
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = element_text(size = txt,  margin = margin(0.1,0,0.1,0, "mm")))

figS3b <- rbind(data.cbn[,c(1,2,5,6,10)],
                  data.cbn.leu[,c(8,11,10,13,14)] %>%
                    filter(condition == 'SDmLeu')) %>%
  filter(base == 'SD', orf_name %in% c('FY4-met15del')) %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  # geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(fill = '#9E9E9E', size = 0.2, outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
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
  labs(x = 'SD-Cys Media', y= 'Relative Colony Size\nat Saturation') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_grid(orf_name ~ carbon, scales = 'free_x', space = 'free_x',
             labeller = labeller(orf_name = c('FY4' = 'FY4',
                                              'FY4-met15del' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.margin = margin(5,5,1,5),
        plot.title = ggtext::element_markdown(size = titles, hjust = 0.5, face = 'bold'),
        panel.background = element_rect(fill = pinningpanelbg, color = 'transparent'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = element_text(size = txt,  margin = margin(0.1,0,0.1,0, "mm")),
        strip.text.y = ggtext::element_markdown(size = txt, colour = 'white',
                                                margin = margin(0.1,1,0.1,0, "mm")))


figS3 <- plot_grid(figS3a,
                   figS3b,
                   labels = c('A','B'),
                   label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
                   ncol = 1, rel_heights = c(1,1))
ggsave(sprintf("%s/FigureS3.jpg",fig_path), figS3,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


##### FIGURE S4: The growth of FY4-met3Δ cells on SD-Met media is generally unaffected by the carbon source
figS4 <- rbind(data.cbn[,c(1,2,5,6,10)],
                     data.cbn.leu[,c(8,11,10,13,14)] %>%
                       filter(condition == 'SDmLeu')) %>%
  filter(base == 'SD', orf_name == 'FY4-met3del') %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  # geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(fill = '#9E9E9E', size = 0.2, outlier.shape = NA) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
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
  labs(title = 'FY4-*met3Δ* in Synthetic Defined Media',
       x = 'SD Media', y= 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_grid(orf_name~carbon, scales = 'free_x', space = 'free_x',
             labeller = labeller(orf_name = c('FY4-met3del' = 'FY4-*met3Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        panel.background = element_rect(fill = pinningpanelbg, color = 'transparent'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,  margin = margin(0.1,0,0.1,0, "mm")),
        strip.text.y = ggtext::element_markdown(size = txt, margin = margin(0.1,1,0.1,0, "mm")))

ggsave(sprintf("%s/FigureS4.jpg",fig_path), figS4,
       height = one.c, width = two.c, units = 'mm',
       dpi = 300)


##### FIGURE S5: Relative H2S levels of mutant strains
figS5a <- readPNG(sprintf('%sFigS5A.png',outsidepanels))
figS5a <- ggplot() + 
  background_image(figS5a) +
  theme(plot.margin = margin(t=10, l=0, r=0, b=10, unit = "mm"),
        plot.background = element_rect(fill = 'white'))

figS5b <- h2s.col %>%
  filter(strain %notin% c('BY4742','BY4741','FY4-*yll058wΔ*')) %>%
  ggplot(aes(x = condition, y = relative_color_intensity)) +
  stat_summary(data = h2s.col.sum %>%
                 filter(strain %notin% c('BY4742','BY4741','FY4-*yll058wΔ*')),
               aes(fill = relative_color_intensity), col = 'transparent', alpha = 1, size = 1,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  stat_compare_means(method = 'kruskal', size = 1.8) +
  scale_fill_gradient2(high = "#5D4037", low = "white", guide = "none") +
  scale_x_discrete(labels = c('BiGGY' = 'BiGGY',
                              'SD-Met-Cys+Glu+Bi' = 'SD-Met\n-Cys+Glu\n+Bi')) +
  scale_y_continuous(breaks = seq(0,1000,150)) +
  labs(y = 'Color Intensity (~H<sub>2</sub>S levels)',
       x = 'Strains') +
  facet_wrap(.~strain, nrow = 2) +
  theme_linedraw() +
  theme(plot.title = element_blank(), #element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(1,0,0.1,0, "mm")))


figS5 <- plot_grid(annotate_figure(figS5a, top = text_grob("SD-Met-Cys+Glu+Bi Media", face = "bold", size = titles)),
                   figS5b,
                   ncol = 1, rel_heights = c(1,2),
                   labels = c('A','B'),
                   label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/FigureS5.jpg",fig_path), figS5,
       height = 200, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


##### FIGURE S6: YLL058W homologs cluster into three distinct classes.
figS6 <- ggplot(df_yll)+
  geom_point(aes(x=PC1,y=PC2,col=km_clust,shape=yll_taxa),
             size = 2)+
  theme(legend.title=element_blank())+
  scale_color_manual(values=c("Cluster 1"="#CC33FF",
                              "Cluster 2"="#99CCFF",
                              "Cluster 3"="#006666"),
                     labels = c('Cluster 1' = 'YLL058W class',
                                'Cluster 2' = 'Ancestral class',
                                'Cluster 3' = '*STR2*/YML082W class'),
                     limits = c('Cluster 2', 'Cluster 1', 'Cluster 3'))+ #
  scale_shape_manual(values=c(20, 0, 2))+
  labs(x="PC1",y="PC2")+
  annotate("segment", x = df_yll[331,]$PC1, xend = df_yll[331,]$PC1, y = df_yll[331,]$PC2, yend = df_yll[331,]$PC2+.7)+
  annotate(geom="text", x=df_yll[331,]$PC1, y=df_yll[331,]$PC2+.9, label="STR2", color="black", size = 2,fontface = "italic")+
  annotate("segment", x = df_yll[330,]$PC1, xend = df_yll[330,]$PC1-.5, y = df_yll[330,]$PC2, yend = df_yll[330,]$PC2+.5)+
  annotate(geom="text", x=df_yll[330,]$PC1-.5, y=df_yll[330,]$PC2+.7, label="YML082W", color="black", size = 2)+
  annotate("segment", x = df_yll[329,]$PC1, xend = df_yll[329,]$PC1+.5, y = df_yll[329,]$PC2, yend = df_yll[329,]$PC2)+
  annotate(geom="text", x=df_yll[329,]$PC1+0.8, y=df_yll[329,]$PC2, label="YLL058W", color="black", size = 2) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_blank(),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/FigureS6.jpg",fig_path), figS6,
       height = one.5c, width = one.5c, units = 'mm',
       bg = 'white',
       dpi = 300)

##### FIGURE S7: Yll058w functions as a monomer, in contrast to the tetrameric conformation seen in canonical homocysteine synthases

##### FIGURE S8: Deletion of YLL058W has negligible effect on fitness and H2S accumulation in the FY4 background
figS8a <- data.jm.2 %>%
  filter(orf_name %in% c('yll')) %>%
  ggplot() +
  geom_boxplot(aes(x = condition, y = relative_fitness),
               fill = '#9E9E9E', outlier.shape = NA, size = 0.2,
               position = position_dodge(width = 0.85)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(limits = c('YPDA', 'SD-Met'),
                   labels = c('YPDA' = 'YPDA',
                              'SD-Met' = 'SD-Met-Cys+Glu')) +
  scale_fill_gradient(name = 'Relative Hydrogen Sulfide',
                      low = "#D7CCC8", high = "#5D4037",
                      limits = c(1,7)) +
  labs(x = 'Strain',
       y = 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        panel.background = element_rect(fill = pinningpanelbg, color = 'transparent'),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt),#, angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        # legend.direction = 'vertical',
        # legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(1,0,0.1,0, "mm")))

figS8b <- h2s.col %>%
  filter(strain %in% c('FY4-*yll058wΔ*')) %>%
  ggplot(aes(x = condition, y = relative_color_intensity)) +
  stat_summary(data = h2s.col.sum %>%
                 filter(strain %in% c('FY4-*yll058wΔ*')),
               aes(fill = relative_color_intensity), col = 'transparent', alpha = 1, size = 1,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  stat_compare_means(method = 'kruskal', size = 2, angle = 270, hjust = 1) +
  scale_fill_gradient2(high = "#5D4037", low = "white", guide = "none") +
  labs(y = 'Color Intensity (~H<sub>2</sub>S levels)',
       x = '') +
  scale_x_discrete(limits = c('BiGGY','SD-Met-Cys+Glu+Bi'),
                   labels = c('BiGGY','SD-Met\n-Cys+Glu\n+Bi'),
                   position = 'top') +
  theme_linedraw() +
  theme(plot.title = element_blank(), #element_text(size = titles, face = 'bold', hjust = 0.5),
        plot.margin = margin(5,0,5,15),
        axis.title.y = element_blank(),
        axis.title.x = ggtext::element_markdown(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'black',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_flip(ylim = c(0,500))

figS8 <- annotate_figure(plot_grid(figS8a, figS8b, NULL,
                   ncol = 3, rel_widths = c(0.9,1,0.5),
                   align = 'h',
                   labels = c('A','B',''),
                   label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                   top = text_grob(expression(paste("FY4-",italic("yll058wΔ"))), size = titles))
ggsave(sprintf("%s/FigureS8.jpg",fig_path), figS8,
       height = one.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


##### FIGURE S9: Coomassie-stained SDS-PAGE gel showing the final fractions that contained the recombinantly-expressed,
##### purified protein used in the in vitro homocysteine biosynthesis assay. 

##### FIGURE S10: Yll058w catalyzes homocysteine biosynthesis but less efficient than Met15. 
data.bioc %>%
  filter(`Time (min)` == 15) %>%
  group_by(Sample) %>%
  summarize(uM = median(uM, na.rm = T), .groups = 'keep')

figS10a <- data.bioc %>%
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
  theme(
    # plot.margin = margin(0,7,5,2),
    plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
    axis.title = element_text(size = titles),
    # axis.title.x = element_text(margin = margin(-2,0,0,0)),
    axis.text = element_text(size = txt),
    legend.title = element_blank(),
    legend.text = element_text(size = txt),
    legend.position = 'bottom',
    legend.key.size = unit(3, "mm"),
    legend.box.spacing = unit(0.5,"mm"),
    legend.margin = margin(0,0,-5,0),
    strip.text = element_text(size = txt,
                              face = 'bold',
                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(xlim = c(0,16),
                  ylim = c(5,125)) +
  guides(fill = guide_legend(override.aes=list(shape = 15, alpha = 1)))

figS10b <- data.mdl.2 %>%
  ggplot(aes(x = log(Yll058w_Met15,10), y = Biomass_Flux)) +
  geom_line(aes(col = Model), lwd = 1) +
  geom_point(size = 2) +
  geom_point(aes(col = Model), size = 0.5) +
  scale_x_continuous(trans = 'reverse',
                     breaks = seq(10,-10,-1),
                     labels = 10^seq(10,-10,-1)) +
  scale_color_manual(name = 'Model with',
                     values = c('C' = '#333333',
                                'D' = '#999999'),
                     labels = c('C'='Hypothesized YLL058W reaction<br/>w/ all *MET15* reactions',
                                'D'='Hypothesized YLL058W reaction<br/>w/o all *MET15* reactions')) +
  labs(x = 'Yll058w/Met15 Kcat Ratio',
       y = 'Simulated Biomass Flux\n("Growth")') +
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

figS10 <- plot_grid(figS10a, figS10b,
                   ncol = 1, rel_heights = c(1,1),
                   align = 'v',
                   labels = c('A','B'),
                   label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/FigureS10.jpg",fig_path), figS10,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)

##### FIGURE S11: A growth defect of met15Δ cells grown with ethanol is seen in the presence and absence of organosulfurs in the media
figS11 <- merge(data.res.gc[str_detect(data.res.gc$orf_name, 'FY'),],
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
                  parse = T, size = 2.2, min.segment.length = 10) +
  facet_grid(~carbon*methionine,
             labeller = labeller(methionine = c('+Met' = 'SD+Met-Cys',
                                                '-Met' = 'SD-Met-Cys'))) +
  scale_color_manual(name = 'Presumed Auxotroph',
                     values = c('Prototroph' = '#FFC107',
                                'Presumed Auxotroph' = '#536DFE',
                                'Uracil-Leucine' = '#E040FB',
                                'Methionine-Uracil-Leucine' = '#FF5722'),
                     limits = c('Presumed Auxotroph','Prototroph'),
                     labels = c('Prototroph'='No',
                                'Presumed Auxotroph'='Yes',
                                'Uracil-Leucine',
                                'Methionine-Uracil-Leucine')) +
  scale_fill_manual(name = 'Auxotrophy',
                    values = c('Prototroph' = '#FFC107',
                               'Presumed Auxotroph' = '#536DFE'),
                    limits = c('Prototroph','Presumed Auxotroph'),
                    guide = 'none') +
  scale_linetype_manual(name = 'Petite',
                        values = c('Yes' = 'dotdash',
                                   'No' = 'solid'),
                        limits = c('Yes','No')) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  labs(x = 'Time (hours)', y= 'Relative Colony Size') +
  coord_cartesian(xlim = c(0,350),
                  ylim = c(0,1.3)) +
  theme_linedraw() +
  theme(plot.margin = margin(1,5,1,5),
        plot.title = element_blank(),
        panel.background = element_rect(fill = pinningpanelbg, color = 'transparent'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.margin = margin(0,0,0,0),
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
         linetype = guide_legend(nrow=1, byrow=TRUE, order = 2, length = 2))

ggsave(sprintf("%s/FigureS11.jpg",fig_path), figS11,
       height = 67, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)

##### FIGURE S12: Titration of the H2S chelator Fe-EDTA.
unique(data.chelator.dose$Media_ID)
figS12 <- data.chelator.dose %>%
  filter(Hours == 53, petite == 'No') %>%
  ggplot(aes(x = Media_ID, y = OD,
             fill = as.factor(Media_ID))) +
  stat_summary(col = 'black', alpha = 0.9, size = 0.3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  scale_fill_manual(name = 'Chelator\nConcentration\n(M)',
                    values = c('SD+Met+Glu' = '#FFFFFF',
                               'SD-Met+Glu' = '#D6E4F4',
                               'SD-Met+Glu+0.024M FeEDTA' = '#A9CAC8',
                               'SD-Met+Glu+0.048M FeEDTA' = '#71B3AF',
                               'SD-Met+Glu+0.096M FeEDTA' = '#6AA3A0',
                               'SD-Met+Glu+0.136M FeEDTA' = '#138F87'),
                    guide = 'none') +
  scale_x_discrete(labels = c('SD+Met+Glu' = 'SD+Met-Cys+Glu',
                              'SD-Met+Glu' = 'SD-Met-Cys+Glu',
                              'SD-Met+Glu+0.024M FeEDTA' = 'SD-Met-Cys+Glu+0.024M FeEDTA',
                              'SD-Met+Glu+0.048M FeEDTA' = 'SD-Met-Cys+Glu+0.048M FeEDTA',
                              'SD-Met+Glu+0.096M FeEDTA' = 'SD-Met-Cys+Glu+0.096M FeEDTA',
                              'SD-Met+Glu+0.136M FeEDTA' = 'SD-Met-Cys+Glu+0.136M FeEDTA')) +
  facet_grid(.~Strain, scales = 'free_x', space = 'free_x',
             labeller = labeller(Strain = c('FY4' = 'FY4',
                                            'FY4-met15del' = 'FY4-*met15Δ*'))) +
  labs(y = 'OD<sub>600</sub> at Saturation') +
  theme_linedraw() +
  theme(plot.margin = margin(5,5,5,15),
        plot.title = element_blank(),
        panel.background = element_rect(fill = flaskpanelbg, color = 'transparent'),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text.x = element_text(size = txt, angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'right',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                                margin = margin(1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,6.5))
ggsave(sprintf("%s/FigureS12.jpg",fig_path), figS12,
       height = one.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)
