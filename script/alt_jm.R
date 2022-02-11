
temp1 <- merge(data.510, strain.labs.510, by = 'orf_name')
# %>%
#   filter(hours == 160, condition %in% c('YPDA','SD-Met'), orf_name %in% c('met5del','met10del'))
temp1 <- temp1[((temp1$hours == 89 & temp1$condition == 'YPDA') |
                 (temp1$hours == 160 & temp1$condition == 'SD-Met')) &
                 temp1$orf_name %in% c('met5del','met10del'),]
temp2 <- merge(data.jm, strain.labs.jm, by = 'orf_name') %>%
  filter(time == 't_final', attempt != 'pilot', condition %in% c('YPDA','SD-Met-Cys+Glu'))
temp2 <- temp2[,c(1:9,13,10,11,15,16)]
colnames(temp2) <- c(colnames(temp2)[1:9],'relative_fitness',colnames(temp2)[11:14])

data.jm.2 <- rbind(temp1, temp2)
data.jm.2$condition[data.jm.2$condition == 'SD-Met-Cys+Glu'] <- 'SD-Met'
data.jm.2$orf_name <- factor(data.jm.2$orf_name, levels = c('FY4','met15','met3','met5del','met10del','met2',
                                                            'met6','met13','cys4','str3','met12','yll'))

head(strain.labs.jm)

strain.labs.jm <- rbind(strain.labs.jm, 
                        data.frame(orf_name = c('met5del','met10del'), labels = c('FY4-*met5Δ*','FY4-*met10Δ*'),
                                   parsed = c('FY4-italic(met5Δ)','FY4-italic(met10Δ)'), auxotrophy = c('Prototroph','Prototroph')))


data.jm.2 %>%
  filter(orf_name %in% c('met5del','met10del','FY4')) %>%
  group_by(condition, orf_name) %>%
  summarize(f = median(relative_fitness, na.rm = T), .groups = 'keep') %>%
  data.frame()

0.1132472 - 1
0.1310366 - 1

data.jm.2 %>%
  group_by(condition, orf_name) %>%
  summarize(max_hrs = max(hours), .groups = 'keep') %>%
  data.frame()

fig3a <- data.jm.2[!(data.jm.2$orf_name %in% c('yll','BY4742','BY4741')),] %>%
  ggplot() +
  geom_boxplot(aes(x = condition, y = relative_fitness),
               fill = '#9E9E9E', outlier.shape = NA, size = 0.2,
               position = position_dodge(width = 0.85)) +
  scale_y_continuous(breaks = seq(-1,2,0.4)) +
  scale_x_discrete(limits = c('YPDA', 'SD-Met'),
                   labels = c('YPDA' = 'YPDA',
                              'SD-Met' = 'SD-Met+Glu')) +
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

fig3a.g <- ggplot_gtable(ggplot_build(fig3a))
stripr <- which(grepl('strip-t', fig3a.g$layout$name))

for (s in stripr) {
  # fig3b.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- 'black'
  l <- fig3a.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[1]]$label
  
  if (length(l) != 0) {
    if (l == 'FY4-') {
      l2 <- fig3a.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$label
      a <- strain.labs.jm$auxotrophy[str_detect(strain.labs.jm$labels, l2)]
    } else {
      a <- strain.labs.jm$auxotrophy[strain.labs.jm$orf_name == l]
    }
    
    if (a == 'Presumed Auxotroph') {
      fig3a.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#536DFE'
    } else if (a == 'Prototroph') {
      fig3a.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107'
      # fig3a.g$grobs[[s]]$grobs[[1]]$children[[2]]$children[[1]]$children[[1]]$children[[2]]$gp$col <- 'black'
    } else {
      fig3a.g$grobs[[s]]$grobs[[1]]$children[[1]]$gp$fill <- '#FFC107' #Unkown = '#BDBDBD'
    }
  }
}



data.bis$Strain <- factor(data.bis$Strain,
                          levels = c('FY4','FY4-met15D','FY4-met3D','FY4-met5D','FY4-met10D','FY4-met2D',
                                     'FY4-met6D','FY4-met13D','FY4-cys4D','FY4-str3D','FY4-met12D','FY4-yllD'))

data.bis.2 <- data.frame(orf_name = c('met10del','met5del','FY4','yll',
                                      'met10del','met5del','FY4','yll'),
                         H2S = c(1,1,4,4,
                                 1,1,4,4))
data.bis <- rbind(data.bis, data.bis[c(64:69),])

fig3c <- data.bis[!(data.bis$Strain %in% c('BY4742','BY4741')),] %>%
  filter(Condition == 'BiGGY', Strain != 'FY4-yllD') %>%
  ggplot(aes(x = Strain, y = HS)) +
  stat_summary(data = data.bis[!(data.bis$Strain %in% c('BY4742','BY4741')),] %>%
                 filter(Condition == 'BiGGY', Strain != 'FY4-yllD') %>%
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
  labs(y = 'Relative\nHydrogen Sulfide',
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
                              'FY4-met5D'='FY4-*met5Δ*',
                              'FY4-met10D'='FY4-*met10Δ*',
                              'FY4-cys4D'='FY4-*cys4Δ*',
                              'FY4-yllD'='FY4-*yll058wΔ*',
                              'BY4742'='BY4742',
                              'BY4741'='BY4741')) +
  theme_linedraw() +
  theme(plot.title = element_blank(),#element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 15, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'black',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(1,7))

fig1a <- readPNG('figures/final/final/MetaPath.png')
fig1a <- ggplot() + 
  background_image(fig1a) +
  theme(plot.margin = margin(t=0, l=35, r=35, b=0, unit = "mm"),
        plot.background = element_blank())

fig3b <- readPNG('figures/final/final/ColonyCrops/BiGGYCrops.png')
fig3b <- ggplot() + 
  background_image(fig3b) +
  theme(plot.margin = margin(t=0, l=10, r=8, b=0, unit = "mm"),
        plot.background = element_blank())

figNSschema <- readPNG('figures/final/final/NoSulfatesSchema.png')
figNSschema <- ggplot() + 
  background_image(figNSschema) +
  theme(plot.margin = margin(t=0, l=35, r=20, b=0, unit = "mm"),
        plot.background = element_blank())


fig3 <- plot_grid(fig1a, fig3a.g, fig3b, fig3c,
                  plot_grid(fig3d.expt, fig3d,
                            ncol = 1, align = 'hv', axis = 'lr',
                            rel_heights = c(0.3,1)) , #fig3d is alt_no_sulfate,
                  ncol = 1, rel_heights = c(1.2,1,0.8,0.7,1.5),
                  labels = c('A','B','C','D','E'),
                  label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure3.jpg",fig_path), fig3,
       height = 250, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 600)
