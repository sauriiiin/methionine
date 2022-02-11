#####

strain.labs.h2s <- read.csv(file = 'data/h2s/FeEDTA_samples.csv')
data.h2s <- read.csv(file = 'data/h2s/FeEDTA_readings.csv')

data.h2s.re <- NULL
for (t in pzfx_tables("/home/sbp29/R/Projects/methionine/data/h2s/chelator_liquid.pzfx")) {
  temp <- read_pzfx("/home/sbp29/R/Projects/methionine/data/h2s/chelator_liquid.pzfx", table = t)
  temp <- melt(temp, id.vars = 'Time (Minutes)', variable.name = 'Strain', value.name = 'OD')
  temp$Set = t
  data.h2s.re <- rbind(data.h2s.re, temp)
}
data.h2s.re <- data.frame(data.h2s.re)
data.h2s.re <- data.h2s.re[data.h2s.re$Set == 'All_Samples',]

temp <- str_split(data.h2s.re$Strain, '_', simplify = T)
temp[,4] <- temp[!str_detect(temp[,3], 'R'),3]
temp[,3] <- temp[str_detect(temp[,2], 'R'),2]
temp[str_detect(temp[,2], 'R'),2] <- ''

temp <- data.frame(temp, stringsAsFactors = F)
colnames(temp) <- c('ORF','FeEDTA','Replicate','Reading')
temp$FeEDTA <- as.character(temp$FeEDTA)
temp$FeEDTA[temp$FeEDTA == 'Fe'] <- 'Present'
temp$FeEDTA[temp$FeEDTA == ''] <- 'Absent'

data.h2s.re <- cbind(data.h2s.re, temp)

strain.labs.h2sr <- read.csv(file = 'data/h2s/FeEDTA_repeat_samples.csv')
data.h2sr <- read.csv(file = 'data/h2s/FeEDTA_repeat_readings.csv')

#####

temp1 <- merge(data.h2s %>%
                 melt(id.vars = 'Time', variable.name = 'Flask', value.name = 'OD'), strain.labs.h2s, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Time == 3075, Expected.Initial.OD == 0.1)
temp1$escapee[temp1$Replicate == 'Rep_1' & temp1$FeEDTA == 'Absent' & temp1$Strain == 'met15del'] <- 'Yes'
temp1$escapee[is.na(temp1$escapee)] <- 'No'
# filter(Time == 3075, Expected.Initial.OD == 0.1, Replicate != 'Rep_1') 

temp2 <- merge(data.h2sr %>%
        melt(id.vars = c('Time', 'Attempt'), variable.name = 'Flask', value.name = 'OD') %>%
        filter(Attempt == 'Original'), strain.labs.h2sr, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Time == 3075, Expected.Initial.OD == 0.1, Strain %in% c('FY4','met15del'))
temp2 <- temp2[,-3]
temp2$escapee <- 'No'

rbind(temp1, temp2) %>%
  filter(escapee == 'No') %>%
  group_by(Strain, FeEDTA) %>%
  summarize(OD_m = mean(OD, na.rm = T), OD_sd = sd(OD, na.rm = T), .groups = 'keep')

fig6b <- rbind(temp1, temp2) %>%
  filter(escapee == 'No') %>%
  ggplot(aes(x = FeEDTA, y = OD)) +
  stat_summary(col = 'black', alpha = 0.9, size = 0.3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_compare_means(method = 't.test', label = 'p.format',
                     paired = F, method.args = list(var.equal = FALSE),
                     label.x = 1.5, label.y = 7.9,
                     hjust = 0.5, size = 2.5) +
  labs(x = 'H<sub>2</sub>S Chelator',
       y = 'OD<sub>600</sub> at Saturation<br/>of Original Cultures') +
  facet_grid(.~Strain, labeller = labeller(Strain = c('FY4' = 'FY4',
                                                      'met15del' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
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



temp3 <- data.h2s.re[data.h2s.re$Time..Minutes. == 3075,] %>%
  group_by(FeEDTA, ORF, Replicate) %>%
  summarize(OD = median(OD, na.rm = T), .groups = 'keep') %>%
  data.frame()
temp3$escapee[temp3$FeEDTA == 'Absent' & temp3$Replicate == 'R2' & temp3$ORF == 'met15D'] <- 'Yes'
temp3$escapee[temp3$FeEDTA == 'Present' & temp3$Replicate == 'R3' & temp3$ORF == 'met15D'] <- 'Yes'
temp3$escapee[is.na(temp3$escapee)] <- 'No'

temp4 <- merge(data.h2sr %>%
        melt(id.vars = c('Time', 'Attempt'), variable.name = 'Flask', value.name = 'OD') %>%
        filter(Attempt == 'Repassage'), strain.labs.h2sr, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Time == 3075, Expected.Initial.OD == 0.1, Strain %in% c('FY4','met15del')) %>%
  group_by(FeEDTA, Strain, Replicate) %>%
  summarise(OD = OD, .groups = 'keep') %>%
  data.frame()
colnames(temp4) <- colnames(temp3)[-5]
temp4$ORF <- as.character(temp4$ORF)
temp4$ORF[temp4$ORF == 'met15del'] <- 'met15D'
temp4$Replicate <- as.character(temp4$Replicate)
temp4$Replicate <- str_replace(temp4$Replicate, 'Rep_','R')
temp4$escapee[temp4$FeEDTA == 'Present' & temp4$Replicate %in% c('R1','R3') & temp4$ORF == 'met15D'] <- 'Yes'
temp4$escapee[is.na(temp4$escapee)] <- 'No'

rbind(temp3, temp4) %>%
  filter(escapee == 'No', ORF == 'FY4', FeEDTA == 'Present') %>%
  group_by(FeEDTA)


fig6c <- rbind(temp3, temp4) %>%
  filter(escapee == 'No') %>%
  ggplot(aes(x = FeEDTA, y = OD)) +
  stat_summary(col = 'black', alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  stat_compare_means(method = 't.test', label = 'p.format',
                     paired = F, method.args = list(var.equal = FALSE),
                     label.x = 1.5, label.y = 7.9,
                     hjust = 0.5, size = 2.5) +
  scale_x_discrete(labels = c('Absent' = 'Absent\nPreviously',
                              'Present' = 'Present\nPreviously')) +
  labs(x = 'H<sub>2</sub>S Chelator',
       y = 'OD<sub>600</sub> at Saturation<br/>of Repassaged Cultures') +
  facet_grid(.~ORF, labeller = labeller(ORF = c('FY4' = 'FY4',
                                                'met15D' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
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

fig2c <- merge(data.res.gc[str_detect(data.res.gc$orf_name, 'FY'),],
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
                  parse = T, size = 2.5, min.segment.length = 10) +
  facet_grid(~carbon*methionine,
             labeller = labeller(methionine = c('+Met' = 'SD + Met',
                                                '-Met' = 'SD - Met'))) +
  scale_color_manual(name = 'Presumed Auxotroph',
                     values = c('Prototroph' = '#FFC107',
                                'Presumed Auxotroph' = '#536DFE',
                                'Uracil-Leucine' = '#E040FB',
                                'Methionine-Uracil-Leucine' = '#FF5722'),
                     limits = c('Prototroph','Presumed Auxotroph'),
                     labels = c('Prototroph'='No',
                                'Presumed Auxotroph'='Yes',
                                'Uracil-Leucine',
                                'Methionine-Uracil-Leucine')) +
  scale_fill_manual(name = 'Auxotrophy',
                    values = c('Prototroph' = '#FFC107',
                               'Presumed Auxotroph' = '#536DFE'),
                    limits = c('Prototroph','Presumed Auxotroph'),
                    guide = F) +
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
        legend.direction = 'horizontal',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=1, byrow=TRUE, order = 1,
                              override.aes = list(size = 3)),
         linetype = guide_legend(nrow=1, byrow=TRUE, order = 2))



fig6a <- readPNG('figures/final/final/Fig4A.png')
fig6a <- ggplot() + 
  background_image(fig6a) +
  theme(plot.margin = margin(t=0, l=10, r=10, b=0, unit = "mm"),
        plot.background = element_blank())

fig6d <- readPNG('figures/final/final/Fig4C.png')
fig6d <- ggplot() + 
  background_image(fig6d) +
  theme(plot.margin = margin(t=0, l=20, r=20, b=0, unit = "mm"),
        plot.background = element_blank())


# fig6 <- plot_grid(plot_grid(fig6a, plot_grid(fig6b, fig6c, ncol = 1, rel_heights = c(1,1.2),
#                                              labels = c('B',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
#                             nrow = 1, rel_widths = c(1.5,1), 
#                             labels = c('A',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
#                   fig6d, ncol = 1, rel_heights = c(1,1.3),
#                   labels = c('','C'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
# ggsave(sprintf("%s/Figure6.jpg",fig_path), fig6,
#        height = two.c, width = two.c, units = 'mm',
#        dpi = 600)


fig5 <- plot_grid(plot_grid(fig6a, plot_grid(fig6b, fig6c, ncol = 1, rel_heights = c(1,1.2),
                                              labels = c('B',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                             nrow = 1, rel_widths = c(1.5,1), 
                             labels = c('A',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                   fig6d,
                   fig2c,
                   ncol = 1, rel_heights = c(1,1.2,0.8),
                   labels = c('','C','D'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure5.jpg",fig_path), fig5,
       height = 230, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 600)



##### SUPPLEMENTARY FIGURE WITH ESCAPEES
figS6a <- rbind(temp1, temp2) %>%
  filter(escapee == 'Yes') %>%
  ggplot(aes(x = FeEDTA, y = OD)) +
  stat_summary(col = 'black', alpha = 0.9, size = 0.3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  labs(x = 'H<sub>2</sub>S Chelator',
       y = 'OD<sub>600</sub> at Saturation') +
  facet_grid(.~Strain, labeller = labeller(Strain = c('FY4' = 'FY4',
                                                      'met15del' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = ggtext::element_markdown(size = txt),
        # axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                                margin = margin(1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,8.2))


figS6b <- rbind(temp3, temp4) %>%
  filter(escapee == 'Yes') %>%
  ggplot(aes(x = FeEDTA, y = OD)) +
  stat_summary(col = 'black', alpha = 0.9, size = .3,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_x_discrete(labels = c('Absent' = 'Absent\nPreviously',
                              'Present' = 'Present\nPreviously')) +
  labs(x = 'H<sub>2</sub>S Chelator',
       y = 'OD<sub>600</sub> at Saturation') +
  facet_grid(.~ORF, labeller = labeller(ORF = c('FY4' = 'FY4',
                                                'met15D' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
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
  coord_cartesian(ylim = c(0,8.2))


figSH2SFlask <- plot_grid(figS6a, figS6b, ncol = 2, rel_widths = c(1,1.8),
                          align = 'hv',
                  labels = c('A','B'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/FigureSH2SFlask.jpg",fig_path), figSH2SFlask,
       height = one.c, width = one.5c, units = 'mm',
       dpi = 600)

