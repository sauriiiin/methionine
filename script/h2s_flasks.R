library(pzfx)
library(growthcurver)

strain.labs.h2s <- read.csv(file = 'data/h2s/FeEDTA_samples.csv')
data.h2s <- read.csv(file = 'data/h2s/FeEDTA_readings.csv')

data.h2s %>%
  filter(Time != 2655) %>%
  melt(id.vars = 'Time', variable.name = 'Flask', value.name = 'OD') %>%
  ggplot(aes(x = Time, y = rollOD, col = Flask)) +
  geom_line()

# data.h2s.smooth <- data.h2s %>%
#   filter(Time != 2655) %>%
#   melt(id.vars = 'Time', variable.name = 'Flask', value.name = 'OD') %>%
#   group_by(Flask) %>%
#   mutate(rollOD = zoo::rollmean(OD, k = 3, fill = 'extend', align = 'center')) %>% 
#   data.frame()
# 
# data.h2s.smooth %>%
#   ggplot(aes(x = Time, y = rollOD, col = Flask)) +
#   geom_line()

data.h2s <- data.h2s %>%
  filter(Time != 2655)

data.pred.h2s <- NULL
t <- data.h2s$Time
for (c in colnames(data.h2s)[2:dim(data.h2s)[2]]) {
  temp <- data.h2s[c]
  temp[temp <= 0] <- temp[temp <= 0] + 0.0001
  lo <- loess.smooth(t, log(temp),
                     span = 0.65, evaluation = 50, degree = 2,
                     family = 'gaussian')
  data.pred.h2s <- cbind(data.pred.h2s,exp(lo$y))
}
data.pred.h2s <- cbind(lo$x, data.pred.h2s)
colnames(data.pred.h2s) <- colnames(data.h2s)
data.pred.h2s <- data.frame(data.pred.h2s)
head(data.pred.h2s)

data.h2s.gc.res <- SummarizeGrowthByPlate(data.pred.h2s)

fig6a <- merge(data.pred.h2s %>%
  melt(id.vars = 'Time', variable.name = 'Flask', value.name = 'OD'), strain.labs.h2s, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Expected.Initial.OD == 0.1) %>%
  ggplot(aes(x = Time, y = OD, col = Strain)) +
  geom_line(aes(linetype = Replicate), size = 0.8) +
  scale_y_continuous(breaks = seq(-2,10,1)) +
  labs(x = 'Time (minutes)',
       y = 'Optical Density (OD) at 600nm') +
  scale_color_manual(name = 'Strain',
                     values = c('FY4' = '#FFC107',
                                'met15del' = '#536DFE'),
                     labels = c('FY4' = 'FY4',
                                'met15del' = 'FY4-*met15Δ*')) +
  scale_linetype_manual(name = 'Biological Replicate',
                        values = c('Rep_1' = 1,
                                   'Rep_2' = 3,
                                   'Rep_3' = 4),
                        labels = c('Rep_1' = '#1',
                                   'Rep_2' = '#2',
                                   'Rep_3' = '#3')) +
  facet_wrap(.~FeEDTA, labeller = labeller(FeEDTA = c('Absent' = 'SD - Met + Glu<br />w/o H<sub>2</sub>S Chelator',
                                                      'Present' = 'SD - Met + Glu<br />w/ H<sub>2</sub>S Chelator'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        # legend.direction = 'vertical',
        # legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(color = guide_legend(nrow=1, byrow=TRUE, order = 1,
                              override.aes = list(size = 3)),
         linetype = guide_legend(nrow=1, byrow=TRUE, order = 2,
                              override.aes = list(size = 0.5)))

fig6b <- merge(data.pred.h2s %>%
        melt(id.vars = 'Time', variable.name = 'Flask', value.name = 'OD'), strain.labs.h2s, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Time == 3075, Expected.Initial.OD == 0.1)  %>%
  ggplot(aes(x = Replicate, y = OD)) +
  stat_summary(col = 'black', alpha = 0.9, size = 1,
               fun = mean, geom = "bar") +
  scale_x_discrete(labels = c('Rep_1' = '#1',
                              'Rep_2' = '#2',
                              'Rep_3' = '#3')) +
  labs(x = 'Biological Replicate',
       y = 'Optical Density (OD) at 600nm') +
  facet_grid(Strain~FeEDTA, labeller = labeller(FeEDTA = c('Absent' = 'SD - Met + Glu<br />w/o H<sub>2</sub>S Chelator',
                                                           'Present' = 'SD - Met + Glu<br />w/ H<sub>2</sub>S Chelator'),
                                                Strain = c('FY4' = 'FY4',
                                                           'met15del' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        # legend.direction = 'vertical',
        # legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                                margin = margin(0.1,0,0.1,0, "mm")),
        strip.text.y = ggtext::element_markdown(size = txt,
                                                margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,8))


#####
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
colnames(temp) <- c('Strain','FeEDTA','Replicate','Reading')
temp$FeEDTA <- as.character(temp$FeEDTA)
temp$FeEDTA[temp$FeEDTA == 'Fe'] <- 'Present'
temp$FeEDTA[temp$FeEDTA == ''] <- 'Absent'

data.h2s.re <- cbind(data.h2s.re, temp)

fig6c <- data.h2s.re[data.h2s.re$Time..Minutes. == 3075,-2] %>%
  group_by(FeEDTA, Strain, Replicate) %>%
  summarize(OD = median(OD, na.rm = T), .groups = 'keep') %>%
  ggplot(aes(x = Replicate, y = OD)) +
  stat_summary(col = 'black', alpha = 0.9, size = 1,
               fun = mean, geom = "bar") +
  scale_x_discrete(labels = c('R1' = '#1',
                              'R2' = '#2',
                              'R3' = '#3')) +
  labs(x = 'Biological Replicate',
       y = 'Optical Density (OD) at 600nm') +
  facet_grid(Strain~FeEDTA, labeller = labeller(FeEDTA = c('Absent' = 'SD - Met + Glu (H<sub>2</sub>S<br />Chelator Previously Absent)',
                                                             'Present' = 'SD - Met + Glu (H<sub>2</sub>S<br />Chelator Previously Present)'),
                                                  Strain = c('FY4' = 'FY4',
                                                             'met15D' = 'FY4-*met15Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        # legend.direction = 'vertical',
        # legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text.x = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")),
        strip.text.y = ggtext::element_markdown(size = txt,
                                                margin = margin(0.1,0,0.1,0, "mm")))  +
  coord_cartesian(ylim = c(0,8))


fig6 <- cowplot::plot_grid(fig6a, cowplot::plot_grid(fig6b, fig6c, nrow = 1, rel_widths = c(1,1),
                                             labels = c('B','C'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                   ncol = 1, rel_heights = c(1,1),
                   labels = c('A',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure6.jpg",fig_path), fig6,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)

