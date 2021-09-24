
data.mdl <- read.csv('/home/sbp29/R/Projects/methionine/data/modeling/YLL_simulation.csv', stringsAsFactors = F)

log10_rev_trans <- trans_new(
  "log10_rev",
  function(x) log10(rev(x)),
  function(x) rev(10 ^ (x)),
  log_breaks(10),
  domain = c(1e-100, Inf)
)

data.mdl$label[data.mdl$Model == 'A'] <- 'A. Default'
data.mdl$label[data.mdl$Model == 'B'] <- 'B. A - All YLL058W reactions'
data.mdl$label[data.mdl$Model == 'C'] <- 'C. B + Hypothesized YLL058W reaction'
data.mdl$label[data.mdl$Model == 'D'] <- 'D. C - All MET15 reactions'

data.mdl %>%
  ggplot(aes(x = Yll058w.Met15, y = Growth)) +
  geom_line(aes(col = Model)) +
  scale_x_continuous(trans = 'log10') +
  scale_color_manual(name = 'Model',
                     breaks = c('A','B','C','D'),
                     lables = c('A. Default',
                                'B. A - All YLL058W reactions',
                                'C. B + Hypothesized YLL058W reaction',
                                'D. C - All MET15 reactions'),
                     values = c('A' = '#212121',
                                'B' = '#5D4037',
                                'C' = '#FF5722',
                                'D' = '#FFC107')) +
  # facet_zoom(ylim = c(0.3225,0.3275), zoom.size = .5) +
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
  guides(fill = guide_legend(override.aes=list(shape = 15, alpha = 1)))
