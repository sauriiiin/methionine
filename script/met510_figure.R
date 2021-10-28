

data.510$orf_name <- factor(data.510$orf_name, levels = c('FY4','BY4741_met10del','BY4741_met5del',
                                                          'met15del','met10del','met5del'))

merge(data.510, strain.labs.510, by = 'orf_name') %>%
  filter(hours == 160) %>%
  ggplot(aes(x = condition, y = relative_fitness)) +
  geom_boxplot(fill = '#616161', outlier.shape = NA, size = 0.2,
               position = position_dodge(width = 0.85)) +
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  scale_x_discrete(limits = c('YPDA', 
                              'SD+Met',
                              'SD-Met',
                              'SD-Ura')) +
  labs(x = 'Strain',
       y = 'Relative Colony Size') +
  coord_cartesian(ylim = c(0,1.3)) +
  facet_wrap(~orf_name, ncol = 3, 
             labeller = labeller(orf_name = c('FY4'='FY4',
                                              'met15del'='FY4-*met15Δ*',
                                              'met5del'='FY4-*met5Δ*',
                                              'met10del'='FY4-*met10Δ*',
                                              'BY4741_met5del'='BY4741-*met5Δ*',
                                              'BY4741_met10del'='BY4741-*met10Δ*'))) +
  theme_linedraw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = titles),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 30, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.direction = 'vertical',
        legend.box = 'vertical',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(0.1,0,0.1,0, "mm")))
