
strain.labs2 <- data.frame(orf_name = c("FY4","FY4-met3del","FY4-met15del","BY4742","BY4741"), 
                           labels = c('FY4','FY4-italic(met3Δ)','FY4-italic(met15Δ)','BY4742','BY4741'),
                           met_aux = c('Prototroph','Presumed Auxotroph','Presumed Auxotroph','Prototroph','Presumed Auxotroph'),
                           ura_aux = c('Prototroph','Prototroph','Prototroph','Presumed Auxotroph','Presumed Auxotroph'),
                           auxotrophy = c('None','Methionine','Methionine','Uracil','Both'),
                           stringsAsFactors = F)
strain.labs2$auxotrophy <- factor(strain.labs2$auxotrophy, levels = c('None','Methionine','Uracil','Both'))

data.sum <- merge(crbn_src_data %>%
                    filter(methionine == '-Met +Ura') %>%
                    group_by(base, carbon, methionine, cysteine, orf_name) %>%
                    summarize(f = median(relative_fitness, na.rm = T), .groups = 'keep') %>% data.frame(),
                  crbn_src_data %>%
                    filter(methionine == '+Met -Ura') %>%
                    group_by(base, carbon, methionine, cysteine, orf_name) %>%
                    summarize(f = median(relative_fitness, na.rm = T), .groups = 'keep') %>% data.frame(),
                  by = c('base','carbon','cysteine','orf_name'), suffixes = c('_MM','_PM'),
                  all = T)
data.sum <- merge(data.sum, strain.labs2, by = 'orf_name')
data <- merge(crbn_src_data, strain.labs2, by = 'orf_name')

fig1B.2 <- data.sum %>%
  filter(base == 'SD', orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = f_PM, y = f_MM)) +
  # geom_abline(linetype = 'dashed', col = '#757575') +
  geom_hline(yintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_vline(xintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_point(aes(shape = carbon), size = 2.5) +
  geom_label_repel(aes(label = labels, fill = auxotrophy),
                   parse = T, size = 2.5) +
  labs(x = 'Median Relative Fitness in SD+Met-Cys-Ura',
       y = 'Median Relative Fitness in SD-Met-Cys+Ura') +
  scale_fill_manual(name = 'Methionine-\nUracil\nAuxotrophy',
                    values = c('None' = '#FFC107',
                               'Methionine' = '#536DFE',
                               'Uracil' = '#E040FB',
                               'Both' = '#D1C4E9')) +
  scale_shape_manual(name = 'Carbon\nSource',
                     values = c('Glucose' = 16,
                                'Galactose' = 17,
                                'Ethanol' = 15)) +
  scale_x_continuous(breaks = seq(-1,2,0.2)) + 
  scale_y_continuous(breaks = seq(-1,2,0.2)) +
  coord_cartesian(xlim = c(0, 1.2),
                  ylim = c(0, 1.2)) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              face = 'bold',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  guides(fill = guide_legend(nrow=2, byrow=TRUE, order = 2),
         color = guide_legend(nrow=2, byrow=TRUE, order = 3),
         shape = guide_legend(nrow=3, byrow=TRUE, order = 1))
# save(fig1B.2, file = 'figures/final/fig1B.2.RData')


##### FIG 1B.2 with BOX PLOTS
data$carbon2 <- factor(data$carbon, levels = c('Ethanol','Galactose','Glucose'))
bp.pm <- data %>%
  filter(methionine == '+Met -Ura', base == 'SD', orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = relative_fitness, y = orf_name)) +
  geom_vline(xintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(aes(fill = auxotrophy), outlier.shape = NA) +
  # geom_violin() +
  # geom_text_repel(data = gc.res.sum2,
  #                 aes(x = rel_auc_PM, y = orf_name, label = labels),
  #                 parse = T, size = 2.5) +
  scale_x_continuous(breaks = seq(-1,2,0.2)) +
  scale_y_discrete(labels = c( 'FY4' = 'FY4',
                               'FY4_pet' = 'FY4 (&rho;-)' ,
                               'FY4-met15del' = 'FY4-*met15Δ*',
                               'FY4-met15del_pet' = 'FY4-*met15Δ* (&rho;-)',
                               'BY4742' = 'BY4742',
                               'BY4742_pet' = 'BY4742 (&rho;-)',
                               'BY4741' = 'BY4741',
                               'BY4741_pet' = 'BY4741 (&rho;-)'),
                   position = 'right') +
  scale_fill_manual(name = 'Methionine-\nUracil\nAuxotrophy',
                    values = c('None' = '#FFC107',
                               'Methionine' = '#536DFE',
                               'Uracil' = '#E040FB',
                               'Both' = '#D1C4E9'),
                    guide = F) +
  facet_wrap(.~carbon2, strip.position="left", nrow = 3) +
  coord_cartesian(xlim = c(0,1.2)) +
  theme_minimal() +
  theme(plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(size = txt),
        axis.text.y.right = ggtext::element_markdown(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")))

bp.mm <- data %>%
  filter(methionine == '-Met +Ura', base == 'SD', orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = relative_fitness, y = orf_name)) +
  geom_vline(xintercept = 0.299, linetype = 'dashed', col = '#009688') +
  geom_boxplot(aes(fill = auxotrophy), outlier.shape = NA) +
  # geom_violin() +
  # geom_text_repel(data = gc.res.sum2,
  #                 aes(x = rel_auc_PM, y = orf_name, label = labels),
  #                 parse = T, size = 2.5) +
  scale_x_continuous(breaks = seq(-1,2,0.2)) +
  scale_y_discrete(labels = c( 'FY4' = 'FY4',
                               'FY4_pet' = 'FY4 (&rho;-)' ,
                               'FY4-met15del' = 'FY4-*met15Δ*',
                               'FY4-met15del_pet' = 'FY4-*met15Δ* (&rho;-)',
                               'BY4742' = 'BY4742',
                               'BY4742_pet' = 'BY4742 (&rho;-)',
                               'BY4741' = 'BY4741',
                               'BY4741_pet' = 'BY4741 (&rho;-)')) +
  scale_fill_manual(name = 'Methionine-\nUracil\nAuxotrophy',
                    values = c('None' = '#FFC107',
                               'Methionine' = '#536DFE',
                               'Uracil' = '#E040FB',
                               'Both' = '#D1C4E9'),
                    guide = F) +
  facet_wrap(.~carbon) +
  coord_flip(xlim = c(0,1.2)) +
  theme_minimal() +
  theme(plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 270, vjust = 0.5, hjust = 0),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")))
fig1B.2.2 <- cowplot::plot_grid(bp.pm, NULL, fig1B.2, bp.mm,
                            nrow = 2, ncol = 2,
                            align = 'vh', axis = 'tblr',
                            rel_widths = c(2,1.2), rel_heights = c(1.2,2))
# save(fig5B, file = 'figures/final/fig1B.2.RData')
ggsave(sprintf("%s/final/Figure1B.2.jpg",fig_path), fig1B.2.2,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)


