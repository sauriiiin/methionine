

strain.labs2 <- data.frame(orf_name = c('FY4', 'FY4_pet', 'FY4-met15del', 'FY4-met15del_pet',
                                        'BY4742', 'BY4742_pet', 'BY4741', 'BY4741_pet'), 
                           labels = c('FY4', 'FY4~(rho)' ,'FY4-italic(met15Δ)', 'FY4-italic(met15Δ)~(rho)',
                                      'BY4742', 'BY4742~(rho)', 'BY4741', 'BY4741~(rho)'),
                           met_aux = c('Prototroph','Prototroph','Presumed Auxotroph','Presumed Auxotroph',
                                       'Prototroph','Prototroph','Presumed Auxotroph','Presumed Auxotroph'),
                           pet = c('No','Yes','No','Yes','No','Yes','No','Yes'))

gc.res.sum2 <- merge(merge(gc.res.sum[str_detect(gc.res.sum$orf_name, o),] %>%
                             filter(condition != 'SD+Met-Ura+Glu', methionine == '-Met') %>%
                             group_by(condition,carbon,methionine,orf_name) %>%
                             summarize(rel_auc = median(rel_auc, na.rm = T), .groups = 'keep'),
                           gc.res.sum[str_detect(gc.res.sum$orf_name, o),] %>%
                             filter(condition != 'SD+Met-Ura+Glu', methionine == '+Met') %>%
                             group_by(condition,carbon,methionine,orf_name) %>%
                             summarize(rel_auc = median(rel_auc, na.rm = T), .groups = 'keep'),
                           by = c('carbon','orf_name'), suffixes = c('_MM','_PM')),
                     strain.labs2, by = 'orf_name')

gc.res.sum2 %>% filter(orf_name %in% c('FY4-met15del_pet','FY4-met15del'), carbon == 'Glucose')

sp.pmmm <- gc.res.sum2 %>%
  ggplot(aes(x = rel_auc_PM, y = rel_auc_MM)) +
  geom_abline(linetype = 'dashed', col = 'black') +
  # geom_segment(aes(x = 0.5377937, y = 0.2756857-0.03, xend = 1.0072735, yend = 0.2756857-0.03), size = 0.3,
  #              arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm"))) +
  # geom_segment(aes(x = 1.0072735+0.03, y = 0.2756857, xend = 1.0072735+0.03, yend = 0.3859428), size = 0.3,
  #              arrow = arrow(ends = "both", angle = 90, length = unit(.2,"cm"))) +
  # annotate('text', x = 0.5377937 + (1.0072735-0.5377937)/2, y = 0.2756857,
  #          label = sprintf('diff. = %0.2f', 1.0072735-0.5377937),
  #           hjust = 0.5, size = 2) +
  # annotate('text', x = 1.0072735, y = 0.2756857 + (0.3859428-0.2756857)/2,
  #          label = sprintf('diff. = %0.2f', 0.3859428-0.2756857),
  #          hjust = 1, size = 2) +
  geom_point(aes(shape = carbon), size = 3) +
  geom_label_repel(aes(label = labels, fill = met_aux),
                   parse = T, size = 2.5) +
  scale_fill_manual(name = 'Methionine\nAuxotrophy',
                    values = c('Prototroph' = '#4CAF50',
                               'Presumed Auxotroph' = '#F44336'),
                    labels = c('Prototroph' = 'Prototroph',
                               'Presumed Auxotroph' = 'Presumed Auxotroph')) +
  scale_shape_manual(name = 'Carbon\nSource',
                     values = c('Glucose' = 16,
                                'Galactose' = 17,
                                'Ethanol' = 15)) +
  labs(x = 'Median Relative AUC in SD+Met-Cys',
       y = 'Median Relative AUC in SD-Met-Cys') +
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
         shape = guide_legend(nrow=2, byrow=TRUE, order = 1))
# save(fig5B, file = 'figures/final/fig5B.RData')


fig5A <- gc.pred[str_detect(gc.pred$orf_name, o),]  %>%
  filter(expt_rep %in% c("1","2"), condition != 'SD+Met-Ura+Glu') %>% 
  ggplot(aes(x = Time, y = rel_cs)) +
  # geom_smooth(aes(group = pos, col = orf_name), method = 'loess', se = F, alpha = 0.2, lwd = 0.2, linetype = 'dashed') +
  # geom_smooth(aes(col = orf_name), method = 'loess', se = F) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               aes(group = orf_name, fill = orf_name), geom="ribbon", alpha = 0.4) +
  stat_summary(aes(group = orf_name, col = orf_name), fun=mean, geom="line", lwd =0.7) +
  scale_color_manual(name = 'Strain',
                     label = c( 'FY4' = 'FY4',
                                'FY4_pet' = 'FY4 (&rho;-)' ,
                                'FY4-met15del' = 'FY4-*met15Δ*', 
                                'FY4-met15del_pet' = 'FY4-*met15Δ* (&rho;-)', 
                                'BY4742' = 'BY4742', 
                                'BY4742_pet' = 'BY4742 (&rho;-)', 
                                'BY4741' = 'BY4741', 
                                'BY4741_pet' = 'BY4741 (&rho;-)'),
                     values = c('FY4' = '#E64A19',
                                'FY4_pet' = '#FF9800' ,
                                'FY4-met15del' = '#673AB7', 
                                'FY4-met15del_pet' = '#E040FB', 
                                'BY4742' = '#E64A19', 
                                'BY4742_pet' = '#FF9800', 
                                'BY4741' = '#673AB7', 
                                'BY4741_pet' = '#E040FB')) +
  scale_fill_discrete(guide = F) +
  labs(y = 'Relative Colony Size',
       x = 'Time (hours)') +
  facet_wrap(.~base*cysteine*carbon*methionine, nrow = 1) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
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
  guides(col = guide_legend(nrow=1)) +
  coord_cartesian(ylim = c(0,1.2))
save(fig5A, file = 'figures/final/fig5A.RData')

# fig5.2 <- ggpubr::ggarrange(fig5A, fig5B, nrow = 2, heights = c(1,1.5),
#                           labels = c('A','B'), font.label = list(face = 'bold', size = lbls, family = "sans"))
# save(fig5.2, file = 'figures/final/fig5.2.RData')

######
head(gc.res)
mm.dis <- gc.res[str_detect(gc.res$orf_name, o),] %>%
  filter(methionine == '-Met', expt_rep %in% c("1","2"), condition != 'SD+Met-Ura+Glu') %>%
  ggplot(aes(x = rel_auc, y = orf_name)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_text_repel(data = gc.res.sum2,
  #                 aes(x = rel_auc_MM, y = orf_name, label = labels),
  #                 parse = T, size = 2.5, angle = 270) +
  facet_grid(carbon~.) +
  scale_x_continuous(breaks = seq(-1,2,0.2)) +
  scale_y_discrete(labels = c( 'FY4' = 'FY4',
                               'FY4_pet' = 'FY4 (&rho;-)' ,
                               'FY4-met15del' = 'FY4-*met15Δ*',
                               'FY4-met15del_pet' = 'FY4-*met15Δ* (&rho;-)',
                               'BY4742' = 'BY4742',
                               'BY4742_pet' = 'BY4742 (&rho;-)',
                               'BY4741' = 'BY4741',
                               'BY4741_pet' = 'BY4741 (&rho;-)')) +
  facet_grid(.~carbon) +
  coord_flip(xlim = c(0,1.2)) +
  theme_minimal() +
  theme(plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 270),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              face = 'bold',
                                              margin = margin(0.1,0,0.1,0, "mm")))

pm.dis <- gc.res[str_detect(gc.res$orf_name, o),] %>%
  filter(methionine == '+Met', expt_rep %in% c("1","2"), condition != 'SD+Met-Ura+Glu') %>%
  ggplot(aes(x = rel_auc, y = orf_name)) +
  geom_boxplot(outlier.shape = NA) +
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
  facet_wrap(.~carbon, strip.position="left", nrow = 2) +
  coord_cartesian(xlim = c(0,1.2)) +
  theme_minimal() +
  theme(plot.title = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y.right = ggtext::element_markdown(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              face = 'bold',
                                              margin = margin(0.1,0,0.1,0, "mm")))


fig5B <- cowplot::plot_grid(pm.dis, NULL, sp.pmmm, mm.dis,
                   nrow = 2, ncol = 2,
                   align = 'vh', axis = 'tblr',
                   rel_widths = c(2, 1), rel_heights = c(1,2))
save(fig5B, file = 'figures/final/fig5B.RData')
ggsave(sprintf("%s/final/Figure5B.jpg",fig_path), fig5B,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)
