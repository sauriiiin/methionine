
temp <- crbn_src_data[str_detect(crbn_src_data$orf_name, 'BY') &
                !(str_detect(crbn_src_data$condition, 'PlM')) &
                crbn_src_data$base != 'SC',] %>%
  group_by(condition, base, cysteine, carbon, methionine, orf_name) %>%
  summarise(f = median(relative_fitness, na.rm = T), .groups = 'keep') %>% data.frame()

temp.es <- data.frame()
for (b in unique(temp$base)) {
  for (c in unique(temp$carbon[temp$base == b])) {
    for (cys in unique(temp$cysteine[temp$carbon == c & temp$base == b])) {
      for (m in unique(temp$methionine[temp$cysteine == cys & temp$carbon == c &temp$base == b])) {
        es <- (temp$f[temp$methionine == m & temp$cysteine == cys & temp$carbon == c &
                       temp$base == b & temp$orf_name == 'BY4741'] - 
                 temp$f[temp$methionine == m & temp$cysteine == cys & temp$carbon == c &
                          temp$base == b & temp$orf_name == 'BY4742'])/
          temp$f[temp$methionine == m & temp$cysteine == cys & temp$carbon == c &
                   temp$base == b & temp$orf_name == 'BY4742']
        temp.es <- rbind(temp.es, data.frame(methionine = m, cysteine = cys, carbon = c, base = b, effsize = es))
      }
    }
  }
}


fig1b <- crbn_src_data[str_detect(crbn_src_data$orf_name, 'BY') &
                !(str_detect(crbn_src_data$condition, 'PlM')) &
                crbn_src_data$base != 'SC',] %>%
  # filter(condition %in% c('MiM_Glu','MiU_Glu')) %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) + 
  geom_boxplot(aes(fill = as.character(bio_rep)), size = .3, outlier.shape = NA) +
  stat_compare_means(aes(label = ..p.signif..), method = 'kruskal',
                     label.x = 1.4, label.y = 1.4, size = 2.2, col = 'red') +
  geom_text(data = temp.es, aes(x = 1.45, y = 1.35, label = sprintf('%0.2f%%',effsize*100)), size = 2) +
  labs(x = 'Strain', y = 'Relative Colony Size') +
  scale_fill_manual(name = 'Biological Replicate',
                    values = c("1" = "#673AB7",
                               "2" = "#009688",
                               "3" = "#607D8B",
                               "4" = "#FFC107")) +
  facet_grid(base*cysteine*carbon~methionine, 
             labeller = labeller(orf_name = c("FY4"="FY4",
                                              "FY4-met3del"="FY4-*met3Δ*",
                                              "FY4-met15del"="FY4-*met15Δ*",
                                              "BY4742"="BY4742",
                                              "BY4741"="BY4741"))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        # axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        # axis.text.x = element_blank(),
        # axis.text.x = ggtext::element_markdown(angle = 40, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              face = 'bold',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))


temp <- crbn_src_data[str_detect(crbn_src_data$orf_name, 'FY') &
                        crbn_src_data$orf_name != "FY4-met3del" &
                        !(str_detect(crbn_src_data$condition, 'PlM')) &
                        crbn_src_data$uracil == '+Ura' &
                        crbn_src_data$base != 'SC',] %>%
  group_by(condition, base, cysteine, carbon, methionine, orf_name) %>%
  summarise(f = median(relative_fitness, na.rm = T), .groups = 'keep') %>% data.frame()

temp.es <- data.frame()
for (b in unique(temp$base)) {
  for (c in unique(temp$carbon[temp$base == b])) {
    for (cys in unique(temp$cysteine[temp$carbon == c & temp$base == b])) {
        es <- (temp$f[temp$cysteine == cys & temp$carbon == c &
                        temp$base == b & temp$orf_name == 'FY4-met15del'] -
                 temp$f[temp$cysteine == cys & temp$carbon == c &
                          temp$base == b & temp$orf_name == 'FY4'])/
          temp$f[temp$cysteine == cys & temp$carbon == c &
                   temp$base == b & temp$orf_name == 'FY4']
        temp.es <- rbind(temp.es, data.frame(methionine = '-Met +Ura', cysteine = cys, carbon = c, base = b, effsize = es))
    }
  }
}


fig1c <- crbn_src_data[str_detect(crbn_src_data$orf_name, 'FY') &
                         crbn_src_data$orf_name != "FY4-met3del" &
                         !(str_detect(crbn_src_data$condition, 'PlM')) &
                         crbn_src_data$uracil == '+Ura' &
                         crbn_src_data$base != 'SC',] %>%
  # filter(condition %in% c('MiM_Glu','MiU_Glu')) %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) + 
  geom_boxplot(aes(fill = as.character(bio_rep)), size = .3, outlier.shape = NA) +
  stat_compare_means(aes(label = ..p.signif..), method = 'kruskal',
                     label.x = 1.4, label.y = 1.4, size = 2.2, col = 'red') +
  geom_text(data = temp.es, aes(x = 1.45, y = 1.35, label = sprintf('%0.2f%%',effsize*100)), size = 2) +
  labs(x = 'Strain', y = 'Relative Colony Size') +
  scale_fill_manual(name = 'Biological Replicate',
                    values = c("1" = "#673AB7",
                               "2" = "#009688",
                               "3" = "#607D8B",
                               "4" = "#FFC107")) +
  scale_x_discrete(labels = c("FY4"="FY4",
                              "FY4-met3del"="FY4-*met3Δ*",
                              "FY4-met15del"="FY4-*met15Δ*",
                              "BY4742"="BY4742",
                              "BY4741"="BY4741")) +
  facet_grid(base*cysteine*carbon~methionine) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        # axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        # axis.text.x = element_blank(),
        axis.text.x = ggtext::element_markdown(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              face = 'bold',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))


temp <- crbn_src_data[str_detect(crbn_src_data$orf_name, 'FY') &
                        crbn_src_data$orf_name != "FY4-met3del" &
                        !(str_detect(crbn_src_data$condition, 'PlM')) &
                        crbn_src_data$uracil == '+Ura' &
                        crbn_src_data$base == 'SC',] %>%
  group_by(condition, base, cysteine, carbon, methionine, orf_name) %>%
  summarise(f = median(relative_fitness, na.rm = T), .groups = 'keep') %>% data.frame()

temp.es <- data.frame()
for (b in unique(temp$base)) {
  for (c in unique(temp$carbon[temp$base == b])) {
    for (cys in unique(temp$cysteine[temp$carbon == c & temp$base == b])) {
      es <- (temp$f[temp$cysteine == cys & temp$carbon == c &
                      temp$base == b & temp$orf_name == 'FY4-met15del'] -
               temp$f[temp$cysteine == cys & temp$carbon == c &
                        temp$base == b & temp$orf_name == 'FY4'])/
        temp$f[temp$cysteine == cys & temp$carbon == c &
                 temp$base == b & temp$orf_name == 'FY4']
      temp.es <- rbind(temp.es, data.frame(methionine = '-Met +Ura', cysteine = cys, carbon = c, base = b, effsize = es))
    }
  }
}


fig1d.a <- crbn_src_data[str_detect(crbn_src_data$orf_name, 'FY') &
                         crbn_src_data$orf_name != "FY4-met3del" &
                         !(str_detect(crbn_src_data$condition, 'PlM')) &
                         crbn_src_data$uracil == '+Ura' &
                         crbn_src_data$base == 'SC',] %>%
  # filter(condition %in% c('MiM_Glu','MiU_Glu')) %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) + 
  geom_boxplot(aes(fill = as.character(bio_rep)), size = .3, outlier.shape = NA) +
  stat_compare_means(aes(label = ..p.signif..), method = 'kruskal',
                     label.x = 1.4, label.y = 1.4, size = 2.2, col = 'red') +
  geom_text(data = temp.es, aes(x = 1.45, y = 1.35, label = sprintf('%0.2f%%',effsize*100)), size = 2) +
  labs(x = 'Strain', y = 'Relative Colony Size') +
  scale_fill_manual(name = 'Biological Replicate',
                    values = c("1" = "#673AB7",
                               "2" = "#009688",
                               "3" = "#607D8B",
                               "4" = "#FFC107")) +
  scale_x_discrete(labels = c("FY4"="FY4",
                              "FY4-met3del"="FY4-*met3Δ*",
                              "FY4-met15del"="FY4-*met15Δ*",
                              "BY4742"="BY4742",
                              "BY4741"="BY4741")) +
  facet_grid(base*cysteine*carbon~methionine) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        # axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        # axis.text.x = element_blank(),
        axis.text.x = ggtext::element_markdown(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              face = 'bold',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))


temp <- crbn_src_data[str_detect(crbn_src_data$orf_name, 'BY') &
                        crbn_src_data$orf_name != "FY4-met3del" &
                        !(str_detect(crbn_src_data$condition, 'PlM')) &
                        crbn_src_data$uracil == '+Ura' &
                        crbn_src_data$base == 'SC',] %>%
  group_by(condition, base, cysteine, carbon, methionine, orf_name) %>%
  summarise(f = median(relative_fitness, na.rm = T), .groups = 'keep') %>% data.frame()

temp.es <- data.frame()
for (b in unique(temp$base)) {
  for (c in unique(temp$carbon[temp$base == b])) {
    for (cys in unique(temp$cysteine[temp$carbon == c & temp$base == b])) {
      es <- (temp$f[temp$cysteine == cys & temp$carbon == c &
                      temp$base == b & temp$orf_name == 'BY4741'] -
               temp$f[temp$cysteine == cys & temp$carbon == c &
                        temp$base == b & temp$orf_name == 'BY4742'])/
        temp$f[temp$cysteine == cys & temp$carbon == c &
                 temp$base == b & temp$orf_name == 'BY4742']
      temp.es <- rbind(temp.es, data.frame(methionine = '-Met +Ura', cysteine = cys, carbon = c, base = b, effsize = es))
    }
  }
}


fig1d.b <- crbn_src_data[str_detect(crbn_src_data$orf_name, 'BY') &
                           crbn_src_data$orf_name != "FY4-met3del" &
                           !(str_detect(crbn_src_data$condition, 'PlM')) &
                           crbn_src_data$uracil == '+Ura' &
                           crbn_src_data$base == 'SC',] %>%
  # filter(condition %in% c('MiM_Glu','MiU_Glu')) %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) + 
  geom_boxplot(aes(fill = as.character(bio_rep)), size = .3, outlier.shape = NA) +
  stat_compare_means(aes(label = ..p.signif..), method = 'kruskal',
                     label.x = 1.4, label.y = 1.4, size = 2.2, col = 'red') +
  geom_text(data = temp.es, aes(x = 1.45, y = 1.35, label = sprintf('%0.2f%%',effsize*100)), size = 2) +
  labs(x = 'Strain', y = 'Relative Colony Size') +
  scale_fill_manual(name = 'Biological Replicate',
                    values = c("1" = "#673AB7",
                               "2" = "#009688",
                               "3" = "#607D8B",
                               "4" = "#FFC107")) +
  scale_x_discrete(labels = c("FY4"="FY4",
                              "FY4-met3del"="FY4-*met3Δ*",
                              "FY4-met15del"="FY4-*met15Δ*",
                              "BY4742"="BY4742",
                              "BY4741"="BY4741")) +
  facet_grid(base*cysteine*carbon~methionine) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        # axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        # axis.text.x = element_blank(),
        axis.text.x = ggtext::element_markdown(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              face = 'bold',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))
fig1d <- ggpubr::ggarrange(fig1d.b, fig1d.a, nrow = 1, common.legend = T, legend = 'bottom')

# fig1e <- crbn_src_data[str_detect(crbn_src_data$condition, 'PlM'),] %>%
#   # filter(condition %in% c('MiM_Glu','MiU_Glu')) %>%
#   ggplot(aes(x = orf_name, y = relative_fitness)) +
#   geom_boxplot(size = .3, outlier.shape = NA, fill = '#9E9E9E') +
#   # stat_compare_means(aes(label = ..p.signif..), method = 'kruskal',
#   #                    label.x = 1.4, label.y = 1.45, size = 2.2, col = 'red') +
#   labs(x = 'Strain', y = 'Relative Colony Size') +
#   scale_fill_manual(name = 'Biological Replicate',
#                     values = c("1" = "#673AB7",
#                                "2" = "#009688",
#                                "3" = "#607D8B",
#                                "4" = "#FFC107")) +
#   scale_x_discrete(labels = c("FY4"="FY4",
#                               "FY4-met3del"="FY4-*met3Δ*",
#                               "FY4-met15del"="FY4-*met15Δ*",
#                               "BY4742"="BY4742",
#                               "BY4741"="BY4741")) +
#   facet_grid(base*cysteine*carbon~methionine) +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         # axis.title.x = element_blank(),
#         axis.text = element_text(size = txt),
#         # axis.text.x = element_blank(),
#         axis.text.x = ggtext::element_markdown(),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = ggtext::element_markdown(size = txt,
#                                               face = 'bold',
#                                               margin = margin(0.1,0,0.1,0, "mm"))) +
#   coord_cartesian(ylim = c(0,1.5))

fig1 <- ggpubr::ggarrange(fig1A, ggpubr::ggarrange(ggpubr::ggarrange(fig1b, fig1c, nrow = 1, widths = c(1.6,1),
                                                                     labels = c('B','C'), common.legend = T, legend = 'none',
                                                                     font.label = list(face = 'bold', size = lbls, family = "sans")),
                                                   fig1d, nrow = 2, heights = c(3,2),
                                                   labels = c('','D'), common.legend = T, legend = 'bottom',font.label = list(face = 'bold', size = lbls, family = "sans")),
                          labels = c('A',''), nrow = 2, heights = c(1,2), font.label = list(face = 'bold', size = lbls, family = "sans"))
ggsave(sprintf("%s/final/Figure1.jpg",fig_path), fig1,
       height = 270, width = two.c, units = 'mm',
       dpi = 600)

# ggpubr::ggarrange(fig1d, fig1e, nrow = 1, widths = c(2,1),
