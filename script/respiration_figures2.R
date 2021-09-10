
source("~/R/Projects/methionine/functions/colorstrip.R")

gc.res$expt_rep <- as.character(gc.res$expt_rep)
head(gc.res.sum)

##### STRAIN LABELES
unique(gc.res.sum$orf_name)
strain.labs <- c('FY4', 'FY4 (&rho;-)' ,'*met15Δ*', '*met15Δ* (&rho;-)', 'BY4742', 'BY4742 (&rho;-)', 'BY4741', 'BY4741 (&rho;-)')

for (c in unique(gc.res.sum$condition)) {
  for (e in unique(gc.res.sum$expt_rep[gc.res.sum$condition == c])) {
    gc.res.sum$rel_auc[gc.res.sum$expt_rep == e & gc.res.sum$condition == c] <-
      gc.res.sum$auc_e[gc.res.sum$expt_rep == e & gc.res.sum$condition == c]/median(gc.res.sum$auc_e[gc.res.sum$expt_rep == e & 
                                                                                                       gc.res.sum$condition == c &
                                                                                                       gc.res.sum$orf_name == 'FY4'])
    gc.res.sum$rel_gr[gc.res.sum$expt_rep == e & gc.res.sum$condition == c] <-
      gc.res.sum$gr[gc.res.sum$expt_rep == e & gc.res.sum$condition == c]/median(gc.res.sum$gr[gc.res.sum$expt_rep == e & 
                                                                                                 gc.res.sum$condition == c &
                                                                                                 gc.res.sum$orf_name == 'FY4'])
    gc.res$rel_auc[gc.res$expt_rep == e & gc.res$condition == c] <-
      gc.res$auc_e[gc.res$expt_rep == e & gc.res$condition == c]/median(gc.res$auc_e[gc.res$expt_rep == e & 
                                                                                                       gc.res$condition == c &
                                                                                                       gc.res$orf_name == 'FY4'])
    gc.res$rel_gr[gc.res$expt_rep == e & gc.res$condition == c] <-
      gc.res$gr[gc.res$expt_rep == e & gc.res$condition == c]/median(gc.res$gr[gc.res$expt_rep == e & 
                                                                                                 gc.res$condition == c &
                                                                                                 gc.res$orf_name == 'FY4'])
  }
  for (o in unique(gc.res.sum$orf_name[gc.res.sum$condition == c & gc.res.sum$orf_name != 'FY4'])) {
    temp.kw <- compare_means(data = gc.res.sum[gc.res.sum$condition == c & gc.res.sum$orf_name %in% c(o, 'FY4'),],
                             method = 'kruskal', formula = rel_auc ~ orf_name) %>% data.frame()
    rauc1 <- median(gc.res.sum$rel_auc[gc.res.sum$condition == c & gc.res.sum$orf_name == o], na.rm = T)
    rauc2 <- median(gc.res.sum$rel_auc[gc.res.sum$condition == c & gc.res.sum$orf_name == 'FY4'], na.rm = T)
    gc.res.sum$kw_p[gc.res.sum$condition == c & gc.res.sum$orf_name == o] <- temp.kw$p
    gc.res.sum$kw_signif[gc.res.sum$condition == c & gc.res.sum$orf_name == o] <- temp.kw$p.signif
    gc.res.sum$emp_effsize[gc.res.sum$condition == c & gc.res.sum$orf_name == o] <- (rauc1 - rauc2)/rauc2
    
    temp.gr.kw <- compare_means(data = gc.res.sum[gc.res.sum$condition == c & gc.res.sum$orf_name %in% c(o, 'FY4'),],
                             method = 'kruskal', formula = gr ~ orf_name) %>% data.frame()
    gr1 <- median(gc.res.sum$gr[gc.res.sum$condition == c & gc.res.sum$orf_name == o], na.rm = T)
    gr2 <- median(gc.res.sum$gr[gc.res.sum$condition == c & gc.res.sum$orf_name == 'FY4'], na.rm = T)
    gc.res.sum$gr_kw_p[gc.res.sum$condition == c & gc.res.sum$orf_name == o] <- temp.gr.kw$p
    gc.res.sum$gr_kw_signif[gc.res.sum$condition == c & gc.res.sum$orf_name == o] <- temp.gr.kw$p.signif
    gc.res.sum$gr_emp_effsize[gc.res.sum$condition == c & gc.res.sum$orf_name == o] <- (gr1 - gr2)/gr2
  }
  gc.res.sum$kw_p[gc.res.sum$condition == c] <- p.adjust(gc.res.sum$kw_p[gc.res.sum$condition == c], method = 'BH')
  gc.res.sum$gr_kw_p[gc.res.sum$condition == c] <- p.adjust(gc.res.sum$gr_kw_p[gc.res.sum$condition == c], method = 'BH')
}

plot_labels <- gc.res.sum %>%
  filter(orf_name != 'FY4') %>%
  group_by(condition, orf_name) %>%
  summarise(effect_size = median(emp_effsize, na.rm = T), signif = max(kw_signif),
            gr_effect_size = median(gr_emp_effsize, na.rm = T), gr_signif = max(gr_kw_signif)) %>%
  data.frame()


##### GROWTH CURVES
gc.pred <- melt(data.pred, id.vars = 'Time', variable.name = 'sample', value.name = 'cs')
head(gc.pred)
temp <- str_split(gc.pred$sample, ',', simplify = T)
colnames(temp) <- c('expt_rep','orf_name','condition','bio_rep','pos')
gc.pred <- cbind(temp, gc.pred)
head(gc.pred)

for (c in unique(gc.pred$condition)) {
  for (e in unique((gc.pred$expt_rep[gc.pred$condition == c]))) {
    max.t <- max(gc.pred$Time[gc.pred$condition == c & gc.pred$expt_rep == e])
    med.cs <- mean(gc.pred$cs[gc.pred$condition == c & gc.pred$expt_rep == e & 
                               gc.pred$Time == max.t & gc.pred$orf_name == 'FY4'], na.rm = T)
    gc.pred$rel_cs[gc.pred$condition == c & gc.pred$expt_rep == e] <- gc.pred$cs[gc.pred$condition == c & gc.pred$expt_rep == e]/med.cs
  }
}
gc.pred$condition <- factor(gc.pred$condition, levels = condition.levels)
gc.pred$orf_name <- factor(gc.pred$orf_name, levels = strain.levels)


#####
conds <- NULL
conds$condition <- unique(gc.res$condition)
conds$base[str_detect(conds$condition, 'SD')] <- 'SD'
conds$methionine[str_detect(conds$condition, '-Met')] <- '-Met'
conds$methionine[is.na(conds$methionine)] <- '+Met'
conds$carbon[str_detect(conds$condition, 'Glu')] <- 'Glucose'
conds$carbon[is.na(conds$carbon)] <- 'Ethanol'
conds$cysteine <- '-Cys'
conds <- data.frame(conds)

gc.res <- merge(gc.res, conds, by = 'condition')
gc.pred <- merge(gc.pred, conds, by = 'condition')
gc.res.sum <- merge(gc.res.sum, conds, by = 'condition')

gc.res$methionine <- factor(gc.res$methionine, levels = c('+Met','-Met'))
gc.pred$methionine <- factor(gc.pred$methionine, levels = c('+Met','-Met'))
gc.pred$carbon <- factor(gc.pred$carbon, levels = c('Glucose','Ethanol'))
gc.res.sum$carbon <- factor(gc.res.sum$carbon, levels = c('Glucose','Ethanol'))

##### FIGURES
for (o in c('FY','BY')) {
  temp <- gc.res[str_detect(gc.res$orf_name, o),] %>%
    filter(condition != 'SD+Met-Ura+Glu') %>% 
    group_by(condition, base, methionine, carbon, cysteine, orf_name) %>%
    summarize(rel_auc = median(rel_auc, na.rm = T), .groups = 'keep')
  temp.es <- data.frame()
  for (b in unique(temp$base)) {
    for (c in unique(temp$carbon[temp$base == b])) {
      for (cys in unique(temp$cysteine[temp$carbon == c & temp$base == b])) {
        for (o1 in unique(temp$orf_name[temp$cysteine == cys & temp$carbon == c &temp$base == b])) {
          es <- (temp$rel_auc[temp$orf_name == o1 & temp$cysteine == cys & temp$carbon == c &
                          temp$base == b & temp$methionine == '-Met'] - 
                   temp$rel_auc[temp$orf_name == o1 & temp$cysteine == cys & temp$carbon == c &
                            temp$base == b & temp$methionine == '+Met'])/
            temp$rel_auc[temp$orf_name == o1 & temp$cysteine == cys & temp$carbon == c &
                     temp$base == b & temp$methionine == '+Met']
          temp.es <- rbind(temp.es, data.frame(orf_name = o1, cysteine = cys, carbon = c, base = b, effsize = es))
        }
      }
    }
  }
  
  plot.rauc.box.kw <- gc.res[str_detect(gc.res$orf_name, o),] %>%
    filter(condition != 'SD+Met-Ura+Glu') %>%
    ggplot(aes(x = methionine, y = rel_auc)) +
    geom_boxplot(aes(fill = bio_rep), outlier.shape = NA) +
    # geom_jitter(aes(col = bio_rep), size = 1) +
    # geom_text(data = plot_labels[str_detect(plot_labels$orf_name, o),] %>% filter(condition != 'SD+Met-Ura+Glu'),
    #           aes(x = orf_name, y = 1.4, label = sprintf('%0.2f%%',effect_size*100)), size = 2.2) +
    # geom_text(data = plot_labels[str_detect(plot_labels$orf_name, o),]  %>% filter(condition != 'SD+Met-Ura+Glu'),
    #           aes(x = orf_name, y = 1.5, label = signif), col = 'red', size = 2.2) +
    stat_compare_means(aes(label = ..p.signif..), method = 'kruskal',
                       label.x = 1.4, label.y = 1.4, size = 2.2, col = 'red') +
    geom_text(data = temp.es, aes(x = 1.45, y = 1.35, label = sprintf('%0.2f%%',effsize*100)), size = 2) +
    scale_fill_manual(name = 'Biological Replicate',
                      values = c("1" = "#673AB7",
                                 "2" = "#009688",
                                 "3" = "#607D8B",
                                 "4" = "#FFC107")) +
    facet_grid(base*cysteine*carbon~orf_name,
               labeller = labeller(orf_name = c( 'FY4' = 'FY4',
                                                 'FY4_pet' = 'FY4 (&rho;-)' ,
                                                 'FY4-met15del' = 'FY4-*met15Δ*',
                                                 'FY4-met15del_pet' = 'FY4-*met15Δ* (&rho;-)',
                                                 'BY4742' = 'BY4742',
                                                 'BY4742_pet' = 'BY4742 (&rho;-)',
                                                 'BY4741' = 'BY4741',
                                                 'BY4741_pet' = 'BY4741 (&rho;-)'))) +
    labs(x = 'Media Condition',
         y = 'Relative Area Under the Curve') +
    # scale_x_discrete(labels = c( 'FY4' = 'FY4',
    #                             'FY4_pet' = 'FY4 (&rho;-)' ,
    #                             'FY4-met15del' = 'FY4-*met15Δ*', 
    #                             'FY4-met15del_pet' = 'FY4-*met15Δ* (&rho;-)', 
    #                             'BY4742' = 'BY4742', 
    #                             'BY4742_pet' = 'BY4742 (&rho;-)', 
    #                             'BY4741' = 'BY4741', 
    #                             'BY4741_pet' = 'BY4741 (&rho;-)')) +
    scale_y_continuous(minor_breaks = seq(-2,2,0.1)) +
    theme_linedraw() +
    theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
          axis.title = element_text(size = titles),
          axis.text = element_text(size = txt),
          # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          # axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
          legend.title = element_text(size = titles),
          legend.text = ggtext::element_markdown(size = txt+2),
          legend.position = 'bottom',
          legend.key.size = unit(3, "mm"),
          legend.box.spacing = unit(0.5,"mm"),
          strip.text = ggtext::element_markdown(size = txt,
                                    face = 'bold',
                                    margin = margin(0.1,0,0.1,0, "mm"))) +
    coord_cartesian(ylim = c(0,1.5))
  # plot.rauc.box.kw <- colorstrip(plot.rauc.box.kw,c("#303F9F","#448AFF","#388E3C","#8BC34A"))
  # ggsave(sprintf("%s/%s/RELATIVE_AUC_BOX_KW_%s.jpg",fig_path, expt.name,o),
  #        plot.rauc.box.kw,
  #        height = 70, width = two.c, units = 'mm',
  #        dpi = 600)
  # 
  # plot.gr.box.kw <- gc.res[str_detect(gc.res$orf_name, o),] %>%
  #   filter(condition != 'SD+Met-Ura+Glu') %>%
  #   ggplot(aes(x = orf_name, y = gr)) +
  #   geom_boxplot(aes(fill = bio_rep), outlier.shape = NA) +
  #   # geom_jitter(aes(col = bio_rep), size = 1) +
  #   geom_text(data = plot_labels[str_detect(plot_labels$orf_name, o),] %>% filter(condition != 'SD+Met-Ura+Glu'),
  #             aes(x = orf_name, y = 0.06, label = sprintf('%0.2f%%',gr_effect_size*100)), size = 2.2) +
  #   geom_text(data = plot_labels[str_detect(plot_labels$orf_name, o),]  %>% filter(condition != 'SD+Met-Ura+Glu'),
  #             aes(x = orf_name, y = 0.065, label = gr_signif), col = 'red', size = 2.2) +
  #   scale_fill_manual(name = 'Biological Replicate',
  #                     values = c("1" = "#673AB7",
  #                                "2" = "#009688",
  #                                "3" = "#607D8B",
  #                                "4" = "#FFC107")) +
  #   facet_wrap(.~condition, nrow = 1) +
  #   labs(x = 'Strain',
  #        y = 'Growth Rate') +
  #   scale_x_discrete(labels = c( 'FY4' = 'FY4',
  #                                'FY4_pet' = 'FY4 (&rho;-)' ,
  #                                'FY4-met15del' = 'FY4-*met15Δ*', 
  #                                'FY4-met15del_pet' = 'FY4-*met15Δ* (&rho;-)', 
  #                                'BY4742' = 'BY4742', 
  #                                'BY4742_pet' = 'BY4742 (&rho;-)', 
  #                                'BY4741' = 'BY4741', 
  #                                'BY4741_pet' = 'BY4741 (&rho;-)')) +
  #   scale_y_continuous(minor_breaks = seq(-2,2,0.005)) +
  #   theme_linedraw() +
  #   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
  #         axis.title = element_text(size = titles),
  #         axis.text = element_text(size = txt+2),
  #         # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  #         axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
  #         legend.title = element_text(size = titles),
  #         legend.text = ggtext::element_markdown(size = txt+2),
  #         legend.position = 'bottom',
  #         legend.key.size = unit(3, "mm"),
  #         legend.box.spacing = unit(0.5,"mm"),
  #         strip.text = element_text(size = txt+2,
  #                                   face = 'bold',
  #                                   margin = margin(0.1,0,0.1,0, "mm"))) +
  #   coord_cartesian(ylim = c(0,0.07))
  # plot.gr.box.kw <- colorstrip(plot.gr.box.kw,c("#303F9F","#448AFF","#388E3C","#8BC34A"))
  # ggsave(sprintf("%s/%s/GROWTHRATE_BOX_KW_%s.jpg",fig_path, expt.name,o),
  #        plot.gr.box.kw,
  #        height = 70, width = two.c, units = 'mm',
  #        dpi = 600)
  # 
  plot.all.gc <- gc.pred[str_detect(gc.pred$orf_name, o),]  %>%
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
    facet_grid(.~base*cysteine*carbon*methionine) +
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
    coord_cartesian(ylim = c(0,1.2))
  # plot.all.gc <- colorstrip(plot.all.gc,c("#303F9F","#448AFF","#388E3C","#8BC34A"))
  # ggsave(sprintf("%s/%s/GROWTH_CURVES_%s.jpg",fig_path, expt.name,o),
  #        plot.all.gc,
  #        height = 70, width = two.c, units = 'mm',
  #        dpi = 600)
  # 
  # plot.growth <- cowplot::plot_grid(plot.all.gc, plot.rauc.box.kw,
  #                             ncol = 1, align = 'v', rel_heights = c(0.9,1),
  #                             labels = c('A','B'),
  #                             label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
  
  if (o == 'FY') {
    fig5 <- ggpubr::ggarrange(plot.all.gc, plot.rauc.box.kw, nrow = 2, heights = c(0.9,1),
                              labels = c('A','B'), font.label = list(face = 'bold', size = lbls, family = "sans"))
    save(fig5, file = 'figures/final/fig5.RData')
    # ggsave(sprintf("%s/%s/GROWTH_AUC_CURVES_%s.jpg",fig_path, expt.name, o),
    #        plot.growth,
    #        height = 140, width = two.c, units = 'mm',
    #        dpi = 600)
  } else {
    figSY <- plot.growth
    save(figSY, file = 'figures/final/figSY.RData')
    # ggsave(sprintf("%s/%s/GROWTH_AUC_CURVES_%s.jpg",fig_path, expt.name, o),
    #        plot.growth,
    #        height = 140, width = two.c, units = 'mm',
    #        dpi = 600)
  }
  # plot.growth2 <- cowplot::plot_grid(plot.all.gc, plot.gr.box.kw,
  #                                    ncol = 1, align = 'v', rel_heights = c(0.9,1),
  #                                    labels = c('A','B'),
  #                                    label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
  # ggsave(sprintf("%s/%s/GROWTH_GR_CURVES_%s.jpg",fig_path, expt.name, o),
  #        plot.growth2,
  #        height = 140, width = two.c, units = 'mm',
  #        dpi = 600)
}

