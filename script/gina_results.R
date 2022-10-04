source('/home/sbp29/R/Projects/methionine/paper/scripts/initialize.R')

pzfx_tables("/home/sbp29/R/Projects/methionine/data/gina/072321_gina_data.pzfx")
data.gina <- read_pzfx("/home/sbp29/R/Projects/methionine/data/gina/072321_gina_data.pzfx", table = "Data 2")
head(data.gina)

#####
data.gina <- melt(data.gina, id.vars = 'Time (min)', variable.name = 'ID', value.name = 'uM')
data.gina <- cbind(data.gina, str_split(data.gina$ID, '_', simplify = T))
colnames(data.gina) <- c(colnames(data.gina)[1:3], 'Sample', 'Replicate')
data.gina$Sample <- factor(data.gina$Sample, levels = c('Met15', 'Yll058w', 'None'))

data.gina$outlier <- NULL
for (t in unique(data.gina$`Time (min)`)) {
  for (s in unique(data.gina$Sample[data.gina$`Time (min)` == t])) {
    data.gina$outlier[data.gina$`Time (min)` == t & data.gina$Sample == s] <-
      isoutlier(data.gina$uM[data.gina$`Time (min)` == t & data.gina$Sample == s],3)
  }
}

#####
ttest.res <- NULL
for (s1 in unique(data.gina$Sample)) {
  for (s2 in unique(data.gina$Sample[data.gina$Sample != s1])) {
    temp.t <- t.test(uM ~ Sample, data = data.gina %>%
                       filter(Sample %in% c(s1, s2)) %>%
                       group_by(`Time (min)`, Sample) %>%
                       summarise(.groups = 'keep', uM = mean(uM, na.rm = T)) %>%
                       data.frame(), paired = TRUE, alternative = 'greater')
    ttest.res <- rbind(ttest.res, data.frame(reference = s1, query = s2, p = temp.t$p.value, stat = temp.t$statistic[[1]]))
  }
}
head(ttest.res)

ttest.res$label[ttest.res$p > 0.05] <- 'ns'
ttest.res$label[ttest.res$p <= 0.05] <- '*'
ttest.res$label[ttest.res$p <= 0.01] <- '**'
ttest.res$label[ttest.res$p <= 0.001] <- '***'
ttest.res$label[ttest.res$p <= 0.0001] <- '****'

ttest.res <- merge(ttest.res, data.gina %>%
                     filter(`Time (min)` == 15) %>%
                     group_by(Sample) %>%
                     summarise(.groups = 'keep', uM = mean(uM, na.rm = T)) %>%
                     data.frame(), by.x = 'query', by.y = 'Sample')

# #####
# for (t in unique(data.gina$`Time (min)`)) {
#   for (s in unique(data.gina$Sample[data.gina$`Time (min)` == t & data.gina$Sample != 'None'])) {
#     temp.kw <- compare_means(data = data.gina[data.gina$`Time (min)` == t & data.gina$Sample %in% c(s,'None') &
#                                                 data.gina$outlier == FALSE,],
#                                 method = 't.test', formula = uM ~ Sample) %>% data.frame()
#     
#   }
# }

#####
plot.gina.res <- data.gina %>%
  # filter(outlier == FALSE) %>%
  ggplot(aes(x = `Time (min)`, y = uM)) +
  stat_summary(aes(group = Sample, col = Sample), fun=mean, geom="line", lwd = 0.7) +
  stat_summary(aes(group = Sample, col = Sample), fun.data = mean_se, geom = "errorbar", lwd = 0.7) +
  stat_summary(aes(group = Sample), fun=mean, geom="point", size =2) +
  stat_summary(aes(group = Sample, col = Sample), fun=mean, geom="point", size = 0.5) +
  labs(y = 'Homocysteine (log2(μM))') +
  scale_color_manual(name = 'Sample',
                     values = c('Met15' = "#D32F2F",
                                'Yll058w' = "#4CAF50",
                                'None' = "#1976D2")) +
  scale_y_continuous(trans = 'log2') +
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
  coord_cartesian(ylim = c(5,125))
fig6C <- plot.gina.res
save(fig6C, file = 'figures/final/fig6C.RData')
# ggsave(sprintf("%s/%s/GINA_RESULTS.jpg",fig_path, expt.name), plot.gina.res,
#        height = one.c, width = one.c, units = 'mm',
#        dpi = 600)



##### ENZYME KINETICS
oah_data <- read.csv('/home/sbp29/R/Projects/methionine/data/gina/OAH_data.csv')
kin_res <- read.csv('/home/sbp29/R/Projects/methionine/data/gina/KmVmax.csv')
plp_res <- read.csv('/home/sbp29/R/Projects/methionine/data/gina/PLP_dependence.csv')


##### FIGURE 3F
fig3f <- oah_data %>%
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


#### PLP DEPENDENCE
head(plp_res)
plp_res$Sample2 <- paste0(plp_res$Group, plp_res$Sample)
plp_res$Sample2 <- factor(plp_res$Sample2, levels = c("wPLPMet15", "woPLPMet15", "wPLPYll058w", "woPLPYll058w", "wPLPNone"))
unique(plp_res$Sample2)

plp_res$Group <- factor(plp_res$Group, levels = c('wPLP','woPLP'))
plp_res$Sample <- factor(plp_res$Sample, levels = c('Met15','Yll058w','None'))


## STATS
plp_res.stats <- NULL
i <- 1
for (s1 in unique(plp_res$Sample2)) {
  for (s2 in unique(plp_res$Sample2[plp_res$Sample2 != s1])) {
    t1 <- plp_res$Rate[plp_res$Sample2 == s1]
    t2 <- plp_res$Rate[plp_res$Sample2 == s2]
    es <- abs(median(t2) - median(t1))/median(t1)
    test <- t.test(t1, t2)
    plp_res.stats$sample[i] <- s2
    plp_res.stats$ref[i] <- s1
    plp_res.stats$es[i] <- es
    plp_res.stats$statistic[i] <- test$statistic
    plp_res.stats$pvalue[i] <- test$p.value
    i <- i + 1
  }
}
plp_res.stats <- data.frame(plp_res.stats)


plpnoplp <- rbind(cbind(Sample = 'Yll058w',
                        Ratio = rbind(transform(plp_res$Rate[plp_res$Sample2 == 'wPLPYll058w']/plp_res$Rate[plp_res$Sample2 == 'woPLPYll058w'][1]),
                                      transform(plp_res$Rate[plp_res$Sample2 == 'wPLPYll058w']/plp_res$Rate[plp_res$Sample2 == 'woPLPYll058w'][2]),
                                      transform(plp_res$Rate[plp_res$Sample2 == 'wPLPYll058w']/plp_res$Rate[plp_res$Sample2 == 'woPLPYll058w'][3]))),
                  cbind(Sample = 'Met15',
                        Ratio = rbind(transform(plp_res$Rate[plp_res$Sample2 == 'wPLPMet15']/plp_res$Rate[plp_res$Sample2 == 'woPLPMet15'][1]),
                                      transform(plp_res$Rate[plp_res$Sample2 == 'wPLPMet15']/plp_res$Rate[plp_res$Sample2 == 'woPLPMet15'][2]),
                                      transform(plp_res$Rate[plp_res$Sample2 == 'wPLPMet15']/plp_res$Rate[plp_res$Sample2 == 'woPLPMet15'][3]))))
colnames(plpnoplp) <- c('Sample','Ratio')


plpnoplp %>%
  mutate(lfc = log2(Ratio)) %>%
  ggplot(aes(x = Sample, y = lfc, fill = Sample)) +
  stat_summary(aes(fill = Sample), size = 1, fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  # scale_y_log10() +
  scale_y_continuous(breaks = seq(0,10,1)) +
  labs(y = 'Log<sub>2</sub>FoldChange(PLP/no PLP)') +
  scale_fill_manual(name = 'Sample',
                    values = c("Met15" = "#330099",
                               "Yll058w" = "#CC33FF",
                               "None" = "#CC99CC"),
                    guide = 'none') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        # axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = titles),
        axis.text = element_text(size = txt),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,6))


## FIGURE
figSPLP.top <- plp_res %>%
  ggplot(aes(x = Group, y = Rate)) +
  stat_summary(aes(fill = Sample), size = 1, fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_y_continuous(breaks = seq(0,100,1)) +
  scale_fill_manual(name = 'Sample',
                    values = c("Met15" = "#330099",
                               "Yll058w" = "#CC33FF",
                               "None" = "#CC99CC"),
                    guide = 'none') +
  # scale_alpha_manual(name = 'Group',
  #                   values = c("wPLP" = 1,
  #                              "woPLP" = 0.8)) +
  scale_x_discrete(labels = c('wPLP' = 'With PLP',
                              'woPLP' = 'Without PLP')) +
  stat_compare_means(method = 't.test', size = 2) +
  labs(x = 'PLP', y = 'Rate (μM/min)') +
  facet_grid(.~Sample, scales = 'free_x', space = 'free') +
  theme_minimal() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title.y = element_text(size = titles, hjust = 0),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text = element_text(size = txt),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(55,61))

figSPLP.bottom <- plp_res %>%
  ggplot(aes(x = Group, y = Rate)) +
  stat_summary(aes(fill = Sample), size = 1, fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_y_continuous(breaks = seq(0,100,1)) +
  scale_fill_manual(name = 'Sample',
                    values = c("Met15" = "#330099",
                               "Yll058w" = "#CC33FF",
                               "None" = "#CC99CC")) +
  # scale_alpha_manual(name = 'Group',
  #                   values = c("wPLP" = 1,
  #                              "woPLP" = 0.8)) +
  scale_x_discrete(labels = c('wPLP' = 'With PLP',
                              'woPLP' = 'Without PLP')) +
  # stat_compare_means(method = 't.test', size = 2) +
  labs(x = 'PLP', y = 'Rate (μM/min)') +
  facet_grid(.~Sample, scales = 'free_x', space = 'free') +
  theme_minimal() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        # axis.title = element_text(size = titles),
        axis.title = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        legend.margin = margin(0,0,-5,0),
        strip.text = element_text(size = txt, color = 'white',
                                  face = 'bold',
                                  margin = margin(-3,0,0,0, "mm"))) +
  guides(fill = guide_legend(override.aes=list(shape = 15, alpha = 1))) +
  coord_cartesian(ylim = c(0,3))

figSPLP <- plot_grid(figSPLP.top, figSPLP.bottom,
          ncol = 1, 
          rel_heights = c(2,1.1),
          align = 'v')

ggsave(sprintf("%s/FigureSPLP_cropped.jpg",fig_path), figSPLP,
       height = one.c, width = one.5c, units = 'mm',
       bg = 'white',
       dpi = 300)  


# figSPLP <- 
plp_res %>%
  ggplot(aes(x = Group, y = Rate)) +
  stat_summary(aes(fill = Sample), size = 1, fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_y_continuous(breaks = seq(0,100,1)) +
  scale_fill_manual(name = 'Sample',
                    values = c("Met15" = "#330099",
                               "Yll058w" = "#CC33FF",
                               "None" = "#CC99CC")) +
  # scale_alpha_manual(name = 'Group',
  #                   values = c("wPLP" = 1,
  #                              "woPLP" = 0.8)) +
  scale_x_discrete(labels = c('wPLP' = 'With PLP',
                              'woPLP' = 'Without PLP')) +
  # stat_compare_means(method = 't.test', size = 2) +
  labs(x = 'PLP', y = 'Rate (μM/min)') +
  facet_grid(.~Sample, scales = 'free_x', space = 'free') +
  theme_minimal() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        # axis.title = element_text(size = titles),
        axis.title = element_blank(),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        legend.margin = margin(0,0,-5,0),
        strip.text = element_text(size = txt, color = 'white',
                                  face = 'bold',
                                  margin = margin(-3,0,0,0, "mm"))) +
  guides(fill = guide_legend(override.aes=list(shape = 15, alpha = 1)))

figSPLP <- plot_grid(figSPLP.top, figSPLP.bottom,
                     ncol = 1, 
                     rel_heights = c(2,1.1),
                     align = 'v')

ggsave(sprintf("%s/FigureSPLP_cropped.jpg",fig_path), figSPLP,
       height = one.c, width = one.5c, units = 'mm',
       bg = 'white',
       dpi = 300)  
  
  
  