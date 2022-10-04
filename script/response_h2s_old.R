
library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)

fig_path <- "~/R/Projects/methionine/paper/figures/"
`%notin%` <- Negate(`%in%`)

##### FIGURE SIZE
one.c <- 85 #single column
one.5c <- 114 #1.5 column
two.c <- 174 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 8

##### LOAD DATA
h2s.col <- read.csv(file = 'paper/data/bismuth_results2.csv')
h2s.col[is.na(h2s.col)] <- 0
h2s.col <- h2s.col %>%
  mutate(color_intensity = (255*2 - (red + green)))
h2s.col$color_intensity[h2s.col$blue == 0 | h2s.col$red == 0 | h2s.col$green == 0] <- 0
head(h2s.col)

h2s.col.sum <- h2s.col %>%
  group_by(condition, strain) %>%
  summarise(red = round(mean(red, na.rm = T)),
            blue = round(mean(blue, na.rm = T)),
            green = round(mean(green, na.rm = T)),
            color_intensity = round(mean(color_intensity, na.rm = T)),
            .groups = 'keep') %>%
  data.frame()

# h2s.col %>%
#   ggplot(aes(x = strain, y = color_intensity)) +
#   # geom_jitter(size = 0.5) +
#   geom_boxplot(aes(fill = as.factor(replicate))) +
#   facet_grid(.~condition) +
#   scale_fill_discrete(name = 'replicate')
# 
# 
# h2s.col %>%
#   ggplot(aes(x = condition, y = color_intensity)) +
#   # geom_jitter(aes(col = as.factor(replicate)), size = 0.5) +
#   geom_boxplot(aes(fill = as.factor(replicate))) +
#   stat_compare_means(method = 'anova') +
#   facet_wrap(.~strain) +
#   scale_fill_discrete(name = 'replicate') +
#   scale_color_discrete(name = 'replicate')


##### SIDE BY SIDE COLOR COMPARISON
h2s.col <- merge(h2s.col, h2s.col.sum %>% filter(strain == 'FY4-*met3Δ*'),
      by = 'condition', suffixes = c('','_control')) %>%
  # mutate(relative_color_intensity = color_intensity/color_intensity_control)
  mutate(relative_color_intensity = color_intensity)
h2s.col$strain <- factor(h2s.col$strain, levels = c('BY4742','BY4741','FY4','FY4-*met15Δ*',
                                                    'FY4-*met3Δ*','FY4-*met5Δ*',
                                                    'FY4-*met10Δ*','FY4-*met2Δ*','FY4-*met6Δ*','FY4-*met13Δ*',
                                                    'FY4-*cys4Δ*','FY4-*str3Δ*','FY4-*met12Δ*','FY4-*yll058wΔ*'))
h2s.col.sum <- merge(h2s.col.sum, h2s.col.sum %>% filter(strain == 'FY4-*met3Δ*'),
      by = 'condition', suffixes = c('','_control')) %>%
  # mutate(relative_color_intensity = color_intensity/color_intensity_control)
  mutate(relative_color_intensity = color_intensity)
h2s.col.sum$strain <- factor(h2s.col.sum$strain, levels = c('BY4742','BY4741','FY4','FY4-*met15Δ*',
                                                    'FY4-*met3Δ*','FY4-*met5Δ*',
                                                    'FY4-*met10Δ*','FY4-*met2Δ*','FY4-*met6Δ*','FY4-*met13Δ*',
                                                    'FY4-*cys4Δ*','FY4-*str3Δ*','FY4-*met12Δ*','FY4-*yll058wΔ*'))

 h2s.col.plot <- h2s.col %>%
  ggplot(aes(x = condition, y = relative_color_intensity)) +
  stat_summary(data = h2s.col.sum,
               aes(fill = relative_color_intensity), col = 'transparent', size = 1,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  # scale_y_continuous(breaks = seq(0,10,1)) +
  scale_fill_gradient2(high = "#5D4037", low = "white", guide = "none") +
  labs(y = 'Color Intensity (~H<sub>2</sub>S levels)',
       x = 'Condition') +
  stat_compare_means(size = 2, method = 'kruskal') +
  facet_wrap(.~strain) +
  theme_linedraw() +
  theme(plot.title = element_blank(), #element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 15, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(0.5,0,0.1,0, "mm"))) +
    coord_cartesian(ylim = c(0,500))
 ggsave(sprintf("%s/h2s_color_comparison.jpg",fig_path), h2s.col.plot,
        height = two.c, width = two.c, units = 'mm',
        bg = 'white',
        dpi = 300)


 h2s.col.plot.met15 <- h2s.col %>%
   filter(strain %in% c('FY4','FY4-*met15Δ*')) %>%
   ggplot(aes(x = condition, y = relative_color_intensity)) +
   stat_summary(data = h2s.col.sum %>%
                  filter(strain %in% c('FY4','FY4-*met15Δ*')),
                aes(fill = relative_color_intensity), col = 'transparent', size = 1,
                fun = mean, geom = "bar") +
   stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
   # scale_y_continuous(breaks = seq(0,10,1)) +
   scale_fill_gradient2(high = "#5D4037", low = "white", guide = "none") +
   labs(y = 'Color Intensity (~H<sub>2</sub>S levels)',
        x = 'Condition') +
   stat_compare_means(size = 2, method = 'kruskal') +
   facet_wrap(.~strain) +
   theme_linedraw() +
   theme(plot.title = element_blank(), #element_text(size = titles, face = 'bold', hjust = 0.5),
         axis.title.x = element_blank(),
         axis.title.y = ggtext::element_markdown(size = titles),
         axis.text = element_text(size = txt),
         axis.text.x = ggtext::element_markdown(size = txt, angle = 15, vjust = 1, hjust = 1),
         legend.title = element_text(size = titles),
         legend.text = element_text(size = txt),
         legend.position = 'bottom',
         legend.key.size = unit(3, "mm"),
         legend.box.spacing = unit(0.5,"mm"),
         strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                               margin = margin(0.5,0,0.1,0, "mm"))) +
   coord_cartesian(ylim = c(0,500))
 
 h2s.col.plot.met15 <- plot_grid(NULL, h2s.col.plot.met15,
                                 nrow = 1, 
                                 rel_widths = c(1,1),
                                 labels = c('A','B'),
                                 label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
 ggsave(sprintf("%s/h2s_color_comparison_met15.jpg",fig_path), h2s.col.plot.met15,
        height = one.c, width = two.c, units = 'mm',
        bg = 'white',
        dpi = 300)


##### FIGURE S5B
figS5b.2 <- h2s.col %>%
  filter(strain %notin% c('BY4742','BY4741','FY4-*yll058wΔ*')) %>%
  ggplot(aes(x = condition, y = relative_color_intensity)) +
  stat_summary(data = h2s.col.sum %>%
                 filter(strain %notin% c('BY4742','BY4741','FY4-*yll058wΔ*')),
               aes(fill = relative_color_intensity), col = 'transparent', alpha = 1, size = 1,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  stat_compare_means(method = 'kruskal', size = 1.8) +
  scale_fill_gradient2(high = "#5D4037", low = "white", guide = "none") +
  scale_x_discrete(labels = c('BiGGY' = 'BiGGY',
                              'SD-Met-Cys+Glu+Bi' = 'SD-Met\n-Cys+Glu\n+Bi')) +
  scale_y_continuous(breaks = seq(0,1000,150)) +
  labs(y = 'Color Intensity (~H<sub>2</sub>S levels)',
       x = 'Strains') +
  facet_wrap(.~strain, nrow = 2) +
  theme_linedraw() +
  theme(plot.title = element_blank(), #element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(1,0,0.1,0, "mm"))) #+
  # coord_flip(ylim = c(0,500))
  

figS5 <- plot_grid(annotate_figure(figS5a, top = text_grob("SD-Met-Cys+Glu+Bi Media", face = "bold", size = titles)),
                   figS5b.2,
                   ncol = 1, rel_heights = c(1,2),
                   labels = c('A','B'),
                   label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/FigureS5.jpg",fig_path), figS5,
       height = 200, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)

##### FIGURE S8B
figS8b <- h2s.col %>%
  filter(strain %in% c('FY4-*yll058wΔ*')) %>%
  ggplot(aes(x = condition, y = relative_color_intensity)) +
  stat_summary(data = h2s.col.sum %>%
                 filter(strain %in% c('FY4-*yll058wΔ*')),
               aes(fill = relative_color_intensity), col = 'transparent', alpha = 1, size = 1,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  stat_compare_means(method = 'kruskal', size = 2, angle = 270, hjust = 1) +
  scale_fill_gradient2(high = "#5D4037", low = "white", guide = "none") +
  labs(y = 'Color Intensity (~H<sub>2</sub>S levels)',
       x = '') +
  scale_x_discrete(limits = c('BiGGY','SD-Met-Cys+Glu+Bi'),
                   labels = c('BiGGY','SD-Met\n-Cys+Glu\n+Bi'),
                   position = 'top') +
  theme_linedraw() +
  theme(plot.title = element_blank(), #element_text(size = titles, face = 'bold', hjust = 0.5),
        plot.margin = margin(5,0,5,15),
        axis.title.y = element_blank(),
        axis.title.x = ggtext::element_markdown(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'black',
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_flip(ylim = c(0,500))



##### LOAD NEW DATA
h2s.col.new <- read.csv(file = 'paper/data/bismuth_results_new.csv')
h2s.col.new[is.na(h2s.col.new)] <- 0
h2s.col.new <- h2s.col.new %>%
  mutate(color_intensity = (255*2 - (red + green)))
h2s.col.new$color_intensity[h2s.col.new$blue == 0 | h2s.col.new$red == 0 | h2s.col.new$green == 0] <- 0
head(h2s.col.new)

h2s.col.new.sum <- h2s.col.new %>%
  group_by(condition, strain) %>%
  summarise(red = round(mean(red, na.rm = T)),
            blue = round(mean(blue, na.rm = T)),
            green = round(mean(green, na.rm = T)),
            color_intensity = round(mean(color_intensity, na.rm = T)),
            cs = mean(cs, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

##### SIDE BY SIDE COLOR COMPARISON
h2s.col.new <- merge(h2s.col.new, h2s.col.new.sum %>% filter(strain == 'FY4'),
                 by = 'condition', suffixes = c('','_control')) %>%
  # mutate(relative_color_intensity = color_intensity/color_intensity_control)
  mutate(relative_color_intensity = color_intensity)
h2s.col.new$strain <- factor(h2s.col.new$strain, levels = c('FY4',
                                                            'FY4-*met15Δ*',
                                                            'FY4-*met3Δ*'))
h2s.col.new$condition <- factor(h2s.col.new$condition, levels = c('SD-Met-Cys+Glu+Bi',
                                                                        'SD+Met(L)-Cys+Glu+Bi',
                                                                        'SD+Met(M)-Cys+Glu+Bi',
                                                                        'SD+Met(H)-Cys+Glu+Bi'))

h2s.col.new.sum <- merge(h2s.col.new.sum, h2s.col.new.sum %>% filter(strain == 'FY4'),
                     by = 'condition', suffixes = c('','_control')) %>%
  # mutate(relative_color_intensity = color_intensity/color_intensity_control)
  mutate(relative_color_intensity = color_intensity)
h2s.col.new.sum$strain <- factor(h2s.col.new.sum$strain, levels = c('FY4',
                                                                    'FY4-*met15Δ*',
                                                                    'FY4-*met3Δ*'))
h2s.col.new.sum$condition <- factor(h2s.col.new.sum$condition, levels = c('SD-Met-Cys+Glu+Bi',
                                                                                'SD+Met(L)-Cys+Glu+Bi',
                                                                                'SD+Met(M)-Cys+Glu+Bi',
                                                                                'SD+Met(H)-Cys+Glu+Bi'))


my_comparisons <- list(c('SD-Met-Cys+Glu+Bi','SD+Met(L)-Cys+Glu+Bi'),
                       c('SD+Met(L)-Cys+Glu+Bi','SD+Met(M)-Cys+Glu+Bi'),
                       c('SD+Met(M)-Cys+Glu+Bi','SD+Met(H)-Cys+Glu+Bi'),
                       c('SD-Met-Cys+Glu+Bi','SD+Met(M)-Cys+Glu+Bi'),
                       c('SD-Met-Cys+Glu+Bi','SD+Met(H)-Cys+Glu+Bi'))


h2s.col.new.plot <- h2s.col.new %>%
  ggplot(aes(x = condition, y = relative_color_intensity)) +
  stat_summary(data = h2s.col.new.sum,
               aes(fill = relative_color_intensity), col = 'transparent', size = 1,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_x_discrete(limits = c('SD-Met-Cys+Glu+Bi',
                              'SD+Met(L)-Cys+Glu+Bi',
                              'SD+Met(M)-Cys+Glu+Bi',
                              'SD+Met(H)-Cys+Glu+Bi')) +
  # scale_y_continuous(breaks = seq(0,1000,150)) +
  scale_fill_gradient2(high = "#5D4037", low = "white", guide = "none") +
  labs(y = 'Color Intensity (~H<sub>2</sub>S levels)',
       x = 'Condition') +
  stat_compare_means(size = 2, method = 'wilcox.test', comparisons = my_comparisons) +
  facet_wrap(.~strain) +
  theme_linedraw() +
  theme(plot.margin = margin(2,2,2,20),
        plot.title = element_blank(), #element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 15,  hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(0.5,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,650))
ggsave(sprintf("%s/h2s_color_new_strainwise.jpg",fig_path), h2s.col.new.plot,
       height = one.5c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


my_comparisons <- list(c('FY4','FY4-*met15Δ*'),
                       c('FY4','FY4-*met3Δ*'),
                       c('FY4-*met3Δ*','FY4-*met15Δ*'))

h2s.col.new.plot2 <- h2s.col.new %>%
  ggplot(aes(x = strain, y = relative_color_intensity)) +
  stat_summary(data = h2s.col.new.sum,
               aes(fill = relative_color_intensity), col = 'transparent', size = 1,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_x_discrete(limits = c('FY4',
                              'FY4-*met15Δ*',
                              'FY4-*met3Δ*')) +
  # scale_y_continuous(breaks = seq(0,1000,150)) +
  scale_fill_gradient2(high = "#5D4037", low = "white", guide = "none") +
  labs(y = 'Color Intensity (~H<sub>2</sub>S levels)',
       x = 'Strain') +
  stat_compare_means(size = 2, method = 'wilcox.test', comparisons = my_comparisons) +
  facet_wrap(.~condition, nrow = 1) +
  theme_linedraw() +
  theme(plot.margin = margin(2,2,2,20),
        plot.title = element_blank(), #element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 15,  hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(0.5,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,650))
ggsave(sprintf("%s/h2s_color_new_conditionwise.jpg",fig_path), h2s.col.new.plot2,
       height = one.5c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


h2s.col.new$condition <- factor(h2s.col.new$condition, levels = c('SD+Met(L)-Cys+Glu+Bi',
                                                                  'SD+Met(M)-Cys+Glu+Bi',
                                                                  'SD+Met(H)-Cys+Glu+Bi',
                                                                  'SD-Met-Cys+Glu+Bi'))
h2s.col.new.sum$condition <- factor(h2s.col.new.sum$condition, levels = c('SD+Met(L)-Cys+Glu+Bi',
                                                                          'SD+Met(M)-Cys+Glu+Bi',
                                                                          'SD+Met(H)-Cys+Glu+Bi',
                                                                          'SD-Met-Cys+Glu+Bi'))
h2s.col.new.plot.met15 <- h2s.col.new %>%
  filter(strain %in% c('FY4','FY4-*met15Δ*'),
         condition %in% c('SD+Met(H)-Cys+Glu+Bi','SD-Met-Cys+Glu+Bi')) %>%
  ggplot(aes(x = condition, y = relative_color_intensity)) +
  stat_summary(data = h2s.col.new.sum %>%
                 filter(strain %in% c('FY4','FY4-*met15Δ*'),
                        condition %in% c('SD+Met(H)-Cys+Glu+Bi','SD-Met-Cys+Glu+Bi')),
               aes(fill = relative_color_intensity), col = 'transparent', size = 1,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  stat_compare_means(size = 2, method = 'kruskal') +
  # scale_y_continuous(breaks = seq(0,10,1)) +
  scale_x_discrete(labels = c('SD+Met(H)-Cys+Glu+Bi' = 'SD+Met-Cys+Glu+Bi',
                              'SD-Met-Cys+Glu+Bi' = 'SD-Met-Cys+Glu+Bi')) +
  scale_fill_gradient2(high = "#5D4037", low = "white", guide = "none") +
  labs(y = 'Color Intensity (~H<sub>2</sub>S levels)',
       x = 'Condition') +
  
  facet_wrap(.~strain) +
  theme_linedraw() +
  theme(plot.title = element_blank(), #element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(size = txt, angle = 15, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(0.5,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,500))

# h2s.col.new.plot.met15 <- plot_grid(NULL, h2s.col.new.plot.met15,
#                                 nrow = 1, 
#                                 rel_widths = c(1,1),
#                                 labels = c('A','B'),
#                                 label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
# ggsave(sprintf("%s/h2s_color_new_comparison_met15.jpg",fig_path), h2s.col.new.plot.met15,
#        height = one.c, width = two.c, units = 'mm',
#        bg = 'white',
#        dpi = 300)


h2s.col.new.plot.met15 <- plot_grid(h2s.col.plot.met15, plot_grid(NULL, h2s.col.new.plot.met15,
                                                                      nrow = 1, 
                                                                      rel_widths = c(1,1),
                                                                      labels = c('C','D'),
                                                                      label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                                    nrow = 2)
ggsave(sprintf("%s/h2s_color_new_comparison_met15.jpg",fig_path), h2s.col.new.plot.met15,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


#### NEW AND OLD TOGETHER
colnames(h2s.col)
colnames(h2s.col.new)

h2s.col.both <- rbind(h2s.col %>%
                          filter(strain %in% c('FY4','FY4-*met15Δ*')),
                        h2s.col.new[,colnames(h2s.col)] %>%
                          filter(strain %in% c('FY4','FY4-*met15Δ*'),
                                 condition %in% c('SD+Met(H)-Cys+Glu+Bi')))
h2s.col.both$condition <- factor(h2s.col.both$condition, levels = c("SD-Met-Cys+Glu+Bi",
                                                                    "BiGGY",
                                                                    "SD+Met(H)-Cys+Glu+Bi"))

h2s.col.both.sum <- rbind(h2s.col.sum %>%
                            filter(strain %in% c('FY4','FY4-*met15Δ*')),
                          h2s.col.new.sum[,colnames(h2s.col.sum)] %>%
                            filter(strain %in% c('FY4','FY4-*met15Δ*'),
                                   condition %in% c('SD+Met(H)-Cys+Glu+Bi')))
h2s.col.both.sum$condition <- factor(h2s.col.both.sum$condition, levels = c("SD-Met-Cys+Glu+Bi",
                                                                    "BiGGY",
                                                                    "SD+Met(H)-Cys+Glu+Bi"))

my_comparisons <- list(c('SD+Met(H)-Cys+Glu+Bi','BiGGY'),
                       c('BiGGY','SD-Met-Cys+Glu+Bi'),
                       c('SD+Met(H)-Cys+Glu+Bi','SD-Met-Cys+Glu+Bi'))


h2s.col.both.plot.met15 <- h2s.col.both %>%
  ggplot(aes(x = condition, y = relative_color_intensity)) +
  stat_summary(data = h2s.col.both.sum,
               aes(fill = relative_color_intensity), col = 'transparent', size = 1,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  # stat_compare_means(size = 2, method = 'kruskal', comparisons = my_comparisons) +
  # scale_y_continuous(breaks = seq(0,10,1)) +
  scale_x_discrete(labels = c('SD+Met(H)-Cys+Glu+Bi' = 'SD+Met-Cys+Glu+Bi\n(++)',
                              'BiGGY' = 'BiGGY\n(+)',
                              'SD-Met-Cys+Glu+Bi' = 'SD-Met-Cys+Glu+Bi\n(-)')) +
  scale_fill_gradient2(high = "#5D4037", low = "white", guide = "none") +
  labs(y = 'Color Intensity (~H<sub>2</sub>S levels)',
       x = 'Condition') +
  facet_wrap(.~strain) +
  theme_linedraw() +
  theme(plot.title = element_blank(), #element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(0.5,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,500))

# h2s.col.both %>%
#   filter(strain == 'FY4-*met15Δ*', condition %in% c('SD+Met(H)-Cys+Glu+Bi','BiGGY')) %>%
#   kruskal_test(relative_color_intensity~condition)

r2r_fig1 <- plot_grid(NULL, h2s.col.both.plot.met15,
                                     nrow = 2, 
                                     rel_heights = c(1,2),
                                     labels = c('A','B'),
                                     label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/h2s_color_both_comparison_met15.jpg",fig_path), r2r_fig1,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)

#### CS vs Color
head(h2s.col.new.sum)

h2s.col.new.sum %>%
  ggplot(aes(x = cs, y = color_intensity)) +
  geom_point(aes(col = strain)) +
  facet_wrap(.~condition)


my_comparisons <- list(c('SD-Met-Cys+Glu+Bi','SD+Met(L)-Cys+Glu+Bi'),
                       c('SD+Met(L)-Cys+Glu+Bi','SD+Met(M)-Cys+Glu+Bi'),
                       c('SD+Met(M)-Cys+Glu+Bi','SD+Met(H)-Cys+Glu+Bi'),
                       c('SD+Met(L)-Cys+Glu+Bi','SD+Met(H)-Cys+Glu+Bi'),
                       c('SD-Met-Cys+Glu+Bi','SD+Met(M)-Cys+Glu+Bi'),
                       c('SD-Met-Cys+Glu+Bi','SD+Met(H)-Cys+Glu+Bi'))

h2s.col.new.cs.plot <- h2s.col.new %>%
  ggplot(aes(x = condition, y = cs)) +
  stat_summary(data = h2s.col.new.sum,
               aes(fill = relative_color_intensity), col = 'transparent', size = 1,
               fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_x_discrete(limits = c('SD-Met-Cys+Glu+Bi',
                              'SD+Met(L)-Cys+Glu+Bi',
                              'SD+Met(M)-Cys+Glu+Bi',
                              'SD+Met(H)-Cys+Glu+Bi')) +
  # scale_y_continuous(breaks = seq(0,1000,150)) +
  scale_fill_gradient2(high = "#5D4037", low = "white", guide = "none") +
  labs(y = 'Colony size (pixels)',
       x = 'Condition') +
  stat_compare_means(size = 2, method = 'wilcox.test', comparisons = my_comparisons) +
  facet_wrap(.~strain) +
  theme_linedraw() +
  theme(plot.margin = margin(2,2,2,20),
        plot.title = element_blank(), #element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = ggtext::element_markdown(size = txt),
        axis.text = element_text(size = txt),
        axis.text.x = element_text(angle = 15,  hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, colour = 'white',
                                              margin = margin(0.5,0,0.1,0, "mm")))
ggsave(sprintf("%s/h2s_color_new_cs_strainwise.jpg",fig_path), h2s.col.new.cs.plot,
       height = one.5c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)
