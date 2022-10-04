source('/home/sbp29/R/Projects/methionine/paper/scripts/initialize.R')
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")


fig1c_pm.dat <- dbGetQuery(conn, 'select b.orf_name, a.hours, avg(a.average) cs
                                  from Branden.carbon_rep1_PlM_Glu_384_CLEAN a, Branden.carbon_rep1_pos2orf_name b
                                  where a.pos = b.pos and b.orf_name != "BOR"
                                  group by a.hours, b.orf_name')


fig1c_mm.dat <- dbGetQuery(conn, 'select b.orf_name, a.hours, avg(a.average) cs
                                  from Branden.carbon_rep1_MiM_Glu_384_CLEAN a, Branden.carbon_rep1_pos2orf_name b
                                  where a.pos = b.pos and b.orf_name != "BOR"
                                  group by a.hours, b.orf_name')

fig1c_mm.dat %>%
  filter(orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = hours, y = cs)) +
  geom_line(aes(col = orf_name)) +
  labs(x = 'time (hrs)', y = 'colony size (pixels)') +
  scale_color_discrete(name = '') +
  theme_linedraw() +
  theme(legend.position = 'bottom')

merge(fig1c_mm.dat,
      fig1c_mm.dat[fig1c_mm.dat$orf_name == 'FY4',], by = c('hours')) %>%
  filter(orf_name.x != 'FY4-met3del') %>%
  mutate(rel_fit = cs.x/cs.y) %>%
  ggplot(aes(x = hours, y = rel_fit)) +
  geom_line(aes(col = orf_name.x)) +
  geom_vline(xintercept = 48, linetype = 'dashed', size = 0.5) +
  labs(x = 'time (hrs)', y = 'relative colony size') +
  scale_color_discrete(name = '') +
  theme_linedraw() +
  theme(legend.position = 'bottom')



fig1c_mu.dat <- dbGetQuery(conn, 'select b.orf_name, a.hours, avg(a.average) cs
                                  from Branden.carbon_rep2_MiU_Glu_384_CLEAN a, Branden.carbon_rep1_pos2orf_name b
                                  where a.pos = b.pos and b.orf_name != "BOR"
                                  group by a.hours, b.orf_name')

fig1c_mu.dat %>%
  filter(orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = hours, y = cs)) +
  geom_line(aes(col = orf_name)) +
  labs(x = 'time (hrs)', y = 'colony size (pixels)') +
  scale_color_discrete(name = '') +
  theme_linedraw() +
  theme(legend.position = 'bottom')

merge(fig1c_mu.dat,
      fig1c_mu.dat[fig1c_mu.dat$orf_name == 'FY4',], by = c('hours')) %>%
  filter(orf_name.x != 'FY4-met3del') %>%
  mutate(rel_fit = cs.x/cs.y) %>%
  ggplot(aes(x = hours, y = rel_fit)) +
  geom_line(aes(col = orf_name.x)) +
  geom_vline(xintercept = 48, linetype = 'dashed', size = 0.5) +
  labs(x = 'time (hrs)', y = 'relative colony size') +
  scale_color_discrete(name = '') +
  theme_linedraw() +
  theme(legend.position = 'bottom')


gc_cs <- rbind(fig1c_mm.dat %>%
        mutate(condition = 'SD-Met-Cys+Glu'),
      fig1c_mu.dat %>%
        mutate(condition = 'SD-Ura+Glu')) %>%
  filter(orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = hours, y = cs)) +
  geom_line(aes(col = orf_name), size = 1) +
  geom_vline(xintercept = 48, linetype = 'dashed', size = 0.5) +
  annotate(geom = 'text', x = 90, y = 11500, label = '2 days', size = 3) +
  labs(x = 'time (hrs)', y = 'colony size (pixels)') +
  scale_color_discrete(name = '') +
  facet_wrap(.~condition) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt, margin = margin(1,0,1,0, "mm")))

gc_rcs <- rbind(merge(fig1c_mm.dat,
            fig1c_mm.dat[fig1c_mm.dat$orf_name == 'FY4',], by = c('hours')) %>%
        filter(orf_name.x != 'FY4-met3del') %>%
        mutate(condition = 'SD-Met-Cys+Glu'),
      
      merge(fig1c_mu.dat,
            fig1c_mu.dat[fig1c_mu.dat$orf_name == 'FY4',], by = c('hours')) %>%
        filter(orf_name.x != 'FY4-met3del') %>%
        mutate(condition = 'SD-Ura+Glu')) %>%
  mutate(rel_fit = cs.x/cs.y) %>%
  ggplot(aes(x = hours, y = rel_fit)) +
  geom_line(aes(col = orf_name.x), size = 1) +
  geom_vline(xintercept = 48, linetype = 'dashed', size = 0.5) +
  annotate(geom = 'text', x = 90, y = 1.43, label = '2 days', size = 3) +
  scale_y_continuous(breaks = seq(-2,2,0.2)) +
  labs(x = 'time (hrs)', y = 'relative colony size') +
  scale_color_discrete(name = '') +
  facet_wrap(.~condition) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt, margin = margin(1,0,1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))


gc_cs_rcs <- ggpubr::ggarrange(gc_cs, gc_rcs,
                  ncol = 1, align = 'hv',
          legend = 'bottom', common.legend = T)

ggsave(sprintf("%s/minus_met_ura_growth_curves.jpg",fig_path), gc_cs_rcs,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


##### COLONY SIZE
fig1c_pm.dat$background[fig1c_pm.dat$orf_name %in% c('BY4741','FY4-met15del')] <- 'met15D'
fig1c_pm.dat$background[is.na(fig1c_pm.dat$background)] <- 'met15'
fig1c_pm.dat$aux[fig1c_pm.dat$orf_name == 'FY4'] <- 'FY4 (Met15 | Ura3)'
fig1c_pm.dat$aux[fig1c_pm.dat$orf_name == 'FY4-met15del'] <- 'FY4-*met15Δ* (*met15Δ* | Ura3)'
fig1c_pm.dat$aux[fig1c_pm.dat$orf_name == 'BY4741'] <- 'BY4741 (*met15Δ* | *ura3Δ*)'
fig1c_pm.dat$aux[fig1c_pm.dat$orf_name == 'BY4742'] <- 'BY4742 (Met15 | *ura3Δ*)'

fig1c_mm.dat$background[fig1c_mm.dat$orf_name %in% c('BY4741','FY4-met15del')] <- 'met15D'
fig1c_mm.dat$background[is.na(fig1c_mm.dat$background)] <- 'met15'
fig1c_mm.dat$aux[fig1c_mm.dat$orf_name == 'FY4'] <- 'FY4 (Met15 | Ura3)'
fig1c_mm.dat$aux[fig1c_mm.dat$orf_name == 'FY4-met15del'] <- 'FY4-*met15Δ* (*met15Δ* | Ura3)'
fig1c_mm.dat$aux[fig1c_mm.dat$orf_name == 'BY4741'] <- 'BY4741 (*met15Δ* | *ura3Δ*)'
fig1c_mm.dat$aux[fig1c_mm.dat$orf_name == 'BY4742'] <- 'BY4742 (Met15 | *ura3Δ*)'

fig1c_mu.dat$background[fig1c_mu.dat$orf_name %in% c('BY4741','BY4742')] <- 'ura3D'
fig1c_mu.dat$background[is.na(fig1c_mu.dat$background)] <- 'ura3'
fig1c_mu.dat$aux[fig1c_mu.dat$orf_name == 'FY4'] <- 'FY4 (Met15 | Ura3)'
fig1c_mu.dat$aux[fig1c_mu.dat$orf_name == 'FY4-met15del'] <- 'FY4-*met15Δ* (*met15Δ* | Ura3)'
fig1c_mu.dat$aux[fig1c_mu.dat$orf_name == 'BY4741'] <- 'BY4741 (*met15Δ* | *ura3Δ*)'
fig1c_mu.dat$aux[fig1c_mu.dat$orf_name == 'BY4742'] <- 'BY4742 (Met15 | *ura3Δ*)'

mmmmu_gc <- rbind(fig1c_mm.dat %>%
        mutate(condition = 'SD-Met-Cys+Glu', saturation = 366),
      fig1c_mu.dat %>%
        mutate(condition = 'SD-Ura+Glu', saturation = 358)) %>%
  filter(orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = hours, y = cs)) +
  geom_line(aes(col = aux), size = 2) +
  geom_vline(xintercept = 48, linetype = 'dashed', size = 0.5) +
  annotate(geom = 'text', x = 90, y = 11500, label = '2 days', size = 3) +
  geom_vline(xintercept = 366, linetype = 'dashed', size = 0.5) +
  annotate(geom = 'text', x = 330, y = 11500, label = '15 days', size = 3) +
  labs(x = 'Time (hrs)', y = 'Colony size (pixels)') +
  scale_color_manual(name = '',
                     values = c('FY4 (Met15 | Ura3)' = '#FFC107',
                                'FY4-*met15Δ* (*met15Δ* | Ura3)' = '#536DFE',
                                'BY4741 (*met15Δ* | *ura3Δ*)' = '#990000',
                                'BY4742 (Met15 | *ura3Δ*)' = '#FF0000')) +
  facet_wrap(condition~background,
             labeller = labeller(background = c('met15' = 'Methionine Prototroph',
                                         'met15D' = 'Presumed Methionine Auxotroph',
                                         'ura3' = 'Uracil Prototroph',
                                         'ura3D' = 'Presumed Uracil Auxotroph'))) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt, margin = margin(1,0,1,0, "mm"))) +
  guides(color = guide_legend(nrow=2, byrow=TRUE, order = 1,
                              override.aes = list(size = 3)))

ggsave(sprintf("%s/minus_met_ura_growth_curves.jpg",fig_path), mmmmu_gc,
       height = two.c, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)


##### NEW FIGURE S3
figS3b

figS3c.right <- rbind(fig1c_mm.dat %>%
        mutate(condition = 'SD-Met-Cys+Glu', saturation = 366),
      fig1c_mu.dat %>%
        mutate(condition = 'SD-Ura+Glu', saturation = 358)) %>%
  filter(orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = hours, y = cs)) +
  geom_line(aes(col = aux), size = 2) +
  geom_vline(xintercept = 48, linetype = 'dashed', size = 0.5) +
  annotate(geom = 'text', x = 90, y = 11500, label = '2 days', size = 3) +
  geom_vline(xintercept = 366, linetype = 'dashed', size = 0.5) +
  annotate(geom = 'text', x = 330, y = 11500, label = '15 days', size = 3) +
  labs(x = 'Time (hours)', y = 'Colony Size (pixels)') +
  scale_color_manual(name = '',
                     values = c('FY4 (Met15 | Ura3)' = '#FFC107',
                                'FY4-*met15Δ* (*met15Δ* | Ura3)' = '#536DFE',
                                'BY4741 (*met15Δ* | *ura3Δ*)' = '#990000',
                                'BY4742 (Met15 | *ura3Δ*)' = '#FF0000')) +
  facet_wrap(condition~background,
             labeller = labeller(background = c('met15' = 'Methionine\nPrototroph',
                                                'met15D' = 'Presumed Methionine\nAuxotroph',
                                                'ura3' = 'Uracil\nPrototroph',
                                                'ura3D' = 'Presumed Uracil\nAuxotroph')),
             nrow = 2) +
  theme_linedraw() +
  theme(panel.background = element_rect(fill = pinningpanelbg, color = 'transparent'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = titles),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt, margin = margin(-0.1,0,0,0,"mm"))) +
  guides(color = guide_legend(nrow=1, byrow=TRUE, order = 1,
                              override.aes = list(size = 3)))

figS3c.left <- fig1c_pm.dat %>%
                        mutate(condition = 'SD+Met-Cys+Glu', saturation = 366) %>%
  filter(orf_name != 'FY4-met3del') %>%
  ggplot(aes(x = hours, y = cs)) +
  geom_line(aes(col = aux), size = 2) +
  geom_vline(xintercept = 48, linetype = 'dashed', size = 0.5) +
  annotate(geom = 'text', x = 90, y = 11500, label = '2 days', size = 3) +
  geom_vline(xintercept = 366, linetype = 'dashed', size = 0.5) +
  annotate(geom = 'text', x = 330, y = 11500, label = '15 days', size = 3) +
  labs(x = 'Time (hours)', y = 'Colony size (pixels)') +
  scale_color_manual(name = '',
                     values = c('FY4 (Met15 | Ura3)' = '#FFC107',
                                'FY4-*met15Δ* (*met15Δ* | Ura3)' = '#536DFE',
                                'BY4741 (*met15Δ* | *ura3Δ*)' = '#990000',
                                'BY4742 (Met15 | *ura3Δ*)' = '#FF0000')) +
  facet_wrap(.~condition,
             labeller = labeller(background = c('met15' = 'Methionine\nPrototroph',
                                                'met15D' = 'Presumed Methionine\nAuxotroph',
                                                'ura3' = 'Uracil\nPrototroph',
                                                'ura3D' = 'Presumed Uracil\nAuxotroph')),
             nrow = 2) +
  theme_linedraw() +
  theme(panel.background = element_rect(fill = pinningpanelbg, color = 'transparent'),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.margin = margin(0,0,0,0),
        legend.title = element_text(size = titles),
        legend.text = ggtext::element_markdown(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt, margin = margin(-0.1,0,0,0,"mm"))) +
  guides(color = guide_legend(nrow=1, byrow=TRUE, order = 1,
                              override.aes = list(size = 3)))


figS3c <- ggpubr::ggarrange(figS3c.left, figS3c.right,
                            nrow = 1, widths = c(1,2),
                            common.legend = T, legend = 'bottom',
                            labels = c('B','C'),
                            font.label = list(size = lbls, face = 'bold', family = 'sans'))

figS3 <- plot_grid(NULL,
                   figS3b,
                   NULL,
                   # plot_grid(NULL,NULL,
                   #           labels = c('C','D'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
                   #           rel_widths = c(1,2), nrow = 1),
                   labels = c('A','B',''),
                   label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
                   ncol = 1, rel_heights = c(1.05,1.25,1.3))
ggsave(sprintf("%s/FigureS3.jpg",fig_path), figS3,
       height = 230, width = two.c, units = 'mm',
       bg = 'white',
       dpi = 300)
