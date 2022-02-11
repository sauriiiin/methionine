
# data.ns2 <- merge(data.ns[!(data.ns$id == 'Agar_Difco_+sulfate_+Met' & data.ns$stage == 'S1'),] %>%
#                     filter(id %in% c('Agar_Difco_+sulfate_+Met',
#                                      'Agarose_Home_+sulfate_-Met',
#                                      'Agarose_Home_-sulfate_-Met'),
#                            orf_name %in% c('FY4', 'FY4_met15del','FY4_met3del')) %>%
#                     # group_by(stage, base, ynb_type, sulfate, methionine, orf_name, hours, cum_hrs, id) %>%
#                     # summarize(average = median(average, na.rm = T), .groups = 'keep') %>%
#                     data.frame(),
#                   data.ns[!(data.ns$id == 'Agar_Difco_+sulfate_+Met' & data.ns$stage == 'S1'),] %>%
#                     filter(id %in% c('Agar_Difco_+sulfate_+Met',
#                                      'Agarose_Home_+sulfate_-Met',
#                                      'Agarose_Home_-sulfate_-Met'),
#                            orf_name %in% c('FY4', 'FY4_met15del','FY4_met3del')) %>%
#                     group_by(stage, id) %>%
#                     summarize(saturation = max(hours, na.rm = T), .groups = 'keep') %>%
#                     data.frame(), by = c('stage','id')) %>%
#   filter(hours == saturation, average != 0)

data.ns2 <- data.ns[!(data.ns$id == 'Agar_Difco_+sulfate_+Met' & data.ns$stage == 'S1'),] %>%
  filter(id %in% c('Agar_Difco_+sulfate_+Met',
                   'Agarose_Home_+sulfate_-Met',
                   'Agarose_Home_-sulfate_-Met'),
         orf_name %in% c('FY4', 'FY4_met15del','FY4_met3del')) %>%
  data.frame()
data.ns2$id <- factor(data.ns2$id, levels = c('Agarose_Home_+sulfate_-Met',
                                              'Agarose_Home_-sulfate_-Met',
                                              'Agar_Difco_+sulfate_+Met'))
# data.ns2 <- data.ns2 %>% filter(average != 0)

fig3d <- data.ns2 %>%
  ggplot(aes(x = cum_hrs, y = average)) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S1'),
               aes(fill = id),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S1'),
               aes(col = id),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S2'),
               aes(fill = id),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S2'),
               aes(col = id),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S3'),
               aes(fill = id),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S3'),
               aes(col = id),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S4'),
               aes(fill = id),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S4'),
               aes(col = id),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S5'),
               aes(fill = id),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S5'),
               aes(col = id),
               fun=mean, geom="line", lwd =1) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S6'),
               aes(fill = id),
               fun.data=mean_sdl, fun.args = list(mult=1), geom="ribbon", alpha = 0.3) +
  stat_summary(data = data.ns2 %>% filter(stage == 'S6'),
               aes(col = id),
               fun=mean, geom="line", lwd =1) +
  scale_fill_manual(name = 'Condition',
                    breaks = c('Agarose_Home_-sulfate_-Met',
                               'Agarose_Home_+sulfate_-Met',
                               'Agar_Difco_+sulfate_+Met'),
                    values = c('Agar_Difco_+sulfate_+Met' = '#9E9E9E',
                               'Agarose_Home_+sulfate_-Met' = '#212121',
                               'Agarose_Home_-sulfate_-Met' = '#FF5722'),
                    labels = c('Agar_Difco_+sulfate_+Met' = 'SD+Met+Glu w/\nInorganic Sulfates',
                               'Agarose_Home_+sulfate_-Met' = 'SD-Met+Glu w/\nInorganic Sulfates',
                               'Agarose_Home_-sulfate_-Met' = 'SD-Met+Glu w/o\nInorganic Sulfates'),
                    guide = F) +
  scale_color_manual(name = 'Condition',
                    breaks = c('Agarose_Home_-sulfate_-Met',
                               'Agarose_Home_+sulfate_-Met',
                               'Agar_Difco_+sulfate_+Met'),
                    values = c('Agar_Difco_+sulfate_+Met' = '#9E9E9E',
                               'Agarose_Home_+sulfate_-Met' = '#212121',
                               'Agarose_Home_-sulfate_-Met' = '#FF5722'),
                    labels = c('Agar_Difco_+sulfate_+Met' = 'SD+Met+Glu w/ Inorganic Sulfates',
                               'Agarose_Home_+sulfate_-Met' = 'SD-Met+Glu w/ Inorganic Sulfates',
                               'Agarose_Home_-sulfate_-Met' = 'SD-Met+Glu w/o Inorganic Sulfates')) +
  # scale_x_discrete(breaks = c('FY4','FY4_met15del','FY4_met3del'),
  #                  labels = c('FY4' = 'FY4',
  #                             'FY4_met15del' = 'FY4-*met15Δ*',
  #                             'FY4_met3del' = 'FY4-*met3Δ*')) +
  scale_x_continuous(breaks = seq(-100,1000,100)) +
  scale_y_continuous(breaks = seq(-4500,13500,4500)) +
  facet_wrap(.~orf_name, ncol = 1,
             labeller = labeller(stage = c('S1' = 'Pin #1',
                                           'S2' = 'Pin #2',
                                           'S3' = 'Pin #3',
                                           'S4' = 'Pin #4',
                                           'S5' = 'Pin #5',
                                           'S6' = 'Pin #6'),
                                 orf_name = c('FY4' = 'FY4',
                                              'FY4_met15del' = 'FY4-*met15Δ*',
                                              'FY4_met3del' = 'FY4-*met3Δ*'))) +
  labs(x = 'Time (hours)',
       y = 'Colony Size (pixels)') +
  theme_minimal() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = ggtext::element_markdown(size = txt, angle = 30, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles, hjust = 0.5),
        legend.text = element_text(size = txt),
        # legend.position = c(0.74, 0.3),
        legend.position = 'bottom',
        legend.box = 'vertical',
        # legend.background = element_rect(fill="white",
        #                                  size=0.4, linetype="solid",
        #                                  colour ="black"),
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        legend.spacing.y = unit(1, "mm"),
        panel.spacing.x = unit(0, "lines"),
        strip.text = ggtext::element_markdown(size = txt, 
                                  margin = margin(0,0,0,0, "mm")),
        strip.text.y = ggtext::element_markdown(size = txt,
                                  margin = margin(0,1,0,0, "mm"))) +
  guides(color = guide_legend(nrow=1, byrow=FALSE, order = 1, 
                             override.aes = list(size = 3))) +
  coord_cartesian(xlim = c(0,910))


data.ns2 %>%
  group_by(stage, condition) %>%
  summarize(max_hrs = max(hours, na.rm = T), .groups = 'keep') %>%
  data.frame()


fig3d.expt <- data.ns2 %>%
  filter(id %in% c('Agarose_Home_-sulfate_-Met',
                           'Agarose_Home_+sulfate_-Met',
                           'Agar_Difco_+sulfate_+Met')) %>%
  # group_by(stage, id, cum_hrs) %>%
  # count() %>%
  group_by(stage, id) %>%
  summarize(med_hrs = median(cum_hrs, na.rm = T), .groups = 'keep') %>%
  data.frame() %>%
  ggplot(aes(x = med_hrs, y = id)) +
  geom_point(aes(col = id), 
             shape = 15, size = 5) +
  geom_point(aes(x = med_hrs - 8, col = id), 
             shape = 15, size = 5) +
  geom_point(aes(x = med_hrs + 8, col = id), 
             shape = 15, size = 5) +
  geom_point(size = .5) +
  scale_y_discrete(breaks = c('Agarose_Home_-sulfate_-Met',
                               'Agarose_Home_+sulfate_-Met',
                               'Agar_Difco_+sulfate_+Met'),
                   limits = c('Agar_Difco_+sulfate_+Met',
                              'Agarose_Home_-sulfate_-Met',
                              'Agarose_Home_+sulfate_-Met')) +
  scale_color_manual(name = 'Condition',
                    breaks = c('Agarose_Home_-sulfate_-Met',
                               'Agarose_Home_+sulfate_-Met',
                               'Agar_Difco_+sulfate_+Met'),
                    values = c('Agar_Difco_+sulfate_+Met' = '#9E9E9E',
                               'Agarose_Home_+sulfate_-Met' = '#212121',
                               'Agarose_Home_-sulfate_-Met' = '#FF5722'),
                    labels = c('Agar_Difco_+sulfate_+Met' = 'SD+Met+Glu w/\nInorganic Sulfates',
                               'Agarose_Home_+sulfate_-Met' = 'SD-Met+Glu w/\nInorganic Sulfates',
                               'Agarose_Home_-sulfate_-Met' = 'SD-Met+Glu w/o\nInorganic Sulfates'),
                    guide = 'none') +
  theme_minimal() +
  theme(plot.title = element_blank(),
        plot.margin = grid::unit(c(0,0,0,0), "mm"),
        panel.spacing = grid::unit(c(0,0,0,0), "mm"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank()) +
  coord_cartesian(xlim = c(0,910))


plot_grid(fig3d.expt, fig3d,
          ncol = 1, align = 'hv', axis = 'lr',
          rel_heights = c(0.3,1.2))  

