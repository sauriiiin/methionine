
##### INITIALIZE
source('/home/sbp29/R/Projects/methionine/paper/scripts/initialize.R')

table <- rbind(c('MM','R1','SD-Met-Cys+Glu'),
               c('MM','R2','SD-Met-Cys+Glu'),
               c('MMEDTA','R1','SD-Met-Cys+Glu+FeEDTA'),
               c('PM','R1','SD+Met-Cys+Glu'),
               c('PM','R2','SD+Met-Cys+Glu'),
               c('PMEDTA','R1','SD+Met-Cys+Glu+FeEDTA')) %>%
  data.frame()

colnames(table) <- c('id','expt_rep','condition')

data.pv3 <- NULL
for (i in seq(1,dim(table)[1])) {
  temp <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate_no, b.plate_row, b.plate_col, c.orf_name
                                  from PV3_FS_%s_%s_384_CLEAN a, PV3_pos2coor b, PV3_pos2orf_name c
                                  where a.pos = b.pos and b.pos = c.pos
                                  order by a.hours, b.plate_no, b.plate_col, b.plate_row',
                                   table$id[i],table$expt_rep[i]))
  
  temp$bio_rep[temp$plate_row%%2==1 & temp$plate_col%%2==1] = '1'
  temp$bio_rep[temp$plate_row%%2==0 & temp$plate_col%%2==1] = '3'
  temp$bio_rep[temp$plate_row%%2==1 & temp$plate_col%%2==0] = '2'
  temp$bio_rep[temp$plate_row%%2==0 & temp$plate_col%%2==0] = '4'
  
  for (h in unique(temp$hours)) {
    for (o in unique(temp$orf_name)) {
      for (b in unique(temp$bio_rep)) {
        temp$average[temp$orf_name == o & temp$hours == h & temp$bio_rep == b][isoutlier(temp$average[temp$orf_name == o & temp$hours == h & temp$bio_rep == b], 2)] <- NA
      }
    }
    temp$relative_fitness[temp$hours == h] <- temp$average[temp$hours == h]/median(temp$average[temp$hours == h & 
                                                                                                  temp$orf_name %in% c('FY4_empty')], na.rm = T)
  }
  
  temp$id <- table$id[i]
  temp$expt_rep <- table$expt_rep[i]
  temp$condition <- table$condition[i]
  data.pv3 <- rbind(data.pv3, temp)
}
head(data.pv3)

data.pv3$average[data.pv3$plate_no == 2 & data.pv3$bio_rep == '4'] <- NA
data.pv3$average[data.pv3$plate_no == 8 & data.pv3$plate_col == 15 & data.pv3$plate_row == 6] <- NA
data.pv3$average[data.pv3$id == 'MMEDTA' & data.pv3$expt_rep == 'R1' & data.pv3$plate_col > 8] <- NA

data.pv3$relative_fitness[is.na(data.pv3$average)] <- NA

data.pv3$orf_name <- factor(data.pv3$orf_name,
                            levels = c("FY4_empty","FY4_ylldel_empty","FY4_met15del_empty","FY4_met15ylldel_empty",
                                       "FY4_ylldel_yll","FY4_met15del_yll","FY4_met15ylldel_yll",
                                       "FY4_met15del_yll_k376a","FY4_met15ylldel_yll_k376a"))

data.pv3 %>%
  filter(hours == 120, id %in% c('MM','PM')) %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  geom_boxplot(aes(fill = orf_name),
               outlier.shape = NA) +
  scale_fill_discrete(guide = 'none') +
  scale_x_discrete(labels = c('FY4_empty' = 'FY4<br />w/ Empty Plasmid',
                              'FY4_ylldel_empty' = 'FY4-*yll058wΔ*<br />w/ Empty Plasmid',
                              'FY4_met15del_empty' = 'FY4-*met15Δ*<br />w/ Empty Plasmid',
                              'FY4_met15ylldel_empty' = 'FY4-*met15Δ yll058wΔ*<br />w/ Empty Plasmid',
                              'FY4_ylldel_yll' = 'FY4-*yll058wΔ*<br />w/ YLL058W Plasmid',
                              'FY4_met15del_yll' = 'FY4-*met15Δ*<br />w/ YLL058W Plasmid',
                              'FY4_met15ylldel_yll' = 'FY4-*met15Δ yll058wΔ*<br />w/ YLL058W Plasmid',
                              'FY4_met15del_yll_k376a' = 'FY4-*met15Δ*<br />w/ YLL058W (K376A) Plasmid',
                              'FY4_met15ylldel_yll_k376a' = 'FY4-*met15Δ yll058wΔ*<br />w/ YLL058W (K376A) Plasmid')) +
  labs(y = 'fitness', x = 'strain') +
  facet_wrap(.~condition,
             nrow = 2, ncol = 2) +
  coord_flip(ylim = c(0,1.8)) +
  theme_linedraw() +
  theme(axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.y = element_blank(),
        axis.text.y = ggtext::element_markdown(size = txt),
        strip.text = element_text(size = txt))


fig3d <- data.pv3 %>%
  filter(hours == 120, id %in% c('MM')) %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  geom_boxplot(outlier.shape = NA, fill = '#9E9E9E', size = 0.3) +
  scale_x_discrete(limits = c(
    'FY4_empty',
    # 'FY4_ylldel_empty',
    'FY4_met15del_empty',
    'FY4_met15ylldel_empty',
    # 'FY4_ylldel_yll',
    'FY4_met15del_yll',
    'FY4_met15ylldel_yll',
    'FY4_met15del_yll_k376a',
    'FY4_met15ylldel_yll_k376a'),
    labels = c('FY4_empty' = 'FY4<br />w/ Empty<br />Plasmid',
               'FY4_ylldel_empty' = 'FY4-*yll058wΔ*<br />w/ Empty<br />Plasmid',
               'FY4_met15del_empty' = 'FY4-*met15Δ*<br />w/ Empty<br />Plasmid',
               'FY4_met15ylldel_empty' = 'FY4-*met15Δ yll058wΔ*<br />w/ Empty<br />Plasmid',
               'FY4_ylldel_yll' = 'FY4-*yll058wΔ*<br />w/ YLL058W<br />Plasmid',
               'FY4_met15del_yll' = 'FY4-*met15Δ*<br />w/ YLL058W<br />Plasmid',
               'FY4_met15ylldel_yll' = 'FY4-*met15Δ yll058wΔ*<br />w/ YLL058W<br />Plasmid',
               'FY4_met15del_yll_k376a' = 'FY4-*met15Δ*<br />w/ YLL058W<br />(K376A) Plasmid',
               'FY4_met15ylldel_yll_k376a' = 'FY4-*met15Δ yll058wΔ*<br />w/ YLL058W<br />(K376A) Plasmid'),
    position = "top") + #,
    # guide = guide_axis(n.dodge=2)) +
  labs(y = 'Relative Colony Size in\nSD-Met-Cys+Gal+G418') +
  theme_linedraw() +
  theme(plot.margin = margin(0,3,10,0),
        plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        panel.background = element_rect(fill = pinningpanelbg),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x.top = ggtext::element_markdown(size = txt-1),
        axis.ticks.x.top = element_line(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.2))


### STATS
data.pv2.temp <- data.pv3 %>%
  filter(orf_name %in% c('FY4_empty',
                         # 'FY4_ylldel_empty',
                         'FY4_met15del_empty',
                         'FY4_met15ylldel_empty',
                         # 'FY4_ylldel_yll',
                         'FY4_met15del_yll',
                         'FY4_met15ylldel_yll',
                         'FY4_met15del_yll_k376a',
                         'FY4_met15ylldel_yll_k376a'),
         id == 'MM', hours == 120) %>%
  group_by(orf_name, expt_rep, bio_rep) %>%
  summarise(relative_fitness = median(relative_fitness, na.rm = T), .groups = 'keep') %>%
  data.frame()

data.pv2.stats <- NULL
for (o1 in unique(data.pv2.temp$orf_name)) {
  for (o2 in unique(data.pv2.temp$orf_name)) {
    if (o1 != o2) {
      temp.stat <- compare_means(relative_fitness ~ orf_name,
                                 data.pv2.temp %>%
                                   filter(orf_name %in% c(o1, o2)) %>%
                                   data.frame(),
                                 method = "kruskal.test", paired = FALSE)
      temp.es <- (median(data.pv2.temp$relative_fitness[data.pv2.temp$orf_name == o1], na.rm = T) - 
                    median(data.pv2.temp$relative_fitness[data.pv2.temp$orf_name == o2], na.rm = T))/
        median(data.pv2.temp$relative_fitness[data.pv2.temp$orf_name == o2], na.rm = T)
      data.pv2.stats <- rbind(data.pv2.stats,
                              data.frame(expt = 'plasmid_validation2', 
                                         ref_orf = o2,
                                         mut_orf = o1,
                                         effect_size = temp.es, p = temp.stat$p))
    }
  }
}
data.pv2.stats$p.adj <- adjust_pvalue(data.pv2.stats$p, method = 'BH')
write.csv(data.pv2.stats, file = 'results/final/stats/plasmid_validation3.csv',row.names = F)



