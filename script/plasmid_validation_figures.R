##### PLASMID VALIDATION FIGURES
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 08/10/2021 

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(ggridges)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(plotly)
library(scales)
library(reshape2)
library(locfit)
library(growthcurver)
library(rstatix)
library(gtools)
library(effsize)

out_path <- "~/R/Projects/methionine/data"
fig_path <- "~/R/Projects/methionine/figures"
res_path <- "~/R/Projects/methionine/results"

expt.name <- "plasmidval"

source("~/R/Projects/methionine/functions/colorstrip.R")
load(file = sprintf('%s/%s/colonysizes.RData', out_path, expt.name))

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

del1.levels <- c('met15','yll','met12')
del2.levels <- c('','yll','met12')
plab.levels <- c('2M','CEN')
plao.levels <- c('empty','yll','met12')

data$deletion1 <- factor(data$deletion1, levels = del1.levels)
data$deletion2 <- factor(data$deletion2, levels = del2.levels)
data$plasmid_backbone <- factor(data$plasmid_backbone, levels = plab.levels)
data$plasmid_orf <- factor(data$plasmid_orf, levels = plao.levels)
data$stage <- factor(data$stage, levels = c('WC','PS1','PS2','FS'))

######
data.sum %>%
  group_by(stage) %>%
  summarize(hours = max(hours))

######
plot.rfit.box <- data.sum %>%
  filter(hours %in% c(48,115,164)) %>%
  ggplot(aes(x = 1, y = relative_fitness)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(size = 1, aes(col = bio_rep, shape = expt_rep)) +
  facet_grid(stage*plasmid_backbone ~ deletion1*deletion2*plasmid_orf,
             labeller = labeller(stage = c('WC' = 'YPDA', 'PS1' = 'SD-Met-Cys+Gal',
                                           'PS2' = 'SD-Met-Cys+Gal', 'FS' = 'SD-Met-Cys+Gal'))) +
  labs(y = 'Relative Colony Size (log10)') +
  scale_y_continuous(trans = 'log10') +
  scale_fill_discrete(name = 'Biological Replicate') +
  theme_light() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.title.x = element_blank(),
        axis.text = element_text(size = txt),
        axis.text.x = element_blank(),
        # axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0.1,20))
# ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX.jpg",fig_path, expt.name),
#        plot.rfit.box,
#        height = two.c*2/3, width = two.c, units = 'mm',
#        dpi = 600)  

##### GROWTH CURVES
data.sum$cum_hrs <- NULL
data.sum$cum_hrs[data.sum$stage == 'WC'] <- data.sum$hours[data.sum$stage == 'WC']
data.sum$cum_hrs[data.sum$stage == 'PS1'] <- data.sum$hours[data.sum$stage == 'PS1'] +
  max(data.sum$cum_hrs[data.sum$stage == 'WC']) 
data.sum$cum_hrs[data.sum$stage == 'PS2'] <- data.sum$hours[data.sum$stage == 'PS2'] +
  max(data.sum$cum_hrs[data.sum$stage == 'PS1'])
data.sum$cum_hrs[data.sum$stage == 'FS'] <- data.sum$hours[data.sum$stage == 'FS'] +
  max(data.sum$cum_hrs[data.sum$stage == 'PS2'])

head(data.sum)
data.sum %>%
  group_by(stage) %>%
  summarise(hours = max(cum_hrs))

plot.rfit.tc <- data.sum %>%
  filter(stage != 'WC') %>%
  ggplot(aes(x = cum_hrs, y = cs)) +
  # stat_summary(data = data.sum %>% filter(stage == 'WC', cum_hrs == 48),
  #              fun = median, geom='point') +
  # geom_vline(xintercept = c(48,164,329,514), linetype = 'dashed', lwd = 0.4) +
  # geom_line(aes(col = bio_rep, linetype = expt_rep)) +
  stat_summary(fun = mean, geom="line", aes(col = stage)) +
  # scale_y_continuous(trans = 'pseudo_log') +
  facet_grid(plasmid_backbone ~ deletion1*deletion2*plasmid_orf) +
  labs(title = 'SD-Met-Cys+Gal',
       x = 'Time (hours)', y = 'Colony Size (pixels)') +
  scale_color_discrete(name = 'Stage') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_TIMECOURSE.jpg",fig_path, expt.name),
       plot.rfit.tc,
       height = one.5c, width = two.c*3/2, units = 'mm',
       dpi = 600)  


##### SUMMARY
data.sum2 <- data.sum %>%
  # filter(plasmid_backbone == '2M') %>%
  group_by(stage, hours, deletion1, deletion2, plasmid_backbone, plasmid_orf) %>%
  summarize(.groups = 'keep', fitness = median(relative_fitness, na.rm = T), cs = median(cs, na.rm = T)) %>%
  data.frame()
