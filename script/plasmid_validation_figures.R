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
data$stage <- factor(data$stage, levels = c('WC','PS1'))

######
plot.rfit.box <- data %>%
  filter(hours %in% c(48,115)) %>%
  ggplot(aes(x = 1, y = relative_fitness)) +
  geom_boxplot(aes(fill = bio_rep), 
               outlier.shape = NA) +
  facet_grid(stage*plasmid_backbone ~ deletion1*deletion2*plasmid_orf,
             labeller = labeller(stage = c('WC' = 'YPDA', 'PS1' = 'SD-Met-Cys+Gal'))) +
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
  coord_cartesian(ylim = c(0.1,10))
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX.jpg",fig_path, expt.name),
       plot.rfit.box,
       height = two.c*2/3, width = two.c, units = 'mm',
       dpi = 600)  

