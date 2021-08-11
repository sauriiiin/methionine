##### NO SULFATE FIGURES
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 08/11/2021 

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

expt.name <- "nosulfate"

source("~/R/Projects/methionine/functions/colorstrip.R")
load(file = sprintf('%s/%s/colonysizes.RData', out_path, expt.name))
load(file = sprintf('%s/%s/stats.RData', out_path, expt.name))

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

#####
head(data)
# temp <- cbind(str_split(unique(data$arm), '_', simplify = T), unique(data$arm))
temp <- data.frame(cbind(str_split(unique(data$condition), '_', simplify = T), unique(data$condition)), stringsAsFactors = F)
colnames(temp) <- c('ynb_type','sulfate','media','base','condition')
temp[temp$condition == 'SD+Met+Glu_Agar',] <- c('Difco','+sulfate','SD+Met','Agar','SD+Met+Glu_Agar')
temp$sulfate[temp$sulfate == 'YNB'] <- '+sulfate'
temp$sulfate <- str_replace(temp$sulfate, 'sulfates', 'sulfate')
temp$methionine <- str_remove(temp$media, 'SD')
data <- merge(data, temp, by = 'condition')
data.sum <- merge(data.sum, temp, by = 'condition')

data$stage[data$stage == 'FS'] <- 'S1'
data$hours[data$stage == 'S1'] <- 149
data.sum$stage[data.sum$stage == 'FS'] <- 'S1'
data.sum$hours[data.sum$stage == 'S1'] <- 149

data$stage <- factor(data$stage, levels = c('S1','S2','S3','Re1'))
data$ynb_type <- factor(data$ynb_type, levels = c('Difco','Home'))
data$sulfate <- factor(data$sulfate, levels = c('+sulfate','-sulfate'))
data$base <- factor(data$base, levels = c('Agar','Agarose'))
data$methionine <- factor(data$methionine, levels = c('+Met','-Met'))
data$orf_name <- factor(data$orf_name, levels = c('FY4','FY4_met3del','FY4_met15del','BY4742','BY4741'))

data.sum$stage <- factor(data.sum$stage, levels = c('S1','S2','S3','Re1'))
data.sum$ynb_type <- factor(data.sum$ynb_type, levels = c('Difco','Home'))
data.sum$sulfate <- factor(data.sum$sulfate, levels = c('+sulfate','-sulfate'))
data.sum$base <- factor(data.sum$base, levels = c('Agar','Agarose'))
data.sum$methionine <- factor(data.sum$methionine, levels = c('+Met','-Met'))
data.sum$orf_name <- factor(data.sum$orf_name, levels = c('FY4','FY4_met3del','FY4_met15del','BY4742','BY4741'))

strain.labs <- c('FY4','FY4-*met3Δ*','FY4-*met15Δ*','BY4742','BY4741')
######
head(data)
plot.cs.box <- data %>%
  filter(hours %in% c(92,149)) %>%
  ggplot(aes(x = orf_name, y = average)) +
  geom_boxplot(aes(fill = orf_name, group = orf_name), size = 0.2,
               outlier.shape = NA) +
  facet_grid(base*ynb_type*sulfate*methionine ~ stage,
             labeller = labeller(stage = c('S1' = 'Pin #1',
                                           'S2' = 'Pin #2',
                                           'S3' = 'Pin #3',
                                           'Re1' = 'Pin #3 - Revival'))) +
  labs(x = 'Strain', y = 'Colony Size') +
  scale_x_discrete(labels = strain.labs) +
  # scale_y_continuous(trans = 'log2') +
  scale_fill_discrete(name = '', guide = F) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX.jpg",fig_path, expt.name),
       plot.cs.box,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)

#####
data %>%
  # filter(stage == 'S3') %>%
  ggplot(aes(x = hours, y = average)) +
  stat_summary(aes(col = orf_name, group = orf_name), fun = 'mean', geom = 'line') +
  facet_grid(base*ynb_type*sulfate*methionine ~ stage)

