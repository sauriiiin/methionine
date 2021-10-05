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
# load(file = sprintf('%s/%s/stats.RData', out_path, expt.name))

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
data.sum$stage[data.sum$stage == 'FS'] <- 'S1'

data$stage <- factor(data$stage, levels = c('S1','S2','S3','S4','S5','Re1','Re2'))
data$ynb_type <- factor(data$ynb_type, levels = c('Difco','Home'))
data$sulfate <- factor(data$sulfate, levels = c('+sulfate','-sulfate'))
data$base <- factor(data$base, levels = c('Agar','Agarose'))
data$methionine <- factor(data$methionine, levels = c('+Met','-Met'))
data$orf_name <- factor(data$orf_name, levels = c('FY4','FY4_met3del','FY4_met15del','BY4742','BY4741'))

data.sum$stage <- factor(data.sum$stage, levels = c('S1','S2','S3','S4','S5','Re1','Re2'))
data.sum$ynb_type <- factor(data.sum$ynb_type, levels = c('Difco','Home'))
data.sum$sulfate <- factor(data.sum$sulfate, levels = c('+sulfate','-sulfate'))
data.sum$base <- factor(data.sum$base, levels = c('Agar','Agarose'))
data.sum$methionine <- factor(data.sum$methionine, levels = c('+Met','-Met'))
data.sum$orf_name <- factor(data.sum$orf_name, levels = c('FY4','FY4_met3del','FY4_met15del','BY4742','BY4741'))

strain.labs <- c('FY4','FY4-*met3Δ*','FY4-*met15Δ*','BY4742','BY4741')
######
# head(data)
# plot.cs.box <- data %>%
#   filter(hours %in% c(92,149)) %>%
#   ggplot(aes(x = orf_name, y = average)) +
#   geom_boxplot(aes(fill = orf_name, group = orf_name), size = 0.2,
#                outlier.shape = NA) +
#   facet_grid(base*ynb_type*sulfate*methionine ~ stage,
#              labeller = labeller(stage = c('S1' = 'Pin #1',
#                                            'S2' = 'Pin #2',
#                                            'S3' = 'Pin #3',
#                                            'Re1' = 'Pin #3 - Revival'))) +
#   labs(x = 'Strain', y = 'Colony Size') +
#   scale_x_discrete(labels = strain.labs) +
#   # scale_y_continuous(trans = 'log2') +
#   scale_fill_discrete(name = '', guide = F) +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
#         axis.ticks.x = element_blank(),
#         legend.title = element_text(size = titles),
#         legend.text = element_text(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm")))
# ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX.jpg",fig_path, expt.name),
#        plot.cs.box,
#        height = two.c, width = two.c, units = 'mm',
#        dpi = 600)

#####
data %>%
  filter(orf_name %in% c('FY4','BY4742'), base == 'Agarose') %>%
  group_by(stage, base, ynb_type, sulfate, methionine, orf_name, hours) %>%
  summarize(average = median(average, na.rm = T), .groups = 'keep') %>%
  ggplot(aes(x = hours, y = average)) +
  geom_line(aes(col = orf_name, linetype = stage)) +
  geom_point(aes(col = orf_name, shape = stage)) +
  facet_wrap(.~base*ynb_type*sulfate*methionine)

data$average[is.na(data$average)] <- 0
# data.tc <- data %>%
#   # filter(orf_name %in% c('FY4','BY4742'), base == 'Agarose') %>%
#   group_by(stage, base, ynb_type, sulfate, methionine, orf_name, hours) %>%
#   summarize(average = median(average, na.rm = T), .groups = 'keep') %>%
#   data.frame()

data$cum_hrs <- NULL
data$cum_hrs[data$stage == 'S1'] <- data$hours[data$stage == 'S1']
data$cum_hrs[data$stage == 'S2'] <- data$hours[data$stage == 'S2'] +
  max(data$cum_hrs[data$stage == 'S1'])
data$cum_hrs[data$stage == 'S3'] <- data$hours[data$stage == 'S3'] +
  max(data$cum_hrs[data$stage == 'S2'])
data$cum_hrs[data$stage == 'Re1'] <- data$hours[data$stage == 'Re1'] +
  max(data$cum_hrs[data$stage == 'S2'])
data$cum_hrs[data$stage == 'S4'] <- data$hours[data$stage == 'S4'] +
  max(data$cum_hrs[data$stage == 'S3'])
data$cum_hrs[data$stage == 'S5'] <- data$hours[data$stage == 'S5'] +
  max(data$cum_hrs[data$stage == 'S4'])
data$cum_hrs[data$stage == 'Re2'] <- data$hours[data$stage == 'Re2'] +
  max(data$cum_hrs[data$stage == 'S5'])

data$id <- paste(data$base, data$ynb_type, data$sulfate, data$methionine, sep = '_')
data$stage <- as.character(data$stage)

data$stage[data$stage == 'Re1'] <- 'S3'
data$stage[data$stage == 'Re2'] <- 'S6'

data$expt_rep[str_detect(data$arm, 'R1')] <- 'R1'
data$expt_rep[str_detect(data$arm, 'R2')] <- 'R2'

# plot.cs.tc <- 
data[!(data$id == 'Agar_Difco_+sulfate_+Met' & data$stage == 'S1'),] %>%
  filter(id %in% c('Agar_Difco_+sulfate_+Met',
                   'Agarose_Home_+sulfate_-Met',
                   'Agarose_Home_-sulfate_-Met')) %>%
  group_by(stage, base, ynb_type, sulfate, methionine, orf_name, hours, cum_hrs, id) %>%
  summarize(average = median(average, na.rm = T), .groups = 'keep') %>%
  data.frame() %>%
  ggplot(aes(x = cum_hrs, y = average)) +
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=.2),
  #              aes(group = id, fill = id), geom="ribbon", alpha = 0.4) +
  # stat_summary(aes(group = id, col = id), fun=mean, geom="line", lwd =0.7) +
  geom_line(aes(col = id)) +
  # geom_point(aes(col = orf_name, shape = stage)) +
  # facet_wrap(.~base*ynb_type*sulfate*methionine*stage) +
  scale_x_continuous(breaks = seq(0,1000,50)) +
  facet_wrap(.~orf_name, nrow = 5) +
  labs(x = 'Time (hours)',
       y = 'Colony Size (pixels)') +
  scale_color_discrete(name = 'Strain') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) #+
  # coord_cartesian(xlim = c(0,750))
# ggsave(sprintf("%s/%s/COLONY_SIZE_TIMECOURSE.jpg",fig_path, expt.name),
#        plot.cs.tc,
#        height = two.c, width = two.c, units = 'mm',
#        dpi = 600)


data[data$stage == 'S1' & data$hours == 165 & data$average != 0,] %>%
  ggplot(aes(x = orf_name, y = average, fill = base, alpha = ynb_type)) +
  geom_boxplot(size = 0.2, outlier.shape = NA,
               position = position_dodge(width = 0.6)) +
  facet_wrap(.~methionine*sulfate, nrow = 3) +
  scale_alpha_manual(name = 'YNB Type',
                     values = c('Difco' = 0.4,
                                'Home' = 1))


