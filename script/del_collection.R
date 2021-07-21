##### DELETION SCREEN
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 07/13/2021 

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
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
library(locfit)
library(growthrates)
library(RMariaDB)

out_path <- "~/R/Projects/methionine/data"
fig_path <- "~/R/Projects/methionine/figures"
res_path <- "~/R/Projects/methionine/results"

expt.name <- "deletion"

source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

del_mm_data <- dbGetQuery(conn, 'select a.*, b.density, b.plate, b.row, b.col
        from MET_DEL_PS2_MM_1536_FITNESS a, MET_DEL_pos2coor b
                          where a.pos = b.pos
                          order by a.hours, b.density, b.plate, b.col, b.row')

del_mm_stats <- dbGetQuery(conn, 'select a.* from MET_DEL_PS2_MM_1536_FITNESS_STATS a')


# del_mm_sig <- dbGetQuery(conn, 'select orf_name from MET_DEL_PS2_MM_1536_PVALUE
# where p <= 0.05 and hours = 69')
# 
# del_mm_data$sig[del_mm_data$orf_name %in% del_mm_sig$orf_name] <- 'Yes'
# del_mm_data$sig[is.na(del_mm_data$sig)] <- 'No'
# 
# del_mm_data %>%
#   filter(hours == max(del_mm_data$hours)) %>%
#   ggplot(aes(x = col, y = row)) +
#   geom_point(aes(size = average, col = sig)) +
#   scale_size_continuous(range = c(0, 2)) +
#   facet_wrap(.~plate)

del_mm_summ <- del_mm_data[!(del_mm_data$plate %in% c(1,2,3,4,5,6,10,13)), ] %>%
# del_mm_summ <- del_mm_data %>%
  group_by(hours, orf_name) %>%
  summarize(cs = median(average, na.rm = T), .groups = 'keep') %>%
  data.frame()
del_mm_summ <- merge(del_mm_summ, del_mm_stats, by = c('hours','orf_name'))

del_mm_lims <- del_mm_data %>%
  filter(orf_name == 'BY4741') %>%
  group_by(hours) %>%
  summarize(cs_ll = quantile(average, 0.025, na.rm = T),
            cs_m = median(average, 0.5, na.rm = T),
            cs_ul = median(average, 0.975, na.rm = T),
            f_ll = quantile(fitness, 0.025, na.rm = T),
            f_m = median(fitness, 0.5, na.rm = T),
            f_ul = median(fitness, 0.975, na.rm = T)) %>%
  data.frame()


del_mm_summ[del_mm_summ$orf_name %in% unique(del_mm_summ$orf_name[del_mm_summ$hours ==  max(del_mm_summ$hours) &
                                                                    del_mm_summ$cs >= 500]) &
              !is.na(del_mm_summ$orf_name),] %>%
  filter(hours > 0) %>%
  ggplot(aes(x = hours, y = cs_mean)) +
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
  #              aes(group = orf_name, fill = orf_name), geom="ribbon", alpha = 0.4) +
  # stat_summary(aes(group = orf_name, col = orf_name), fun=mean, geom="line", lwd =0.7)
  geom_line(aes(group = orf_name, col = orf_name), lwd =0.7) + 
  geom_line(data = del_mm_lims, aes(x = hours, y = f_ll), col = 'red', linetype = 'dashed', lwd = 0.5) +
  geom_line(data = del_mm_lims, aes(x = hours, y = f_ul), col = 'red', linetype = 'dashed', lwd = 0.5) +
  scale_color_discrete(guide = F)

del_mm_data[del_mm_data$orf_name %in% unique(del_mm_summ$orf_name[del_mm_summ$fitness > 100 & !is.na(del_mm_summ$fitness)]) & del_mm_data$hours == 68,]
del_mm_data[del_mm_data$orf_name == 'YNL275W' & del_mm_data$hours == 68,]

####
del_mm_data$orf_type[del_mm_data$orf_name == 'BY4741'] <- 'Reference'
del_mm_data$orf_type[del_mm_data$orf_name != 'BY4741'] <- 'Mutants'

del_mm_data[!is.na(del_mm_data$orf_type),] %>%
  filter(hours >= 60) %>%
  ggplot(aes(x = fitness)) +
  geom_line(aes(col = orf_type), stat = 'density') +
  facet_wrap(.~hours, scales = 'free') +
  coord_cartesian(xlim = c(0,5))


