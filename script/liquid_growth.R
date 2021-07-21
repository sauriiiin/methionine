##### LIQUID GROWTH EXPERIMENT
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 07/08/2021 

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

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")

out_path <- "~/R/Projects/methionine/data"
fig_path <- "~/R/Projects/methionine/figures"
res_path <- "~/R/Projects/methionine/results"

expt.name <- "liquidgrowth"

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

##### GATHER DATA
map <- read.csv(file = 'data/liquidgrowth/platemap.csv', stringsAsFactors = F)
map$wellID <- paste(map$Row_alpha, map$Column, sep = '')
map$Strain[map$Strain == 'FY4_met13-15'] <- 'BY4741_met3'
data <- read.csv(file = 'data/liquidgrowth/GrowDen.csv')

##### PLATE MAP
map$Starting.OD <- factor(map$Starting.OD, levels = c(0.05,0.1,0.3))
head(map)

map %>%
  group_by(Starting.OD, Strain, Replicate) %>%
  count() %>%
  data.frame()

map %>%
  ggplot(aes(x = Column, y = Row)) +
  geom_tile(aes(fill = Starting.OD), col = 'black') +
  geom_text(aes(label = Strain), size = 2) +
  scale_y_continuous(breaks = seq(1,8,1)) +
  scale_x_continuous(breaks = seq(1,12,1)) +
  coord_cartesian(xlim = c(1,12), ylim = c(8,1)) +
  theme_linedraw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank())

##### GROWTH CURVES
data$Time <- as.numeric(data$Time)
data$Time <- seq(0,15*(length(data$Time)-1),15)

gcdat <- melt(data, id.vars = c('Time'), variable.name = 'wellID', value.name = 'OD')
# gcdat <- merge(gcdat, map, by = 'wellID')

gc.res <- SummarizeGrowthByPlate(data)
gc.res <- merge(gc.res, map, by.x = 'sample', by.y = 'wellID')
head(gc.res)

gc.gr.res <- NULL
for (i in 2:dim(data)[2]) {
  if (colnames(data[i]) %in% map$wellID) {
    fit0 <- fit_easylinear(data$Time[2:dim(data)[1]][data[2:dim(data)[1],i] > 0],
                           data[2:dim(data)[1],i][data[2:dim(data)[1],i] > 0],
                           h = 20, quota = 1);
    
    temp_res <- data.frame(colnames(data[i]), maxgr = coef(fit0)[[3]],
                           dtime = log(2)/coef(fit0)[[3]], ltime = coef(fit0)[[4]])
    gc.gr.res <- rbind(gc.gr.res, temp_res)
  }
}
gc.gr.res <- data.frame(gc.gr.res)
colnames(gc.gr.res) <- c('sample','gr','dtime','lag')
head(gc.gr.res)

gc.res <- merge(gc.res, gc.gr.res, by = 'sample')
gc.res <- merge(gc.res, gcdat[gcdat$Time == 15,c(2,3)], by.x = 'sample', by.y = 'wellID')

gc.res %>%
  # filter(Strain == 'FY4') %>%
  ggplot(aes(x = Strain, y = OD)) +
  geom_violin() +
  geom_jitter()

gc.res$category[gc.res$Strain == 'BY4741' & gc.res$OD <= 0.006] <- 'Low'
gc.res$category[gc.res$Strain == 'BY4741' & gc.res$OD >= 0.03] <- 'High'
gc.res$category[gc.res$Strain == 'BY4741' & is.na(gc.res$category)] <- 'Medium'

gc.res$category[gc.res$Strain == 'BY4741_met3' & gc.res$sample %in% c('B11','F2','F7')] <- 'Low'
gc.res$category[gc.res$Strain == 'BY4741_met3' & gc.res$sample %in% c('B6','C10','D4')] <- 'High'
gc.res$category[gc.res$Strain == 'BY4741_met3' & is.na(gc.res$category)] <- 'Medium'

gc.res$category[gc.res$Strain == 'FY4' & gc.res$OD <= 0.001] <- 'Low'
gc.res$category[gc.res$Strain == 'FY4' & gc.res$OD >= 0.01] <- 'High'
gc.res$category[gc.res$Strain == 'FY4' & is.na(gc.res$category)] <- 'Medium'

gc.res$category[gc.res$Strain == 'FY4_met15' & gc.res$OD <= 0.004] <- 'Low'
gc.res$category[gc.res$Strain == 'FY4_met15' & gc.res$OD >= 0.016] <- 'High' 
gc.res$category[gc.res$Strain == 'FY4_met15' & is.na(gc.res$category)] <- 'Medium'

gc.res$category[gc.res$Strain == 'FY4_met3' & gc.res$OD >= 0.009] <- 'High'
gc.res$category[gc.res$Strain == 'FY4_met3' & gc.res$OD <= 0.003] <- 'Low'
gc.res$category[gc.res$Strain == 'FY4_met3' & is.na(gc.res$category)] <- 'Medium'

gc.res$category <- factor(gc.res$category, levels = c('Low','Medium','High'))

gc.res$outlier <- FALSE
for (s in unique(gc.res$Strain)) {
  for (c in unique(gc.res$category[gc.res$Strain == s])) {
    gc.res$outlier[gc.res$Strain == s & gc.res$category == c] <- isoutlier(gc.res$auc_e[gc.res$Strain == s & gc.res$category == c],2) | 
                        isoutlier(gc.res$gr[gc.res$Strain == s & gc.res$category == c],2)
  }
}



gc.res.sum <- gc.res %>% 
  filter(outlier == FALSE, sample %in% gcdat$wellID[gcdat$Time == max(gcdat$Time) & gcdat$OD >= 0.1]) %>%
  group_by(Strain, category) %>%
  summarise(gr = mean(gr, na.rm = T), auc = mean(auc_e, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

plot.gc <- merge(gcdat, gc.res, by.x = 'wellID', by.y = 'sample') %>%
  filter(wellID != 'B3') %>%
  ggplot(aes(x = Time, y = OD.x)) +
  # geom_point(size = 1) +
  geom_smooth(method = 'loess') +
  labs(x = 'Time (minutes)', y = 'OD') +
  geom_text(data = gc.res.sum, aes(x = 4000, y = 1.35,
                                   label = sprintf('GR = %0.3f | AUC = %0.2f', gr, auc)),
            size = 2) +
  facet_wrap(.~Strain*category, ncol = 3) +
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
  coord_cartesian(ylim = c(0, 1.4))
ggsave(sprintf("%s/%s/GROWTH_CURVES.jpg",fig_path, expt.name), plot.gc,
       height = 200, width = 150, units = 'mm',
       dpi = 600)


