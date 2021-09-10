
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
library(ggforce)
library(plotly)
library(scales)
library(reshape2)

data <- read.csv(file = 'data/h2s/pilot/210909_nocell_pilot.csv', stringsAsFactors = F)

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
data <- melt(data, id.vars = c('strain','condition','sOD','time'),
             variable.name = 'wavelength', value.name = 'OD')

data %>%
  ggplot(aes(x = time, y = OD)) +
  geom_point() +
  geom_line(aes(col = strain, linetype = condition)) +
  facet_grid(wavelength~sOD) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = txt),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = ggtext::element_markdown(size = txt,
                                              margin = margin(0.1,0,0.1,0, "mm")))
  