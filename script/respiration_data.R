##### RESPIRATION EXPERIMENT - DATA
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 05/03/2021 

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
library(RMariaDB)

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- "~/R/Projects/methionine/data"
fig_path <- "~/R/Projects/methionine/figures"

##### GATHER AND CLEAN DATA FROM SQL

tables <- rbind(c('respiration_exp1R_MiM_Et_Rep2_JX__BinaThreZoom_cleaned_384_JPEG', 'MET-', 'URA+', 'EtOH', 'SD-MET+EtOH'),
  c('respiration_exp1R_MiM_Glu_Rep2_JX__BinaThreZoom_cleaned_384_JPEG', 'MET-', 'URA+', 'Glucose', 'SD-MET+Glucose'),
  c('respiration_exp1R_PlM_Et_Rep2_JX__BinaThreZoom_cleaned_384_JPEG', 'MET+', 'URA+', 'EtOH', 'SD+MET+EtOH'),
  c('respiration_exp1R_PlM_Glu_Rep2_JX__BinaThreZoom_cleaned_384_JPEG', 'MET+', 'URA+', 'Glucose', 'SD+MET+Glucose'),
  c('respiration_exp1R_Ura_Rep2_JX__binary_Zoom_cleaned_384_JPEG', 'MET+', 'URA-', 'Glucose', 'SD+MET-URA+Glucose'))

data <- NULL
for (i in seq(1,dim(tables)[1])) {
  temp <- dbGetQuery(conn, sprintf('select a.pos, a.hours, a.average, b.density, b.plate, b.row, b.col, c.orf_name
                                          from Branden.%s a,
                                   Branden.respiration_exp_1_JX_pos2coor b,
                                   Branden.respiration_exp_1_JX_pos2orf_name c
                                   where a.pos = b.pos and b.pos = c.pos
                                   order by a.hours, b.plate, b.col, b.row',
                                   tables[i,1]))
  for (h in unique(temp$hours)) {
    for (o in unique(temp$orf_name)) {
      temp$average[temp$orf_name == o & temp$hours == h][isoutlier(temp$average[temp$orf_name == o & temp$hours == h], 2)] <- NA
    }
  }

  temp$methionine <- tables[i,2]
  temp$uracil <- tables[i,3]
  temp$carbon <- tables[i,4]
  temp$condition <- tables[i,5]
  
  data <- rbind(data,temp)
}
data <- data.frame(data)
head(data)

data$bio_rep[data$row%%2==1 & data$col%%2==1] = '1'
data$bio_rep[data$row%%2==0 & data$col%%2==1] = '3'
data$bio_rep[data$row%%2==1 & data$col%%2==0] = '2'
data$bio_rep[data$row%%2==0 & data$col%%2==0] = '4'

##### SAVE COLONY SIZE DATA
save(data, file = sprintf("%s/methionine_respiration_colonysizes.RData", out_path))


##### DATA FOR GROWTH CURVE ANALYSIS
data <- data[!(data$hours %in% c(46.23, 87.35, 25.88, 16.75, 54.36, 37.13)) &
               data$orf_name != 'BOR',]

data.gc <- data %>%
  group_by(condition, orf_name, bio_rep, hours) %>%
  summarise(average = median(average, na.rm = T)) %>%
  data.frame()

data.pred <- NULL
col.names <- NULL
for (o in unique(data.gc$orf_name)) {
  for (c in unique(data.gc$condition[data.gc$orf_name == o])) {
    for (b in unique(data.gc$bio_rep[data.gc$orf_name == o & data.gc$condition == c])) {
      temp <- data.gc[data.gc$orf_name == o & data.gc$condition == c & data.gc$bio_rep == b,]
      # if (sum(log(temp$average) < 5) <= 5) {
        lo <- loess.smooth(temp$hours, log(temp$average),
                           span = 0.6, evaluation = 500, degree = 2,
                           family = 'gaussian')
        data.pred <- cbind(data.pred,exp(lo$y))
        col.names <- cbind(col.names,paste(o,c,b,sep = ','))
        
        temp.plot <- ggplot() +
          geom_line(data = data.frame(lo), aes(x = x, y = y)) +
          geom_point(data = temp, aes(x = hours, y = log(average)), shape = 1) +
          labs(y = 'Log( Colony Size (pixels) )',
               x = 'Time (hours)',
               title = paste(o,c,b,sep = ' | ')) +
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
          coord_cartesian(ylim = c(3,9))
        ggsave(sprintf("%s/growthcurves/%s_%s_%s.jpg",fig_path, o, c, b), temp.plot,
               height = 100, width = 100, units = 'mm',
               dpi = 600)
      # }
    }
  }
}
data.pred <- cbind(lo$x, data.pred)
data.pred <- data.frame(data.pred)
colnames(data.pred) <- c('Time',col.names)
head(data.pred)

##### SAVE GROWTH CURVE DATA
save(data.pred, file = sprintf("%s/methionine_respiration_growthcurve.RData", out_path))

#####
# END
#####



