
library(RMariaDB)
library(readxl)
library(stringr)
library(dplyr)

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")


tables <- data.frame(
  expt_id = c('SCmLeu', 'SCmLeu', 'SDmLeu', 'SDmLeu', 'YPDA', 'SDpMET', 'SDpMET'),
  arm = c('R1', 'R2', 'R1', 'R2', 'R1', 'R1', 'R2')
)

##### GATHER AND CLEAN DATA FROM SQL
data.cbn.leu <- NULL
for (e in unique(tables$expt_id)) {
  for (a in unique(tables$arm[tables$expt_id == e])) {
    temp <- dbGetQuery(conn, sprintf('select a.pos, a.hours, a.average, b.density, b.plate, b.row, b.col, c.orf_name
                                       from CARS_%s_%s_384_CLEAN a, CARS_mL_pos2coor b, CARS_mL_pos2orf_name c
                                       where a.pos = b.pos and b.pos = c.pos
                                       order by a.hours, b.plate, b.col, b.row', e, a))
    temp$bio_rep[temp$row%%2==1 & temp$col%%2==1] = '1'
    temp$bio_rep[temp$row%%2==0 & temp$col%%2==1] = '3'
    temp$bio_rep[temp$row%%2==1 & temp$col%%2==0] = '2'
    temp$bio_rep[temp$row%%2==0 & temp$col%%2==0] = '4'
    
    for (h in unique(temp$hours)) {
      for (o in unique(temp$orf_name)) {
        for (b in unique(temp$bio_rep)) {
          temp$average[temp$orf_name == o & temp$hours == h & temp$bio_rep == b][isoutlier(temp$average[temp$orf_name == o & temp$hours == h & temp$bio_rep == b], 2)] <- NA
        }
      }
      temp$relative_fitness[temp$hours == h] <- temp$average[temp$hours == h]/median(temp$average[temp$hours == h & 
                                                                                                    temp$orf_name == 'FY4'], na.rm = T)
    }
    temp$condition <- e
    temp$expt_rep <- a
    
    data.cbn.leu <- rbind(data.cbn.leu,temp)
  }
}
data.cbn.leu <- data.frame(data.cbn.leu)
data.cbn.leu$carbon <- 'Glucose'
head(data.cbn.leu)


data.cbn.leu %>%
  ggplot(aes(x = arm, y = relative_fitness)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(.~orf_name)


data.cbn.leu %>%
  ggplot(aes(x = col, y = row)) +
  geom_tile(aes(fill = average)) +
  scale_y_reverse() +
  facet_wrap(.~arm*expt_rep*orf_name)


####
data.cbn.leu$base[str_detect(data.cbn.leu$condition, 'SC')] <- 'SC'
data.cbn.leu$base[str_detect(data.cbn.leu$condition, 'SD')] <- 'SD'
data.cbn.leu$base[is.na(data.cbn.leu$base)] <- 'YPDA'

head(data.cbn)
head(data.cbn.leu)





