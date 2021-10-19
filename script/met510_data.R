##### JAKES MUTANATS EXPERIMENT - data.510
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 05/05/2021 

##### INITIALIZE
library(RMariaDB)

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### GATHER AND CLEAN DATA FROM SQL
tables <- rbind(c('MM', 'MET510_pos2coor', 'MET510_pos2orf_name', 'attempt 1', 'SD-Met'),
                c('PM', 'MET510_pos2coor', 'MET510_pos2orf_name', 'attempt 1', 'SD-Ura'), # PM and MU plates were mislabled in the INFO file
                c('MU', 'MET510_pos2coor', 'MET510_pos2orf_name', 'attempt 1', 'SD+Met'),
                c('YPD', 'MET510_YPD_pos2coor', 'MET510_YPD_pos2orf_name', 'attempt 1', 'YPDA'))

data.510 <- NULL
for (i in seq(1,dim(tables)[1])) {
  temp <- dbGetQuery(conn, sprintf('select a.pos, a.hours, a.average, b.density, b.plate, b.row, b.col, c.orf_name
                                   from MET510_FS_%s_384_CLEAN a, %s b, %s c
                                   where a.pos = b.pos and b.pos = c.pos
                                   order by a.hours, b.plate, b.col, b.row',
                                   tables[i,1], tables[i,2], tables[i,3]))
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
    temp$relative_fitness[temp$hours == h] <- temp$average[temp$hours == h]/median(temp$average[temp$hours == h & temp$orf_name == 'FY4'], na.rm = T)
  }
  temp$attempt <- tables[i,4]
  temp$condition <- tables[i,5]
  
  data.510 <- rbind(data.510,temp)
}
data.510 <- data.frame(data.510)
head(data.510)

##### LABELS
strain.labs.510 <- data.frame(orf_name = c('FY4','met15del','met5del','met10del','BY4741_met5del','BY4741_met10del'),
                              labels = c('FY4', 'FY4-*met15Δ*', 'FY4-*met5Δ*', 'FY4-*met10Δ*', 'BY4741-*met5Δ*', 'BY4741-*met10Δ*'),
                              parsed = c('FY4', 'FY4-italic(met15Δ)', 'FY4-italic(met5Δ)','FY4-italic(met10Δ)',
                                         'BY4741-italic(met5Δ)', 'BY4741-italic(met10Δ)'))

##### SAVE COLONY SIZE DATA
save(data.510, strain.labs.510, file = "~/R/Projects/methionine/figures/final/met510.RData")

#####
# END
#####