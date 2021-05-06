##### JAKES MUTANATS EXPERIMENT - DATA
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 05/05/2021 

##### INITIALIZE
library(RMariaDB)

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- "~/R/Projects/methionine/data"
expt.name <- "jakes"

##### GATHER AND CLEAN DATA FROM SQL
tables <- rbind(c('jakes_mutants_pilot_MiM_Glu_384_JPEG_edge', 'jakes_mutants_pilot_pos2coor', 'jakes_mutants_pilot_pos2orf_name', 'pilot', 'SD-MET+Glucose'),
                c('jakes_mutants_c1_MiM_Glu_Binar_384_CLEAN_edge', 'jakes_mutants_pos2coor', 'jakes_mutants_pos2orf_name', 'attempt1', 'SD-MET+Glucose'),
                c('jakes_mutants_c1_YPDA_Binar_384_CLEAN_edge', 'jakes_mutants_pos2coor', 'jakes_mutants_pos2orf_name', 'attempt1', 'YPDA'),
                c('jakes_mutants_c2_MiM_Glu_Binar_384_CLEAN_edge', 'jakes_mutants_pos2coor', 'jakes_mutants_pos2orf_name', 'attempt2', 'SD-MET+Glucose'),
                c('jakes_mutants_c2_YPDA_Binar_384_CLEAN_edge', 'jakes_mutants_pos2coor', 'jakes_mutants_pos2orf_name', 'attempt2', 'YPDA'))

data <- NULL
for (i in seq(1,dim(tables)[1])) {
  temp <- dbGetQuery(conn, sprintf('select a.pos, a.hours, a.average, b.density, b.plate, b.row, b.col, c.orf_name
                                   from Branden.%s a, Branden.%s b, Branden.%s c
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
  }
  
  temp$attempt <- tables[i,4]
  temp$condition <- tables[i,5]
  
  data <- rbind(data,temp)
}
data <- data.frame(data)
head(data)

##### SAVE COLONY SIZE DATA
save(data, file = sprintf("%s/%s/colonysizes.RData", out_path, expt.name))

#####
# END
#####