##### RESPIRATION EXPERIMENT - DATA
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 05/03/2021 

##### INITIALIZE
library(RMariaDB)

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- "~/R/Projects/methionine/data"
expt.name <- "resp_1_2"

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
save(data, file = sprintf("%s/%s/colonysizes.RData", out_path, expt.name))

#####
# END
#####