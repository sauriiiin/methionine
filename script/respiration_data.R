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
expt.name <- "respiration"

##### GATHER AND CLEAN DATA FROM SQL
tables <- rbind(c('respiration_exp1R_MiM_Et_Rep2_JX__BinaThreZoom_cleaned_384_JPEG', 'MET-', 'URA+', 'EtOH', 'SD-Met-Cys+EtOH', '1'),
  c('respiration_exp1R_MiM_Glu_Rep2_JX__BinaThreZoom_cleaned_384_JPEG', 'MET-', 'URA+', 'Glucose', 'SD-Met-Cys+Glu', '1'),
  c('respiration_exp1R_PlM_Et_Rep2_JX__BinaThreZoom_cleaned_384_JPEG', 'MET+', 'URA+', 'EtOH', 'SD+Met-Cys+EtOH', '1'),
  c('respiration_exp1R_PlM_Glu_Rep2_JX__BinaThreZoom_cleaned_384_JPEG', 'MET+', 'URA+', 'Glucose', 'SD+Met-Cys+Glu', '1'),
  c('respiration_exp1R_Ura_Rep2_JX__binary_Zoom_cleaned_384_JPEG', 'MET+', 'URA-', 'Glucose', 'SD+Met-Ura+Glu', '1'),
  
  c('respiration_exp1R_MiM_Et_Rep1_BVO_binary_384_CLEAN', 'MET-', 'URA+', 'EtOH', 'SD-Met-Cys+EtOH', '2'),
  c('respiration_exp1R_MiM_Glu_Rep1_BVO_binary_384_CLEAN', 'MET-', 'URA+', 'Glucose', 'SD-Met-Cys+Glu', '2'),
  c('respiration_exp1R_PlM_Et_Rep1_BVO_binary_384_CLEAN', 'MET+', 'URA+', 'EtOH', 'SD+Met-Cys+EtOH', '2'),
  c('respiration_exp1R_PlM_Glu_Rep1b_BVO_binary_384_CLEAN', 'MET+', 'URA+', 'Glucose', 'SD+Met-Cys+Glu', '2'),
  c('respiration_exp1R_Ura_Rep1_BVO_384_binary_384_CLEAN', 'MET+', 'URA-', 'Glucose', 'SD+Met-Ura+Glu', '2'),
  
  c('respiration_exp2R_MiM_Et_Rep1_SBP_binary_384_CLEAN', 'MET-', 'URA+', 'EtOH', 'SD-Met-Cys+EtOH', '3'),
  c('respiration_exp2R_MiM_Glu_Rep1_SBP_binary_384_CLEAN', 'MET-', 'URA+', 'Glucose', 'SD-Met-Cys+Glu', '3'),
  c('respiration_exp2R_PlM_Et_Rep1_BVO_binary_384_CLEAN', 'MET+', 'URA+', 'EtOH', 'SD+Met-Cys+EtOH', '3'),
  c('respiration_exp2R_PlM_Glu_Rep1_BVO_binary_384_CLEAN', 'MET+', 'URA+', 'Glucose', 'SD+Met-Cys+Glu', '3'),
  c('respiration_exp2R_Ura_Rep1_SBP_binary_384_CLEAN', 'MET+', 'URA-', 'Glucose', 'SD+Met-Ura+Glu', '3'))

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
    temp$relative_cs[temp$hours == h] <-
      temp$average[temp$hours == h]/median(temp$average[temp$hours == h & temp$orf_name == 'FY4'], na.rm = T)
  }

  temp$methionine <- tables[i,2]
  temp$uracil <- tables[i,3]
  temp$carbon <- tables[i,4]
  temp$condition <- tables[i,5]
  temp$expt_rep <- tables[i,6]
  
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