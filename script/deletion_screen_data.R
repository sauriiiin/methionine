##### DELETION SCREEN DATA
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 07/15/2021 

##### INITIALIZE
library(RMariaDB)

out_path <- "~/R/Projects/methionine/data"
expt.name <- "deletion"

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

##### GATHER DATA
stages <- c('Pre-Screen #1','Pre-Screen #2','Final Screen')
arms <- c('SD-Met-Cys+Gal','SD+Met-Cys+Gal')
phenotypes <- c('Beneficial','Neutral','Deleterious')

info <- data.frame(rbind(c('PS1','MM',1536,'MET_DEL'), c('PS1','PM',1536,'MET_DEL'),
                         c('PS2','MM',1536,'MET_DEL'), c('PS2','PM',1536,'MET_DEL'),
                         c('FS','MM',6144,'MET_DEL_FS_MM'), c('FS','PM',6144,'MET_DEL')))
colnames(info) <- c('stage','arm','density','p2c')

pgs <- dbGetQuery(conn, 'select orf_name from PROTOGENES where pg_2012 = 1')

data <- NULL
data.sum <- NULL
for (i in seq(1,dim(info)[1])) {
  temp <- dbGetQuery(conn, sprintf('select a.*, b.density, b.plate, b.row, b.col
                                   from MET_DEL_%s_%s_%s_FITNESS a,
                                   %s_pos2coor b
                                   where a.pos = b.pos
                                   order by a.hours, b.plate, b.col, b.row',
                                   info[i,1],info[i,2],info[i,3],
                                   info[i,4]))
  temp <- temp[!is.na(temp$orf_name) & temp$orf_name != 'NULL',]
  
  temp2 <- dbGetQuery(conn, sprintf('select a.*, b.p
                                    from MET_DEL_%s_%s_%s_FITNESS_STATS a
                                    left join
                                    MET_DEL_%s_%s_%s_PVALUE b
                                    on a.hours = b.hours and a.strain_id = b.strain_id
                                    order by a.hours',
                                    info[i,1],info[i,2],info[i,3],
                                    info[i,1],info[i,2],info[i,3]))
  
  for (h in unique(temp$hours)) {
    for (o in unique(temp$orf_name)) {
      temp$average[temp$orf_name == o & temp$hours == h][isoutlier(temp$average[temp$orf_name == o & temp$hours == h], 2) |
                                                           isoutlier(temp$fitness[temp$orf_name == o & temp$hours == h], 2)] <- NA
      temp$fitness[temp$orf_name == o & temp$hours == h][isoutlier(temp$average[temp$orf_name == o & temp$hours == h], 2) |
                                                           isoutlier(temp$fitness[temp$orf_name == o & temp$hours == h], 2)] <- NA
    }
  }
  temp$stage <- info[i,1]
  temp$arm <- info[i,2]
  
  temp2$stage <- info[i,1]
  temp2$arm <- info[i,2]
  
  data <- rbind(data,temp)
  data.sum <- rbind(data.sum,temp2)
}
data <- data.frame(data)
data$orf_type[data$orf_name %in% pgs$orf_name] <- 'Proto-gene'
data$orf_type[is.na(data$orf_type)] <- 'Gene'
data$arm <- as.character(data$arm)
data$stage <- as.character(data$stage)
data$arm[data$arm == 'MM'] <- 'SD-Met-Cys+Gal'
data$arm[data$arm == 'PM'] <- 'SD+Met-Cys+Gal'
data$stage[data$stage == 'PS1'] <- 'Pre-Screen #1'
data$stage[data$stage == 'PS2'] <- 'Pre-Screen #2'
data$stage[data$stage == 'FS'] <- 'Final Screen'

data$stage <- factor(data$stage, levels = stages)
data$arm <- factor(data$arm, levels = arms)


data.sum <- data.frame(data.sum)
data.sum$orf_type[data.sum$orf_name %in% pgs$orf_name] <- 'Proto-gene'
data.sum$orf_type[is.na(data.sum$orf_type)] <- 'Gene'
data.sum$arm <- as.character(data.sum$arm)
data.sum$stage <- as.character(data.sum$stage)
data.sum$arm[data.sum$arm == 'MM'] <- 'SD-Met-Cys+Gal'
data.sum$arm[data.sum$arm == 'PM'] <- 'SD+Met-Cys+Gal'
data.sum$stage[data.sum$stage == 'PS1'] <- 'Pre-Screen #1'
data.sum$stage[data.sum$stage == 'PS2'] <- 'Pre-Screen #2'
data.sum$stage[data.sum$stage == 'FS'] <- 'Final Screen'
data.sum$phenotype[data.sum$p <= 0.05 & data.sum$cs_mean > 1] <- 'Beneficial'
data.sum$phenotype[data.sum$p <= 0.05 & data.sum$cs_mean < 1] <- 'Deleterious'
data.sum$phenotype[is.na(data.sum$phenotype)] <- 'Neutral'

data.sum$stage <- factor(data.sum$stage, levels = stages)
data.sum$arm <- factor(data.sum$arm, levels = arms)
data.sum$phenotype <- factor(data.sum$phenotype, levels = phenotypes)

for (a in arms) {
  for (s in stages) {
    data.sum$saturation[data.sum$arm == a & data.sum$stage == s] <- 
      max(data.sum$hours[data.sum$arm == a & data.sum$stage == s])
  }
}

save(data, data.sum, file = sprintf('%s/%s/data.RData',out_path,expt.name))
