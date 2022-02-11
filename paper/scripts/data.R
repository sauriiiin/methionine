
source('/home/sbp29/R/Projects/methionine/paper/scripts/initialize.R')

##### FIGURE 1C
info <- data.frame(rbind(c('rep1','MiM_Glu',366),
                         c('rep1','MiU_Glu',366),
                         c('rep1','MiM_Gal',366),
                         c('rep1','MiU_Gal',366),
                         c('rep1','MiM_Et',366),
                         c('rep1','MiU_Et',366),
                         c('rep1','PlM_Glu',366),
                         c('rep2','MiM_Glu',358),
                         c('rep2','MiU_Glu',358),
                         c('rep2','MiM_Gal',358),
                         c('rep2','MiU_Gal',358),
                         c('rep2','MiM_Et',358),
                         c('rep2','MiU_Et',358),
                         c('rep2','PlM_Glu',358)))
colnames(info) <- c('rep','condition','time')
data.cbn <- NULL
for (r in unique(info$rep)) {
  for (c in unique(info$condition[info$rep == r])) {
    for (t in unique(info$time[info$rep == r & info$condition == c])) {
      temp <- dbGetQuery(conn, sprintf('select a.pos, c.orf_name, a.hours, a.average, b.plate, b.col, b.row
                         from Branden.carbon_%s_%s_384_CLEAN a, Branden.carbon_%s_pos2coor b, Branden.carbon_%s_pos2orf_name c
                         where a.pos = b.pos and b.pos = c.pos
                         and c.orf_name != "BOR" and a.hours = %s
                         order by a.hours, b.plate, b.col, b.row', r, c, r, r, t))
      
      for (o in unique(temp$orf_name)) {
        temp$average[temp$orf_name == o][isoutlier(temp$average[temp$orf_name == o], 2)] <- NA
      }
      temp$relative_fitness <- temp$average/median(temp$average[temp$orf_name == 'FY4'], na.rm = T)
      
      temp$expt_rep <- r
      temp$condition <- c
      
      data.cbn <- rbind(data.cbn, temp)
    }
  }
}

info <- data.frame(rbind(c('R1','GLU'),
                         c('R2','GLU'),
                         c('R3','GLU'),
                         c('R1','GAL'),
                         c('R2','GAL'),
                         c('R3','GAL')))
colnames(info) <- c('rep','condition')
for (r in unique(info$rep)) {
  for (c in unique(info$condition[info$rep == r])) {
    temp <- dbGetQuery(conn, sprintf('select a.pos, c.orf_name, a.hours, a.average, b.plate, b.col, b.row
                                     from CARS_FS_MM_%s_%s_384_CLEAN a, MET_NS_pos2coor b, MET_NS_pos2orf_name c
                                     where a.pos = b.pos and b.pos = c.pos
                                     order by a.hours, b.plate, b.col, b.row', c, r))
    
    for (o in unique(temp$orf_name)) {
      temp$average[temp$orf_name == o][isoutlier(temp$average[temp$orf_name == o], 2)] <- NA
    }
    temp$relative_fitness <- temp$average/median(temp$average[temp$orf_name == 'FY4'], na.rm = T)
    
    temp$expt_rep <- r
    temp$condition <- c
    
    data.cbn <- rbind(data.cbn, temp)
  }
}
data.cbn$expt_rep[data.cbn$expt_rep == 'R1'] <- 'rep1'
data.cbn$expt_rep[data.cbn$expt_rep == 'R2'] <- 'rep2'
data.cbn$expt_rep[data.cbn$expt_rep == 'R3'] <- 'rep3'
data.cbn$orf_name[data.cbn$orf_name == 'FY4_met15del'] <- 'FY4-met15del'
data.cbn$orf_name[data.cbn$orf_name == 'FY4_met3del'] <- 'FY4-met3del'
data.cbn$bio_rep[data.cbn$row%%2==1 & data.cbn$col%%2==1] = '1'
data.cbn$bio_rep[data.cbn$row%%2==0 & data.cbn$col%%2==1] = '3'
data.cbn$bio_rep[data.cbn$row%%2==1 & data.cbn$col%%2==0] = '2'
data.cbn$bio_rep[data.cbn$row%%2==0 & data.cbn$col%%2==0] = '4'

data.cbn$condition <- factor(data.cbn$condition, levels = c('MiM_Glu','MiU_Glu',
                                                                      'MiM_Gal','MiU_Gal',
                                                                      'MiM_Et','MiU_Et',
                                                                      'GLU','GAL','PlM_Glu'))
strain.labs <- c('FY4','FY4-*met3Δ*','FY4-*met15Δ*','BY4742','BY4741')
data.cbn$orf_name <- factor(data.cbn$orf_name, levels = c("FY4","FY4-met3del","FY4-met15del","BY4742","BY4741"))
data.cbn <- data.cbn[data.cbn$expt_rep == 'rep1',c('condition','expt_rep','bio_rep','orf_name','relative_fitness')] 

head(data.cbn)
unique(data.cbn$condition)
data.cbn$carbon[str_detect(data.cbn$condition,'Glu')] <- 'Glucose'
data.cbn$carbon[str_detect(data.cbn$condition,'GLU')] <- 'Glucose'
data.cbn$carbon[str_detect(data.cbn$condition,'Gal')] <- 'Galactose'
data.cbn$carbon[str_detect(data.cbn$condition,'GAL')] <- 'Galactose'
data.cbn$carbon[str_detect(data.cbn$condition,'Et')] <- 'Ethanol'
data.cbn$carbon <- factor(data.cbn$carbon, levels = c('Glucose','Galactose','Ethanol'))

data.cbn$methionine[str_detect(data.cbn$condition,'MiM')] <- '-Met +Ura'
data.cbn$methionine[str_detect(data.cbn$condition,'PlM')] <- '+Met +Ura'
data.cbn$methionine[str_detect(data.cbn$condition,'MiU')] <- '+Met -Ura'
data.cbn$methionine[str_detect(data.cbn$condition,'GLU')] <- '-Met +Ura'
data.cbn$methionine[str_detect(data.cbn$condition,'GAL')] <- '-Met +Ura'
data.cbn$uracil[str_detect(data.cbn$condition,'MiM')] <- '+Ura'
data.cbn$uracil[str_detect(data.cbn$condition,'PlM')] <- '+Ura'
data.cbn$uracil[str_detect(data.cbn$condition,'MiU')] <- '-Ura'
data.cbn$uracil[str_detect(data.cbn$condition,'GLU')] <- '+Ura'
data.cbn$uracil[str_detect(data.cbn$condition,'GAL')] <- '+Ura'
data.cbn$cysteine[str_detect(data.cbn$condition,'MiM')] <- '-Cys'
data.cbn$cysteine[str_detect(data.cbn$condition,'PlM')] <- '-Cys'
data.cbn$cysteine[str_detect(data.cbn$condition,'MiU')] <- '-Cys'
data.cbn$cysteine[str_detect(data.cbn$condition,'GLU')] <- '-Cys'
data.cbn$cysteine[str_detect(data.cbn$condition,'GAL')] <- '-Cys'

data.cbn$base[data.cbn$condition %in% c('GLU','GAL')] <- 'SC'
data.cbn$base[is.na(data.cbn$base)] <- 'SD'

data.cbn[data.cbn$condition == 'MiU_Gal',]

strain.labs.cbn <- data.frame(orf_name = c("FY4","FY4-met3del","FY4-met15del","BY4742","BY4741"), 
                           labels = c('FY4','FY4-italic(met3Δ)','FY4-italic(met15Δ)','BY4742','BY4741'),
                           met_aux = c('Prototroph','Presumed Auxotroph','Presumed Auxotroph','Prototroph','Presumed Auxotroph'),
                           ura_aux = c('Prototroph','Prototroph','Prototroph','Presumed Auxotroph','Presumed Auxotroph'),
                           auxotrophy = c('None','Methionine','Methionine','Uracil','Both'),
                           stringsAsFactors = F)
strain.labs.cbn$auxotrophy <- factor(strain.labs.cbn$auxotrophy, levels = c('None','Methionine','Uracil','Both'))

data.cbn <- merge(data.cbn, strain.labs.cbn, by = 'orf_name')

## LEU DATA
tables <- data.frame(
  expt_id = c('SCmLeu', 'SCmLeu', 'SDmLeu', 'SDmLeu', 'YPDA', 'SDpMET', 'SDpMET'),
  arm = c('R1', 'R2', 'R1', 'R2', 'R1', 'R1', 'R2')
)

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
data.cbn.leu$base[str_detect(data.cbn.leu$condition, 'SC')] <- 'SC'
data.cbn.leu$base[str_detect(data.cbn.leu$condition, 'SD')] <- 'SD'
data.cbn.leu$base[is.na(data.cbn.leu$base)] <- 'YPDA'

# write.csv(rbind(data.cbn[,c(1,2,5,6,10)],
#                 data.cbn.leu[,c(8,11,10,13,14)] %>%
#                   filter(condition == 'SDmLeu')) %>%
#             filter(base == 'SD', orf_name %in% c('BY4741','BY4742')),
#           file = '/home/sbp29/R/Projects/methionine/paper/data/Figure1C.csv')


##### FIGURE 1E: SEE BELOW


##### FIGURE 1G
tables <- dbGetQuery(conn, 'SELECT table_name
                     FROM information_schema.tables
                     WHERE table_schema = "Branden"
                     and table_name like "Repeat_%_CLEAN" and table_rows < 10000')
temp <- str_split(tables$TABLE_NAME, '_', simplify = T) %>% data.frame()

tables <- cbind(tables, temp[,c(2:4)])
colnames(tables) <- c('table_name','aux','carbon','pin')

data.rp <- NULL
for (t in unique(tables$table_name)) {
  temp <- dbGetQuery(conn, sprintf("select a.*, b.density, b.plate, b.row, b.col, c.orf_name
                                   from Branden.%s a, Repeat_11_pos2coor b, Repeat_11_pos2orf_name c
                                   where a.pos = b.pos and b.pos = c.pos
                                   order by a.hours, a.pos", t))
  temp$aux <- tables$aux[tables$table_name == t]
  temp$carbon <- tables$carbon[tables$table_name == t]
  temp$pin <- tables$pin[tables$table_name == t]
  temp$table_name <- tables$table_name[tables$table_name == t]
  
  for (h in unique(temp$hours)) {
    for (o in unique(temp$orf_name)) {
      temp$average[temp$orf_name == o & temp$hours == h][isoutlier(temp$average[temp$orf_name == o & temp$hours == h], 2) |
                                                           isoutlier(temp$fitness[temp$orf_name == o & temp$hours == h], 2)] <- NA
    }
    temp$fitness[temp$hours == h] <- temp$average[temp$hours == h]/
      median(temp$average[temp$hours == h & temp$orf_name == 'FY4'], na.rm = T)
  }
  data.rp <- rbind(data.rp, temp)
}

temp <- data.rp %>%
  group_by(table_name) %>%
  summarise(.groups = 'keep', max_hrs = max(hours, na.rm = T)) %>%
  data.frame()
data.rp <- merge(data.rp, temp, by = 'table_name')

strain.labs.rp <- data.frame(orf_name = c("FY4","FY4-met3del","FY4-met15del","BY4742","BY4741"), 
                           labels = c('FY4','FY4-italic(met3Δ)','FY4-italic(met15Δ)','BY4742','BY4741'),
                           met_aux = c('Prototroph','Presumed Auxotroph','Presumed Auxotroph','Prototroph','Presumed Auxotroph'),
                           ura_aux = c('Prototroph','Prototroph','Prototroph','Presumed Auxotroph','Presumed Auxotroph'),
                           auxotrophy = c('None','Methionine','Methionine','Uracil','Both'),
                           stringsAsFactors = F)
strain.labs.rp$auxotrophy <- factor(strain.labs.rp$auxotrophy, levels = c('None','Methionine','Uracil','Both'))
strain.labs.rp$orf_name[strain.labs.rp$orf_name == 'FY4-met15del'] <- 'FY4-met15D'

data.rp$aux <- factor(data.rp$aux, levels = c('Ura', 'Met'))
data.rp$carbon <- factor(data.rp$carbon, levels = c('Glu', 'Gal'))
data.rp$cum_hrs <- data.rp$hours
data.rp$cum_hrs[data.rp$pin == 2] <- data.rp$cum_hrs[data.rp$pin == 2] +
  max(data.rp$cum_hrs[data.rp$pin == 1], na.rm = T)
data.rp$cum_hrs[data.rp$pin == 3] <- data.rp$cum_hrs[data.rp$pin == 3] +
  max(data.rp$cum_hrs[data.rp$pin == 2], na.rm = T)
data.rp$cum_hrs[data.rp$pin == 4] <- data.rp$cum_hrs[data.rp$pin == 4] +
  max(data.rp$cum_hrs[data.rp$pin == 3], na.rm = T)
data.rp$cum_hrs[data.rp$pin == 5] <- data.rp$cum_hrs[data.rp$pin == 5] +
  max(data.rp$cum_hrs[data.rp$pin == 4], na.rm = T)
data.rp$cum_hrs[data.rp$pin == 6] <- data.rp$cum_hrs[data.rp$pin == 6] +
  max(data.rp$cum_hrs[data.rp$pin == 5], na.rm = T)
data.rp$cum_hrs[data.rp$pin == 7] <- data.rp$cum_hrs[data.rp$pin == 7] +
  max(data.rp$cum_hrs[data.rp$pin == 6], na.rm = T)

data.rp.rf <- merge(data.rp %>%
                      filter(aux == 'Met', orf_name %in% c('FY4'), carbon == 'Glu', hours == max_hrs) %>%
                      group_by(pin, orf_name, hours, cum_hrs) %>%
                      summarize(average = median(average, na.rm = T), .groups = 'keep'),
                    
                    data.rp %>%
                      filter(aux == 'Met', orf_name %in% c('FY4-met15D'), carbon == 'Glu', hours == max_hrs) %>%
                      group_by(pin, orf_name, hours, cum_hrs) %>%
                      summarize(average = median(average, na.rm = T), .groups = 'keep'),
                    by = c('pin','hours','cum_hrs'), suffixes = c('_ref',''))
data.rp.rf$relative_fitness <- data.rp.rf$average/data.rp.rf$average_ref

# write.csv(merge(data.rp, strain.labs.rp,
#                 by = 'orf_name') %>%
#             filter(aux == 'Met', orf_name %in% c('FY4','FY4-met15D'), carbon == 'Glu'),
#           file = '/home/sbp29/R/Projects/methionine/paper/data/Figure1G.csv')

##### FIGURE 2B
tables <- rbind(c('jakes_mutants_pilot_MiM_Glu_384_JPEG_edge', 'jakes_mutants_pilot_pos2coor', 'jakes_mutants_pilot_pos2orf_name', 'Replicate 0', 'SD-Met'),
                c('jakes_mutants_c1_MiM_Glu_Binar_384_CLEAN_edge', 'jakes_mutants_pos2coor', 'jakes_mutants_pos2orf_name', 'Replicate 1', 'SD-Met'),
                c('jakes_mutants_c1_YPDA_Binar_384_CLEAN_edge', 'jakes_mutants_pos2coor', 'jakes_mutants_pos2orf_name', 'Replicate 1', 'YPDA'),
                c('jakes_mutants_c2_MiM_Glu_Binar_384_CLEAN_edge', 'jakes_mutants_pos2coor', 'jakes_mutants_pos2orf_name', 'Replicate 2', 'SD-Met'),
                c('jakes_mutants_c2_YPDA_Binar_384_CLEAN_edge', 'jakes_mutants_pos2coor', 'jakes_mutants_pos2orf_name', 'Replicate 2', 'YPDA'))

data.jm <- NULL
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
  
  data.jm <- rbind(data.jm,temp)
}
data.jm <- data.frame(data.jm)
head(data.jm)
data.jm <- data.jm[data.jm$orf_name != 'BOR',]

attempt.levels <- c('Replicate 0', 'Replicate 1', 'Replicate 2')
condition.levels <- c('YPDA', 'SD-Met')
data.jm$attempt <- factor(data.jm$attempt, levels = attempt.levels)
data.jm$condition <- factor(data.jm$condition, levels = condition.levels)

data.jm <- data.jm[!(data.jm$condition == 'YPDA' & data.jm$hours > 100),]
data.jm$time[data.jm$hours < 10] <- 't_0'
data.jm$time[data.jm$hours > 10] <- 't_final'

for (a in unique(data.jm$attempt)) {
  for (c in unique(data.jm$condition[data.jm$attempt == a])) {
    for (h in unique(data.jm$hours[data.jm$attempt == a & data.jm$condition == c])) {
      data.jm$relative_cs[data.jm$attempt == a & data.jm$hours == h & data.jm$condition == c] <- 
        data.jm$average[data.jm$attempt == a & data.jm$hours == h & data.jm$condition == c]/
        median(data.jm$average[data.jm$attempt == a & data.jm$hours == h & data.jm$condition == c & data.jm$orf_name == 'FY4'], na.rm = T)
      
    }
  }
}
strain.levels <- c('BY4742', 'BY4741', 'FY4', 'met15', 'met3', 'met2', 'met6', 'met13', 'cys4', 'met12', 'str3', 'yll')
strain.labs.jm <- c('BY4742', 'BY4741', 'FY4', 'FY4-*met15Δ*', 'FY4-*met3Δ*', 'FY4-*met2Δ*', 'FY4-*met6Δ*', 'FY4-*met13Δ*','FY4-*cys4Δ*',
                    'FY4-*met12Δ*', 'FY4-*str3Δ*', 'FY4-*yll058wΔ*')
strain.labs.jm.parse <- c('BY4742', 'BY4741', 'FY4', 'FY4-italic(met15Δ)', 'FY4-italic(met3Δ)', 'FY4-italic(met2Δ)',
                          'FY4-italic(met6Δ)', 'FY4-italic(met13Δ)','FY4-italic(cys4Δ)',
                          'FY4-italic(met12Δ)', 'FY4-italic(str3Δ)', 'FY4-italic(yll058wΔ)')
strain.labs.jm <- data.frame(orf_name = strain.levels, labels = strain.labs.jm, parsed = strain.labs.jm.parse)
strain.labs.jm$auxotrophy[strain.labs.jm$orf_name %in% c('BY4742','FY4')] <- 'Prototroph'
strain.labs.jm$auxotrophy[strain.labs.jm$orf_name %in% c('BY4741','met15','met3','met2','met6','met13','cys4')] <- 'Presumed Auxotroph'
strain.labs.jm$auxotrophy[strain.labs.jm$orf_name %in% c('met12')] <- 'Prototroph'
strain.labs.jm$auxotrophy[strain.labs.jm$orf_name %in% c('str3','yll')] <- 'Unknown'
strain.labs.jm$auxotrophy <- factor(strain.labs.jm$auxotrophy, levels = c('Prototroph', 'Presumed Auxotroph', 'Unknown'))

## MET5/10 DATA
tables <- rbind(c('MM', 'MET510_pos2coor', 'MET510_pos2orf_name', 'Replicate 1', 'SD-Met'),
                c('PM', 'MET510_pos2coor', 'MET510_pos2orf_name', 'Replicate 1', 'SD-Ura'), # PM and MU plates were mislabled in the INFO file
                c('MU', 'MET510_pos2coor', 'MET510_pos2orf_name', 'Replicate 1', 'SD+Met'),
                c('YPD', 'MET510_YPD_pos2coor', 'MET510_YPD_pos2orf_name', 'Replicate 1', 'YPDA'))

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

strain.labs.510 <- data.frame(orf_name = c('FY4','met15del','met5del','met10del','BY4741_met5del','BY4741_met10del'),
                              labels = c('FY4', 'FY4-*met15Δ*', 'FY4-*met5Δ*', 'FY4-*met10Δ*', 'BY4741-*met5Δ*', 'BY4741-*met10Δ*'),
                              parsed = c('FY4', 'FY4-italic(met15Δ)', 'FY4-italic(met5Δ)','FY4-italic(met10Δ)',
                                         'BY4741-italic(met5Δ)', 'BY4741-italic(met10Δ)'))


temp1 <- merge(data.510, strain.labs.510, by = 'orf_name')
temp1 <- temp1[((temp1$hours == 89 & temp1$condition == 'YPDA') |
                  (temp1$hours == 160 & temp1$condition == 'SD-Met')) &
                 temp1$orf_name %in% c('met5del','met10del'),]
temp2 <- merge(data.jm, strain.labs.jm, by = 'orf_name') %>%
  filter(time == 't_final', attempt != 'Replicate 0', condition %in% c('YPDA','SD-Met'))
temp2 <- temp2[,c(1:9,13,10,11,14,15)]
colnames(temp2) <- c(colnames(temp2)[1:9],'relative_fitness',colnames(temp2)[11:14])
data.jm.2 <- rbind(temp1, temp2)
data.jm.2$orf_name <- factor(data.jm.2$orf_name, levels = c('FY4','met15','met3','met5del','met10del','met2',
                                                            'met6','met13','cys4','str3','met12','yll',
                                                            'BY4742','BY4741'))
strain.labs.jm <- rbind(strain.labs.jm,
                        data.frame(orf_name = c('met5del','met10del'), labels = c('FY4-*met5Δ*','FY4-*met10Δ*'),
                                   parsed = c('FY4-italic(met5Δ)','FY4-italic(met10Δ)'), auxotrophy = c('Prototroph','Prototroph')))

# write.csv(data.jm.2[!(data.jm.2$orf_name %in% c('BY4742','BY4741')),],
#           file = '/home/sbp29/R/Projects/methionine/paper/data/Figure2B.csv')

##### FIGURE 2D
data.bis <- NULL
for (t in pzfx_tables("/home/sbp29/R/Projects/methionine/data/bismuth/bismuth_graphs_2.pzfx")) {
  temp <- read_pzfx("/home/sbp29/R/Projects/methionine/data/bismuth/bismuth_graphs_2.pzfx", table = t)
  temp <- melt(temp, variable.name = 'Strain', value.name = 'HS')
  temp$Condition = t
  data.bis <- rbind(data.bis, temp)
}
data.bis <- data.frame(data.bis)
data.bis$Strain <- as.character(data.bis$Strain)
data.bis$Strain <- factor(data.bis$Strain,
                          levels = c('FY4','FY4-met15D','FY4-met3D','FY4-met5D','FY4-met10D','FY4-met2D',
                                     'FY4-met6D','FY4-met13D','FY4-cys4D','FY4-str3D','FY4-met12D','FY4-yllD',
                                     'BY4742','BY4741'))
data.bis.2 <- data.frame(Strain = c('FY4-met5D','FY4-met10D','FY4-yllD',
                                    'FY4-met5D','FY4-met10D','FY4-yllD',
                                    'FY4-met5D','FY4-met10D','FY4-yllD',
                                    'FY4-met5D','FY4-met10D','FY4-yllD'),
                         HS = c(1,1,4,1,1,4,
                                1,1,4,1,1,4),
                         Condition = c('BiGGY','BiGGY','BiGGY','BiGGY','BiGGY','BiGGY',
                                       'SD-Met-Cys+Bi','SD-Met-Cys+Bi','SD-Met-Cys+Bi',
                                       'SD-Met-Cys+Bi','SD-Met-Cys+Bi','SD-Met-Cys+Bi'))
data.bis <- rbind(data.bis, data.bis.2)

# write.csv(data.bis,
#           file = '/home/sbp29/R/Projects/methionine/paper/data/Figure2D.csv')

##### FIGURE 2E
tables <- read_excel('/home/sbp29/RAW_Data/Methionine/NoSulfate_/NO_SUL_INFO.xlsx')
tables <- tables[order(tables$expt_id, tables$condition, tables$arm),] %>% filter(expt_id == 'SUL')
# head(tables)

data.ns <- NULL
for (s in unique(tables$stage_id)) {
  if (s %in% c('S3','S4','S5','Re2')) {
    p2c <- 'MET_NS_S3_pos2coor'
  } else if (s == 'Re1') {
    p2c <- 'MET_NS_Re1_pos2coor'
  } else {
    p2c <- 'MET_NS_pos2coor'
  }
  for (t in unique(tables$arm[tables$stage_id == s])) {
    temp <- dbGetQuery(conn, sprintf('select a.pos, a.hours, a.average, b.density, b.plate, b.row, b.col, c.orf_name
                                     from SUL_%s_%s_384_CLEAN a, %s b, MET_NS_pos2orf_name c
                                     where a.pos = b.pos and b.pos = c.pos
                                     order by a.hours, b.plate, b.col, b.row', s, t, p2c))
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
      if (sum(unique(temp$orf_name) == 'FY4') > 0) {
        temp$relative_fitness[temp$hours == h] <- temp$average[temp$hours == h]/
          median(temp$average[temp$hours == h & temp$orf_name == 'FY4'], na.rm = T)
      } else {
        temp$relative_fitness <- 0
      }
    }
    temp$arm <- t
    temp$stage <- s
    temp$condition <- unique(tables$condition[tables$arm == t])
    
    data.ns <- rbind(data.ns,temp)
  }
}
data.ns <- data.frame(data.ns)
# head(data.ns)

temp <- data.frame(cbind(str_split(unique(data.ns$condition), '_', simplify = T), unique(data.ns$condition)), stringsAsFactors = F)
colnames(temp) <- c('ynb_type','sulfate','media','base','condition')
temp[temp$condition == 'SD+Met+Glu_Agar',] <- c('Difco','+sulfate','SD+Met','Agar','SD+Met+Glu_Agar')
temp$sulfate[temp$sulfate == 'YNB'] <- '+sulfate'
temp$sulfate <- str_replace(temp$sulfate, 'sulfates', 'sulfate')
temp$methionine <- str_remove(temp$media, 'SD')
data.ns <- merge(data.ns, temp, by = 'condition')

data.ns$stage[data.ns$stage == 'FS'] <- 'S1'
data.ns$stage <- factor(data.ns$stage, levels = c('S1','S2','S3','S4','S5','Re1','Re2'))
data.ns$ynb_type <- factor(data.ns$ynb_type, levels = c('Difco','Home'))
data.ns$sulfate <- factor(data.ns$sulfate, levels = c('+sulfate','-sulfate'))
data.ns$base <- factor(data.ns$base, levels = c('Agar','Agarose'))
data.ns$methionine <- factor(data.ns$methionine, levels = c('+Met','-Met'))
data.ns$orf_name <- factor(data.ns$orf_name, levels = c('FY4','FY4_met15del','FY4_met3del','BY4742','BY4741'))

strain.labs.ns <- c('FY4','FY4-*met3Δ*','FY4-*met15Δ*','BY4742','BY4741')

data.ns$cum_hrs <- NULL
data.ns$cum_hrs[data.ns$stage == 'S1'] <- data.ns$hours[data.ns$stage == 'S1']
data.ns$cum_hrs[data.ns$stage == 'S2'] <- data.ns$hours[data.ns$stage == 'S2'] +
  max(data.ns$cum_hrs[data.ns$stage == 'S1'])
data.ns$cum_hrs[data.ns$stage == 'S3'] <- data.ns$hours[data.ns$stage == 'S3'] +
  max(data.ns$cum_hrs[data.ns$stage == 'S2'])
data.ns$cum_hrs[data.ns$stage == 'Re1'] <- data.ns$hours[data.ns$stage == 'Re1'] +
  max(data.ns$cum_hrs[data.ns$stage == 'S2'])
data.ns$cum_hrs[data.ns$stage == 'S4'] <- data.ns$hours[data.ns$stage == 'S4'] +
  max(data.ns$cum_hrs[data.ns$stage == 'S3'])
data.ns$cum_hrs[data.ns$stage == 'S5'] <- data.ns$hours[data.ns$stage == 'S5'] +
  max(data.ns$cum_hrs[data.ns$stage == 'S4'])
data.ns$cum_hrs[data.ns$stage == 'Re2'] <- data.ns$hours[data.ns$stage == 'Re2'] +
  max(data.ns$cum_hrs[data.ns$stage == 'S5'])

data.ns$id <- paste(data.ns$base, data.ns$ynb_type, data.ns$sulfate, data.ns$methionine, sep = '_')
data.ns$stage <- as.character(data.ns$stage)

data.ns$stage[data.ns$stage == 'Re1'] <- 'S3'
data.ns$stage[data.ns$stage == 'Re2'] <- 'S6'

data.ns$expt_rep[str_detect(data.ns$arm, 'R1')] <- 'R1'
data.ns$expt_rep[str_detect(data.ns$arm, 'R2')] <- 'R2'

data.ns2 <- data.ns[!(data.ns$id == 'Agar_Difco_+sulfate_+Met' & data.ns$stage == 'S1'),] %>%
  filter(id %in% c('Agar_Difco_+sulfate_+Met',
                   'Agarose_Home_+sulfate_-Met',
                   'Agarose_Home_-sulfate_-Met'),
         orf_name %in% c('FY4', 'FY4_met15del','FY4_met3del')) %>%
  data.frame()
data.ns2$id <- factor(data.ns2$id, levels = c('Agarose_Home_+sulfate_-Met',
                                              'Agarose_Home_-sulfate_-Met',
                                              'Agar_Difco_+sulfate_+Met'))
data.ns2$orf_name <- factor(data.ns2$orf_name, levels = c('FY4', 'FY4_met15del','FY4_met3del'))
# write.csv(data.ns2,
#           file = '/home/sbp29/R/Projects/methionine/paper/data/Figure2E.csv')


##### FIGURE 3B
## clusters
yll_align<-import.fasta("/home/acwach/Homology/yll_align.fas")
yll_align<-c(yll_align[-c(35,53,63)],yll_align[c(35,53,63)])
yll_dist<-mat.dis(yll_align, yll_align)

rnames<-paste(colnames(yll_dist),runif(length(colnames(yll_dist))),sep="_")
colnames(yll_dist)<-rnames
rownames(yll_dist)<-rnames

yll_mmds<-mmds(yll_dist, pc = 2, group.file = NULL)

km1<-kmeans.run(yll_mmds$coord, nb.clus = 3)
clus1<-substr(names(unlist(km1[[2]][1])),3,1000)
clus2<-substr(names(unlist(km1[[2]][2])),3,1000)
clus3<-substr(names(unlist(km1[[2]][3])),3,1000)

yll058w_colors<-rep("black",length(rnames))
yll058w_colors[which(rnames %in% clus1)]<-"Cluster 1"
yll058w_colors[which(rnames %in% clus2)]<-"Cluster 2"
yll058w_colors[which(rnames %in% clus3)]<-"Cluster 3"

yll_taxa<-array()
yll_taxa<-rep("Ascomycota",length(rnames))
yll_taxa[c(3,10,13,103)]<-"Bacteria"
yll_taxa[c(109,113,114,287,288,289,297,298,299)]<-"Fungal outgroups"

df_yll<-data.frame(PC1=yll_mmds$coord$PC1,PC2=yll_mmds$coord$PC2,km_clust=yll058w_colors,spp_name=rnames,taxa=yll_taxa)
df_yll$km_clust<-factor(df_yll$km_clust,levels=c("Cluster 1","Cluster 2","Cluster 3"))#,"Cluster 4"))

unique(df_yll$km_clust)

## tree
tree0<-read.tree(text='(Lipomycetaceae,(Trigonopsidaceae,(Dipodascaceae,
((CUG-Ser1 clade,Pichiaceae),((Saccharomycopsis,Ascoidea),((Cyberlindnera,Barnettozyma),(Hanseniaspora,((Lachancea,(Eremothecium,Kluyveromyces)),((Torulaspora,Zygotorulaspora),(Tetrapisispora,((Naumovozyma,Kazachstania),(Nakaseomyces,Saccharomyces))))))))))));
')
tree1<-read.tree(text='(Lipomycetaceae,(Trigonopsidaceae,(Dipodascaceae,
(((Cephaloascus,((Meyerozyma,Yamadazyma),((Metschnikowia,Hyphopichia),(Debaryomyces,(Priceomyces,((Spathaspora,Scheffersomyces),(Teunomyces,Suhomyces) ) )))))
,Pichiaceae),((Saccharomycopsis,Ascoidea),((Cyberlindnera,Barnettozyma),(Hanseniaspora,((Lachancea,(Eremothecium,Kluyveromyces)),((Torulaspora,Zygotorulaspora),(Tetrapisispora,((Naumovozyma,Kazachstania),(Nakaseomyces,Saccharomyces))))))))))));
')

cluster2_present<-c("Lachancea","Zygotorulaspora","Torulaspora","Saccharomyces")
cluster1_present<-c("Lachancea","Tetrapisispora","Nakaseomyces","Kazachstania","Hanseniaspora","Torulaspora","Zygotorulaspora","Kluyveromyces","Eremothecium","Naumovozyma","Saccharomyces")
cluster3_present<-c("Trigonopsidaceae","Lipomycetaceae","Dipodascaceae","Pichiaceae","Saccharomycopsis","Debaryomyces","Yamadazyma",
                    "Sporopachydermia","Cyberlindnera","Alloascoidea","Starmera","Nakazawaea","Cephaloascus","Priceomyces","Barnettozyma","Metschnikowia","Peterozyma",
                    "Cephaloascus","Meyerozyma","Spathaspora","Scheffersomyces","Hyphopichia")

clus1colors<-rep("transparent",30)
clus1colors[c(1,2,3,4,5,6,7,8,9,10,11)]<-"#006666"
clus1colors2<-rep("transparent",30)
clus1colors2[c(1,2,3,4,5,6,7,8,9,10,11)]<-"#9E9E9E"

clus2colors<-rep("transparent",30)
clus2colors[c(1,6,7,10)]<-"#CC33FF"
clus2colors2<-rep("transparent",30)
clus2colors2[c(1,6,7,10)]<-"#212121"

clus3colors<-rep("transparent",30)
clus3colors[c(12,13,15,18,19,20,21,22,23,24,25,26,27,28,29,30)]<-"#99CCFF"
clus3colors2<-rep("transparent",30)
clus3colors2[c(12,13,15,18,19,20,21,22,23,24,25,26,27,28,29,30)]<-"#9E9E9E"


##### FIGURE 3D
tables <- data.frame(
  expt_id = c('PV_FY_MM', 'PV_FY_MM', 'PV_FY_PM', 'PV_FY_PM'),
  stage_id = c('PS1_2', 'PS1_2', 'PS1_2', 'PS1_2'),
  arm = c('R1', 'R2', 'R1', 'R2')
)

data.pv2 <- NULL
for (e in unique(tables$expt_id)) {
  for (s in unique(tables$stage_id[tables$expt_id == e])) {
    for (a in unique(tables$arm[tables$expt_id == e & tables$stage_id == s])) {
      temp <- dbGetQuery(conn, sprintf('select a.pos, a.hours, a.average, b.density, b.plate, b.row, b.col, c.orf_name
                                       from %s_%s_%s_384_CLEAN a, PV2_FY_pos2coor b, PV2_FY_pos2orf_name c
                                       where a.pos = b.pos and b.pos = c.pos
                                       order by a.hours, b.plate, b.col, b.row', e, s, a))
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
                                                                                                      temp$orf_name %in% c('Plasmid_1')], na.rm = T)
      }
      temp$arm <- e
      temp$stage <- s
      temp$expt_rep <- a
      # temp$condition <- unique(tables$condition[tables$expt_id == e & tables$stage == s & tables$arm == a])
      
      data.pv2 <- rbind(data.pv2,temp)
    }
  }
}
data.pv2 <- data.frame(data.pv2)

# write.csv(data.pv2,
#           file = '/home/sbp29/R/Projects/methionine/paper/data/Figure3D.csv')



##### FIGURE 3E
data.bioc <- read_pzfx("/home/sbp29/R/Projects/methionine/data/gina/072321_gina_data.pzfx", table = "Data 2")
data.bioc <- melt(data.bioc, id.vars = 'Time (min)', variable.name = 'ID', value.name = 'uM')
data.bioc <- cbind(data.bioc, str_split(data.bioc$ID, '_', simplify = T))
colnames(data.bioc) <- c(colnames(data.bioc)[1:3], 'Sample', 'Replicate')
data.bioc$Sample <- factor(data.bioc$Sample, levels = c('Met15', 'Yll058w', 'None'))

# write.csv(data.bioc,
#           file = '/home/sbp29/R/Projects/methionine/paper/data/Figure3E.csv')



##### FIGURE 4B
strain.labs.h2s <- read.csv(file = 'data/h2s/FeEDTA_samples.csv')
data.h2s <- read.csv(file = 'data/h2s/FeEDTA_readings.csv')

data.h2s.re <- NULL
for (t in pzfx_tables("/home/sbp29/R/Projects/methionine/data/h2s/chelator_liquid.pzfx")) {
  temp <- read_pzfx("/home/sbp29/R/Projects/methionine/data/h2s/chelator_liquid.pzfx", table = t)
  temp <- melt(temp, id.vars = 'Time (Minutes)', variable.name = 'Strain', value.name = 'OD')
  temp$Set = t
  data.h2s.re <- rbind(data.h2s.re, temp)
}
data.h2s.re <- data.frame(data.h2s.re)
data.h2s.re <- data.h2s.re[data.h2s.re$Set == 'All_Samples',]

temp <- str_split(data.h2s.re$Strain, '_', simplify = T)
temp[,4] <- temp[!str_detect(temp[,3], 'R'),3]
temp[,3] <- temp[str_detect(temp[,2], 'R'),2]
temp[str_detect(temp[,2], 'R'),2] <- ''

temp <- data.frame(temp, stringsAsFactors = F)
colnames(temp) <- c('ORF','FeEDTA','Replicate','Reading')
temp$FeEDTA <- as.character(temp$FeEDTA)
temp$FeEDTA[temp$FeEDTA == 'Fe'] <- 'Present'
temp$FeEDTA[temp$FeEDTA == ''] <- 'Absent'

data.h2s.re <- cbind(data.h2s.re, temp)

strain.labs.h2sr <- read.csv(file = 'data/h2s/FeEDTA_repeat_samples.csv')
data.h2sr <- read.csv(file = 'data/h2s/FeEDTA_repeat_readings.csv')

## TOP
temp1 <- merge(data.h2s %>%
                 melt(id.vars = 'Time', variable.name = 'Flask', value.name = 'OD'), strain.labs.h2s, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Time == 3075, Expected.Initial.OD == 0.1)
temp1$escapee[temp1$Replicate == 'Rep_1' & temp1$FeEDTA == 'Absent' & temp1$Strain == 'met15del'] <- 'Yes'
temp1$escapee[is.na(temp1$escapee)] <- 'No'
# filter(Time == 3075, Expected.Initial.OD == 0.1, Replicate != 'Rep_1') 

temp2 <- merge(data.h2sr %>%
                 melt(id.vars = c('Time', 'Attempt'), variable.name = 'Flask', value.name = 'OD') %>%
                 filter(Attempt == 'Original'), strain.labs.h2sr, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Time == 3075, Expected.Initial.OD == 0.1, Strain %in% c('FY4','met15del'))
temp2 <- temp2[,-3]
temp2$escapee <- 'No'

rbind(temp1, temp2) %>%
  filter(escapee == 'No') %>%
  group_by(Strain, FeEDTA) %>%
  summarize(OD_m = mean(OD, na.rm = T), OD_sd = sd(OD, na.rm = T), .groups = 'keep')

# write.csv(rbind(temp1, temp2) %>% filter(escapee == 'No'),
#           file = '/home/sbp29/R/Projects/methionine/paper/data/Figure4B_OriginalCultures.csv')

## BOTTOM
temp3 <- data.h2s.re[data.h2s.re$Time..Minutes. == 3075,] %>%
  group_by(FeEDTA, ORF, Replicate) %>%
  summarize(OD = median(OD, na.rm = T), .groups = 'keep') %>%
  data.frame()
temp3$escapee[temp3$FeEDTA == 'Absent' & temp3$Replicate == 'R2' & temp3$ORF == 'met15D'] <- 'Yes'
temp3$escapee[temp3$FeEDTA == 'Present' & temp3$Replicate == 'R3' & temp3$ORF == 'met15D'] <- 'Yes'
temp3$escapee[is.na(temp3$escapee)] <- 'No'

temp4 <- merge(data.h2sr %>%
                 melt(id.vars = c('Time', 'Attempt'), variable.name = 'Flask', value.name = 'OD') %>%
                 filter(Attempt == 'Repassage'), strain.labs.h2sr, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Time == 3075, Expected.Initial.OD == 0.1, Strain %in% c('FY4','met15del')) %>%
  group_by(FeEDTA, Strain, Replicate) %>%
  summarise(OD = OD, .groups = 'keep') %>%
  data.frame()
colnames(temp4) <- colnames(temp3)[-5]
temp4$ORF <- as.character(temp4$ORF)
temp4$ORF[temp4$ORF == 'met15del'] <- 'met15D'
temp4$Replicate <- as.character(temp4$Replicate)
temp4$Replicate <- str_replace(temp4$Replicate, 'Rep_','R')
temp4$escapee[temp4$FeEDTA == 'Present' & temp4$Replicate %in% c('R1','R3') & temp4$ORF == 'met15D'] <- 'Yes'
temp4$escapee[is.na(temp4$escapee)] <- 'No'

# write.csv(rbind(temp3, temp4) %>% filter(escapee == 'No'),
#           file = '/home/sbp29/R/Projects/methionine/paper/data/Figure4B_RepassagedCultures.csv')

##### FIGURE 1E
temp5 <- merge(data.h2sr %>%
                 melt(id.vars = c('Time', 'Attempt'), variable.name = 'Flask', value.name = 'OD'), 
               strain.labs.h2sr, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Time == 3075, Expected.Initial.OD == 0.1, Strain %in% c('BY4742','BY4741'),
         FeEDTA == 'Absent') %>%
  group_by(Attempt, Strain, Replicate) %>%
  summarise(OD = OD, .groups = 'keep') %>%
  data.frame()
temp5$Strain <- factor(temp5$Strain, levels = c('BY4742','BY4741'))

# write.csv(temp5 %>% filter(Attempt == 'Original'),
#           file = '/home/sbp29/R/Projects/methionine/paper/data/Figure1E.csv')


##### FIGURE 4D
strain.labs.res <- data.frame(orf_name = c('FY4', 'FY4_pet', 'FY4-met15del', 'FY4-met15del_pet',
                                           'BY4742', 'BY4742_pet', 'BY4741', 'BY4741_pet'), 
                              labels = c('FY4', 'FY4~(rho)' ,'FY4-italic(met15Δ)', 'FY4-italic(met15Δ)~(rho)',
                                         'BY4742', 'BY4742~(rho)', 'BY4741', 'BY4741~(rho)'),
                              met_aux = c('Prototroph','Prototroph','Presumed Auxotroph','Presumed Auxotroph',
                                          'Prototroph','Prototroph','Presumed Auxotroph','Presumed Auxotroph'),
                              pet = c('No','Yes','No','Yes','No','Yes','No','Yes'))
strain.levels = c('FY4', 'FY4_pet', 'FY4-met15del', 'FY4-met15del_pet', 'BY4742', 'BY4742_pet', 'BY4741', 'BY4741_pet')
condition.levels = c('SD+Met-Cys+Glu','SD-Met-Cys+Glu','SD+Met-Cys+EtOH','SD-Met-Cys+EtOH','SD+Met-Ura+Glu')


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

data.res <- NULL
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
  
  data.res <- rbind(data.res,temp)
}
data.res <- data.frame(data.res)
head(data.res)

data.res$bio_rep[data.res$row%%2==1 & data.res$col%%2==1] = '1'
data.res$bio_rep[data.res$row%%2==0 & data.res$col%%2==1] = '3'
data.res$bio_rep[data.res$row%%2==1 & data.res$col%%2==0] = '2'
data.res$bio_rep[data.res$row%%2==0 & data.res$col%%2==0] = '4'
data.res <- data.res[data.res$orf_name != 'BOR',]
data.res$orf_name <- factor(data.res$orf_name, levels = strain.levels)
data.res$condition <- factor(data.res$condition, levels = condition.levels)

data.res <- data.res[!(data.res$hours %in% c(46.23, 87.35, 25.88, 16.75, 54.36, 37.13)),]

data.gc <- data.res %>%
  group_by(expt_rep, condition, orf_name, bio_rep, pos, hours) %>%
  summarise(average = median(average, na.rm = T), .groups = 'keep') %>%
  data.frame()

data.res.gc <- NULL
col.names <- NULL
for (e in unique(data.gc$expt_rep)) {
  for (o in unique(data.gc$orf_name[data.gc$expt_rep == e])) {
    for (c in unique(data.gc$condition[data.gc$expt_rep == e & data.gc$orf_name == o])) {
      for (b in unique(data.gc$bio_rep[data.gc$expt_rep == e & data.gc$orf_name == o & data.gc$condition == c])) {
        for (p in unique(data.gc$pos[data.gc$expt_rep == e & data.gc$orf_name == o & data.gc$condition == c & data.gc$bio_rep == b])) {
          temp <- data.gc[data.gc$expt_rep == e &
                            data.gc$orf_name == o & data.gc$condition == c &
                            data.gc$bio_rep == b & data.gc$pos == p,]
          if (sum(is.na(temp$average)) <= 5) {
            lo <- loess.smooth(temp$hours, log(temp$average),
                               span = 0.6, evaluation = 100, degree = 2,
                               family = 'gaussian')
            data.res.gc <- cbind(data.res.gc,exp(lo$y))
            col.names <- cbind(col.names,paste(e,o,c,b,p,sep = ','))
          }
        }
      }
    }
  }
}
data.res.gc <- cbind(lo$x, data.res.gc)
data.res.gc <- data.frame(data.res.gc)
colnames(data.res.gc) <- c('Time',col.names)

data.res.gc <- melt(data.res.gc, id.vars = 'Time', variable.name = 'sample', value.name = 'cs')
temp <- str_split(data.res.gc$sample, ',', simplify = T)
colnames(temp) <- c('expt_rep','orf_name','condition','bio_rep','pos')
data.res.gc <- cbind(temp, data.res.gc)

for (c in unique(data.res.gc$condition)) {
  for (e in unique((data.res.gc$expt_rep[data.res.gc$condition == c]))) {
    max.t <- max(data.res.gc$Time[data.res.gc$condition == c & data.res.gc$expt_rep == e])
    med.cs <- mean(data.res.gc$cs[data.res.gc$condition == c & data.res.gc$expt_rep == e & 
                                    data.res.gc$Time == max.t & data.res.gc$orf_name == 'FY4'], na.rm = T)
    data.res.gc$rel_cs[data.res.gc$condition == c & data.res.gc$expt_rep == e] <- data.res.gc$cs[data.res.gc$condition == c & data.res.gc$expt_rep == e]/med.cs
  }
}
data.res.gc$condition <- factor(data.res.gc$condition, levels = condition.levels)
data.res.gc$orf_name <- factor(data.res.gc$orf_name, levels = strain.levels)

conds <- NULL
conds$condition <- unique(data.res.gc$condition)
conds$base[str_detect(conds$condition, 'SD')] <- 'SD'
conds$methionine[str_detect(conds$condition, '-Met')] <- '-Met'
conds$methionine[is.na(conds$methionine)] <- '+Met'
conds$carbon[str_detect(conds$condition, 'Glu')] <- 'Glucose'
conds$carbon[is.na(conds$carbon)] <- 'Ethanol'
conds$cysteine <- '-Cys'
conds <- data.frame(conds)

data.res.gc <- merge(data.res.gc, conds, by = 'condition')
data.res.gc$methionine <- factor(data.res.gc$methionine, levels = c('+Met','-Met'))
data.res.gc$carbon <- factor(data.res.gc$carbon, levels = c('Glucose','Ethanol'))

# write.csv(merge(data.res.gc[str_detect(data.res.gc$orf_name, 'FY'),],
#                 strain.labs.res, by = 'orf_name') %>%
#             filter(expt_rep %in% c("1","2"), condition != 'SD+Met-Ura+Glu'),
#           file = '/home/sbp29/R/Projects/methionine/paper/data/Figure4D.csv')


##### SUPPLEMENTARY FIGURES #####
##### FIGURE S3
data.cbn.leu$orf_name[data.cbn.leu$orf_name == 'met15'] <- 'FY4-met15del'
# write.csv(rbind(data.cbn[,c(1,2,5,6,10)],
#                 data.cbn.leu[,c(8,11,10,13,14)] %>%
#                   filter(condition == 'SCmLeu')) %>%
#             filter(base == 'SC', orf_name != 'FY4-met3del'),
#           file = '/home/sbp29/R/Projects/methionine/paper/data/FigureS3.csv')


##### FIGURE S4
# write.csv(rbind(data.cbn[,c(1,2,5,6,10)],
#                 data.cbn.leu[,c(8,11,10,13,14)] %>%
#                   filter(condition == 'SDmLeu')) %>%
#             filter(base == 'SD', orf_name %in% c('FY4','FY4-met15del')),
#           file = '/home/sbp29/R/Projects/methionine/paper/data/FigureS4.csv')


##### FIGURE S5
# write.csv(rbind(data.cbn[,c(1,2,5,6,10)],
#                 data.cbn.leu[,c(8,11,10,13,14)] %>%
#                   filter(condition == 'SDmLeu')) %>%
#             filter(base == 'SD', orf_name == 'FY4-met3del'),
#           file = '/home/sbp29/R/Projects/methionine/paper/data/FigureS5.csv')

##### FIGURE S9
data.mdl <- read.csv('/home/sbp29/R/Projects/methionine/data/modeling/YLL_simulation.csv', stringsAsFactors = F)
data.mdl$label[data.mdl$Model == 'A'] <- 'A. Default'
data.mdl$label[data.mdl$Model == 'B'] <- 'B. A - All YLL058W reactions'
data.mdl$label[data.mdl$Model == 'C'] <- 'C. B + Hypothesized YLL058W reaction'
data.mdl$label[data.mdl$Model == 'D'] <- 'D. C - All MET15 reactions'

# write.csv(data.mdl, file = '/home/sbp29/R/Projects/methionine/paper/data/FigureS9.csv')


