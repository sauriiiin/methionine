##### PLASMID VALIDATION EXPERIMENT - DATA
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 08/10/2021 

##### INITIALIZE
library(RMariaDB)
library(readxl)
library(stringr)
library(dplyr)

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- "~/R/Projects/methionine/data"
expt.name <- "plasmidval"

tables <- read_excel('/home/sbp29/RAW_Data/Methionine/PlasmidValidation2/PV2_FY_INFO.xlsx')
tables <- tables[order(tables$expt_id, tables$condition, tables$arm),] %>% data.frame()
head(tables)

strain.labs.pv2 <- read.csv('/home/sbp29/R/Projects/methionine/data/plasmidval/PV2_Strains.csv', stringsAsFactors = F)
head(strain.labs.pv2)

##### GATHER AND CLEAN DATA FROM SQL
data <- NULL
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
      
      if (e == 'PV_FY_MM') {
        if (s == 'PS2') {
          if (a == 'R2') {
            temp$average[temp$plate == 10] <- NA
          }
        }
      }

      for (h in unique(temp$hours)) {
        for (o in unique(temp$orf_name)) {
          for (b in unique(temp$bio_rep)) {
            temp$average[temp$orf_name == o & temp$hours == h & temp$bio_rep == b][isoutlier(temp$average[temp$orf_name == o & temp$hours == h & temp$bio_rep == b], 2)] <- NA
          }
        }
        temp$relative_fitness[temp$hours == h] <- temp$average[temp$hours == h]/median(temp$average[temp$hours == h & 
                                                                                                      temp$orf_name %in% c('Plasmid_1','Plasmid_2')], na.rm = T)
      }
      temp$arm <- e
      temp$stage <- s
      temp$expt_rep <- a
      # temp$condition <- unique(tables$condition[tables$expt_id == e & tables$stage == s & tables$arm == a])

      data <- rbind(data,temp)
    }
  }
}
data <- data.frame(data)
head(data)

# data.sum <- data %>%
#   group_by(stage,condition,expt_rep,hours,orf_name,deletion1,deletion2,plasmid_backbone,plasmid_orf,bio_rep) %>%
#   summarize(relative_fitness = median(relative_fitness, na.rm = T),
#             cs = median(average, na.rm = T), .groups = 'keep') %>%
#   data.frame()

##### SAVE COLONY SIZE DATA
# save(data, data.sum, file = sprintf("%s/%s/colonysizes.RData", out_path, expt.name))
# save(anova.res, file = sprintf("%s/%s/stats.RData", out_path, expt.name))

#####
# END
#####


##### PS1_2
tables <- data.frame(
  expt_id = c('PV_FY_MM', 'PV_FY_MM', 'PV_FY_PM', 'PV_FY_PM'),
  stage_id = c('PS1_2', 'PS1_2', 'PS1_2', 'PS1_2'),
  arm = c('R1', 'R2', 'R1', 'R2')
)

##### GATHER AND CLEAN DATA FROM SQL
data2 <- NULL
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
                                                                                                      temp$orf_name %in% c('Plasmid_1','Plasmid_2')], na.rm = T)
      }
      temp$arm <- e
      temp$stage <- s
      temp$expt_rep <- a
      # temp$condition <- unique(tables$condition[tables$expt_id == e & tables$stage == s & tables$arm == a])
      
      data2 <- rbind(data2,temp)
    }
  }
}
data2 <- data.frame(data2)
data <- rbind(data, data2)


##### PLOTTING THE COMBINED RESULTS
strain.labs.pv2$Strain <- factor(strain.labs.pv2$Strain, levels = strain.labs.pv2$Strain)
data$orf_name <- factor(data$orf_name, levels = c('Plasmid_1','Plasmid_2','Plasmid_3','Plasmid_4',
                                                  'Plasmid_9','Plasmid_10','Plasmid_5','Plasmid_6',
                                                  'Plasmid_11','Plasmid_12','Plasmid_7','Plasmid_8',
                                                  'Plasmid_13','Plasmid_14','Plasmid_15','Plasmid_16',
                                                  'Plasmid_17'))
data$stage <- factor(data$stage, levels = c('PS1','PS2','PS1_2'))

data %>%
  ggplot(aes(x = arm, y = relative_fitness)) +
  geom_boxplot(aes(fill = stage), outlier.shape = NA) +
  facet_wrap(.~orf_name, ncol = 6,
             labeller = labeller(orf_name = c('Plasmid_1' = 'FY4\nempty_KAN',
                                              'Plasmid_2' = 'FY4\nempty_HYG',
                                              'Plasmid_3' = 'FY4 yll\nempty_KAN',
                                              'Plasmid_4' = 'FY4 yll\nempty_HYG',
                                              'Plasmid_5' = 'FY4 met15\nempty_KAN',
                                              'Plasmid_6' = 'FY4 met15\nempty_HYG',
                                              'Plasmid_7' = 'FY4 met15 yll\nempty_KAN',
                                              'Plasmid_8' = 'FY4 met15 yll\nempty_HYG',
                                              'Plasmid_9' = 'FY4 yll\nyll_KAN',
                                              'Plasmid_10' = 'FY4 yll\nmet15_HYG',
                                              'Plasmid_11' = 'FY4 met15\nyll_KAN',
                                              'Plasmid_12' = 'FY4 met15\nmet15_HYG',
                                              'Plasmid_13' = 'FY4 met15 yll\nyll_KAN',
                                              'Plasmid_14' = 'FY4 met15 yll\nmet15_HYG',
                                              'Plasmid_15' = 'FY4 met15 yll\nyll_KAN + met15_HYG',
                                              'Plasmid_16' = 'FY4 met15 yll\nyll_KAN + yll_CEN_HYG',
                                              'Plasmid_17' = 'FY4 met15 yll\nempty_KAN + empty_HYG'))) +
  scale_x_discrete(limits = c('PV_FY_PM','PV_FY_MM'),
                   labels = c('+ Met', '- Met')) +
  labs(y = 'Colony Size (pixels)') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        axis.title.x = element_blank(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))

