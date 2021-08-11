##### NO SULFATE EXPERIMENT - DATA
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 08/04/2021 

##### INITIALIZE
library(RMariaDB)
library(readxl)
library(dplyr)
library(rstatix)

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

out_path <- "~/R/Projects/methionine/data"
expt.name <- "nosulfate"

tables <- read_excel(sprintf('%s/%s/NO_SUL_INFO.xlsx',out_path,expt.name))
tables <- tables[order(tables$expt_id, tables$condition, tables$arm),] %>% filter(expt_id == 'SUL')
head(tables)

##### GATHER AND CLEAN DATA FROM SQL
data <- NULL
for (s in unique(tables$stage_id)) {
  if (s == 'S3') {
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

    data <- rbind(data,temp)
  }
}
data <- data.frame(data)
head(data)

data.sum <- data %>%
  group_by(condition,stage,arm,orf_name,bio_rep) %>%
  summarize(cs = median(average, na.rm = T), 
            relative_fitness = median(relative_fitness, na.rm = T), .groups = 'keep') %>%
  data.frame()

##### ANOVA ANALYSIS
anova.res <- NULL
for (s in unique(data$stage)) {
  for (c in unique(data$condition[data$stage == s])) {
    for (s1 in unique(data$orf_name[data$stage == s & data$condition == c])) {
      for (s2 in unique(data$orf_name[data$stage == s & data$condition == c & data$orf_name != s1])) {
        res.rcs.aov <- data[data$stage == s & data$condition == c & data$orf_name %in% c(s2, s1),] %>%
          data.frame() %>%
          anova_test(relative_fitness ~ orf_name * bio_rep)
        
        res.rcs.kw <- data.sum[data.sum$stage == s & data.sum$condition == c & data.sum$orf_name %in% c(s2, s1),] %>%
          data.frame() %>%
          kruskal_test(relative_fitness ~ orf_name)
        
        cs1 <- data$relative_fitness[data$stage == s & data$condition == c & data$orf_name == s1]
        cs2 <- data$relative_fitness[data$stage == s & data$condition == c & data$orf_name == s2]
        emp_effsize <- (median(cs2, na.rm = T) - median(cs1, na.rm = T))/median(cs1, na.rm = T)
        
        anova.res <- rbind(anova.res, cbind(s, c, s1, s2,
                                            res.rcs.aov$p[res.rcs.aov$Effect == 'orf_name'],
                                            res.rcs.aov$p[res.rcs.aov$Effect == 'orf_name:bio_rep'],
                                            res.rcs.kw$p,
                                            emp_effsize))
      }
    }
  }
}
colnames(anova.res) <- c('stage','condition','reference','query','rcs_between','rcs_within','kw','effect_size')
anova.res <- data.frame(anova.res, stringsAsFactors = F)
anova.res$rcs_between <- as.numeric(anova.res$rcs_between)
anova.res$rcs_within <- as.numeric(anova.res$rcs_within)
anova.res$kw <- as.numeric(anova.res$kw)
anova.res$effect_size <- as.numeric(anova.res$effect_size)
anova.res$rcs_between <- p.adjust(anova.res$rcs_between, method = 'BH')
anova.res$rcs_within <- p.adjust(anova.res$rcs_within, method = 'BH')
anova.res$kw <- p.adjust(anova.res$kw, method = 'BH')
anova.res$label_aov[anova.res$rcs_between > 0.05] <- 'ns'
anova.res$label_aov[anova.res$rcs_between <= 0.05] <- '*'
anova.res$label_aov[anova.res$rcs_between <= 0.01] <- '**'
anova.res$label_aov[anova.res$rcs_between <= 0.001] <- '***'
anova.res$label_aov[anova.res$rcs_between <= 0.0001] <- '****'
anova.res$label_kw[anova.res$kw > 0.05] <- 'ns'
anova.res$label_kw[anova.res$kw <= 0.05] <- '*'
anova.res$label_kw[anova.res$kw <= 0.01] <- '**'
anova.res$label_kw[anova.res$kw <= 0.001] <- '***'
anova.res$label_kw[anova.res$kw <= 0.0001] <- '****'

##### SAVE COLONY SIZE DATA
save(data, data.sum, file = sprintf("%s/%s/colonysizes.RData", out_path, expt.name))
save(anova.res, file = sprintf("%s/%s/stats.RData", out_path, expt.name))

#####
# END
#####