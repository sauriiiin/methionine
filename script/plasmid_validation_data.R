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

tables <- read_excel(sprintf('%s/%s/PV_INFO.xlsx',out_path,expt.name))
tables <- tables[order(tables$expt_id, tables$condition, tables$arm),] %>% data.frame()
head(tables)
##### GATHER AND CLEAN DATA FROM SQL
data <- NULL
for (s in unique(tables$stage_id)) {
  for (a in unique(tables$arm[tables$stage_id == s])) {
    temp <- dbGetQuery(conn, sprintf('select a.pos, a.hours, a.average, b.density, b.plate, b.row, b.col, c.orf_name
                                     from PV_%s_%s_384_CLEAN a, PV_pos2coor b, PV_pos2orf_name c
                                     where a.pos = b.pos and b.pos = c.pos
                                     order by a.hours, b.plate, b.col, b.row', s, a))
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
      temp$relative_fitness[temp$hours == h & str_detect(temp$orf_name,'_CEN_')] <- temp$average[temp$hours == h & str_detect(temp$orf_name,'_CEN_')]/
        median(temp$average[temp$hours == h & temp$orf_name == 'met15del_CEN_empty'], na.rm = T)
      temp$relative_fitness[temp$hours == h & str_detect(temp$orf_name,'_2M_')] <- temp$average[temp$hours == h & str_detect(temp$orf_name,'_2M_')]/
        median(temp$average[temp$hours == h & temp$orf_name == 'met15del_2M_empty'], na.rm = T)
    }
    temp$stage <- s
    temp$expt_rep <- a
    temp$condition <- unique(tables$condition[tables$stage == s & tables$arm == a])

    data <- rbind(data,temp)
  }
}
data <- data.frame(data)
head(data)

temp <- cbind(str_split(unique(data$orf_name), '_', simplify = T), unique(data$orf_name))
temp[temp[,4] == "",4] <- temp[temp[,4] == "",3]
temp[temp[,4] == temp[,3],3] <- temp[temp[,4] == temp[,3],2]
temp[temp[,3] == temp[,2],2] <- ""
colnames(temp) <- c('deletion1','deletion2','plasmid_backbone','plasmid_orf','orf_name')
temp <- data.frame(temp)                        
temp$deletion1 <- str_remove(temp$deletion1, 'del')
temp$deletion2 <- str_remove(temp$deletion2, 'del')

data <- merge(data, temp, by= 'orf_name')

data$hours[data$stage == 'WC' & data$hours == 0] <- 44
data$hours[data$stage == 'WC' & data$hours == 4] <- 48

data.sum <- data %>%
  group_by(condition,expt_rep,orf_name,deletion1,deletion2,plasmid_backbone,plasmid_orf,bio_rep) %>%
  summarize(relative_fitness = median(relative_fitness, na.rm = T), .groups = 'keep') %>%
  data.frame()


##### ANOVA ANALYSIS
# anova.res <- NULL
# for (c in unique(data$condition)) {
#   for (s1 in unique(data$orf_name)) {
#     for (s2 in unique(data$orf_name[data$orf_name != s1])) {
#       res.rcs.aov <- data[data$condition == c & data$orf_name %in% c(s2, s1),] %>%
#         data.frame() %>%
#         anova_test(relative_fitness ~ orf_name * bio_rep)
#       
#       res.rcs.kw <- data.sum[data.sum$condition == c & data.sum$orf_name %in% c(s2, s1),] %>%
#         data.frame() %>%
#         kruskal_test(relative_fitness ~ orf_name)
#       
#       cs1 <- data$relative_fitness[data$condition == c & data$orf_name == s1]
#       cs2 <- data$relative_fitness[data$condition == c & data$orf_name == s2]
#       emp_effsize <- (median(cs2, na.rm = T) - median(cs1, na.rm = T))/median(cs1, na.rm = T)
#       
#       anova.res <- rbind(anova.res, cbind(c, s1, s2,
#                                           res.rcs.aov$p[res.rcs.aov$Effect == 'orf_name'],
#                                           res.rcs.aov$p[res.rcs.aov$Effect == 'orf_name:bio_rep'],
#                                           res.rcs.kw$p,
#                                           emp_effsize))
#     }
#   }
# }
# colnames(anova.res) <- c('condition','reference','query','rcs_between','rcs_within','kw','effect_size')
# anova.res <- data.frame(anova.res, stringsAsFactors = F)
# anova.res$rcs_between <- as.numeric(anova.res$rcs_between)
# anova.res$rcs_within <- as.numeric(anova.res$rcs_within)
# anova.res$kw <- as.numeric(anova.res$kw)
# anova.res$effect_size <- as.numeric(anova.res$effect_size)
# anova.res$rcs_between <- p.adjust(anova.res$rcs_between, method = 'BH')
# anova.res$rcs_within <- p.adjust(anova.res$rcs_within, method = 'BH')
# anova.res$kw <- p.adjust(anova.res$kw, method = 'BH')
# anova.res$label_aov[anova.res$rcs_between > 0.05] <- 'ns'
# anova.res$label_aov[anova.res$rcs_between <= 0.05] <- '*'
# anova.res$label_aov[anova.res$rcs_between <= 0.01] <- '**'
# anova.res$label_aov[anova.res$rcs_between <= 0.001] <- '***'
# anova.res$label_aov[anova.res$rcs_between <= 0.0001] <- '****'
# anova.res$label_kw[anova.res$kw > 0.05] <- 'ns'
# anova.res$label_kw[anova.res$kw <= 0.05] <- '*'
# anova.res$label_kw[anova.res$kw <= 0.01] <- '**'
# anova.res$label_kw[anova.res$kw <= 0.001] <- '***'
# anova.res$label_kw[anova.res$kw <= 0.0001] <- '****'

##### SAVE COLONY SIZE DATA
save(data, data.sum, file = sprintf("%s/%s/colonysizes.RData", out_path, expt.name))
# save(anova.res, file = sprintf("%s/%s/stats.RData", out_path, expt.name))

#####
# END
#####