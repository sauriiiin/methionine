##### FINAL STATS FOR METHIONINE PAPER
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 09/24/2021 

##### INITIALIZE
library(dplyr)
library(gtools)
library(ggpubr)
library(RMariaDB)

source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

load(file = 'figures/final/data.RData')

##### CARBON SOURCE
data.cbn$id <- paste(data.cbn$condition, data.cbn$orf_name, sep = '_')
data.cbn.pmt <- permutations(length(unique(data.cbn$id)),2,unique(data.cbn$id))

data.cbn.stats <- NULL
for (i in seq(1,dim(data.cbn.pmt)[1])) {

  temp <- data.cbn[data.cbn$id %in% data.cbn.pmt[i,],] %>%
    group_by(id, expt_rep, bio_rep) %>%
    summarise(relative_fitness = median(relative_fitness, na.rm = T), .groups = 'keep') %>%
    data.frame()
  temp.p <- kruskal.test(relative_fitness ~ id, data = temp)
  
  data.cbn.stats$orf_name_1[i] <- as.character(unique(data.cbn$orf_name[data.cbn$id == data.cbn.pmt[i,1]]))
  data.cbn.stats$base_1[i] <- as.character(unique(data.cbn$base[data.cbn$id == data.cbn.pmt[i,1]]))
  data.cbn.stats$methionine_1[i] <- as.character(unique(data.cbn$methionine[data.cbn$id == data.cbn.pmt[i,1]]))
  data.cbn.stats$uracil_1[i] <- as.character(unique(data.cbn$uracil[data.cbn$id == data.cbn.pmt[i,1]]))
  data.cbn.stats$cysteine_1[i] <- as.character(unique(data.cbn$cysteine[data.cbn$id == data.cbn.pmt[i,1]]))
  data.cbn.stats$carbon_1[i] <- as.character(unique(data.cbn$carbon[data.cbn$id == data.cbn.pmt[i,1]]))
  data.cbn.stats$met_aux_1[i] <- as.character(unique(data.cbn$met_aux[data.cbn$id == data.cbn.pmt[i,1]]))
  data.cbn.stats$ura_aux_1[i] <- as.character(unique(data.cbn$ura_aux[data.cbn$id == data.cbn.pmt[i,1]]))
  data.cbn.stats$auxotrophy_1[i] <- as.character(unique(data.cbn$auxotrophy[data.cbn$id == data.cbn.pmt[i,1]]))
  data.cbn.stats$labels_1[i] <- as.character(unique(data.cbn$labels[data.cbn$id == data.cbn.pmt[i,1]]))
  
  data.cbn.stats$orf_name_2[i] <- as.character(unique(data.cbn$orf_name[data.cbn$id == data.cbn.pmt[i,2]]))
  data.cbn.stats$base_2[i] <- as.character(unique(data.cbn$base[data.cbn$id == data.cbn.pmt[i,2]]))
  data.cbn.stats$methionine_2[i] <- as.character(unique(data.cbn$methionine[data.cbn$id == data.cbn.pmt[i,2]]))
  data.cbn.stats$uracil_2[i] <- as.character(unique(data.cbn$uracil[data.cbn$id == data.cbn.pmt[i,2]]))
  data.cbn.stats$cysteine_2[i] <- as.character(unique(data.cbn$cysteine[data.cbn$id == data.cbn.pmt[i,2]]))
  data.cbn.stats$carbon_2[i] <- as.character(unique(data.cbn$carbon[data.cbn$id == data.cbn.pmt[i,2]]))
  data.cbn.stats$met_aux_2[i] <- as.character(unique(data.cbn$met_aux[data.cbn$id == data.cbn.pmt[i,2]]))
  data.cbn.stats$ura_aux_2[i] <- as.character(unique(data.cbn$ura_aux[data.cbn$id == data.cbn.pmt[i,2]]))
  data.cbn.stats$auxotrophy_2[i] <- as.character(unique(data.cbn$auxotrophy[data.cbn$id == data.cbn.pmt[i,2]]))
  data.cbn.stats$labels_2[i] <- as.character(unique(data.cbn$labels[data.cbn$id == data.cbn.pmt[i,2]]))
  
  data.cbn.stats$method[i] <- temp.p$method
  data.cbn.stats$statistic[i] <- temp.p$statistic
  data.cbn.stats$p[i] <- temp.p$p.value
  
  temp <- data.cbn[data.cbn$id %in% data.cbn.pmt[i,],] %>%
    group_by(id) %>%
    summarise(relative_fitness = median(relative_fitness, na.rm = T), .groups = 'keep') %>%
    data.frame()
  
  data.cbn.stats$effsize[i] <- (temp$relative_fitness[temp$id == data.cbn.pmt[i,2]] - temp$relative_fitness[temp$id == data.cbn.pmt[i,1]])/
    temp$relative_fitness[temp$id == data.cbn.pmt[i,1]]
}
data.cbn.stats <- data.frame(data.cbn.stats)
# dbWriteTable(conn, 'MET_CARS_STATS', data.cbn.stats, overwrite = T)
write.csv(data.cbn.stats, file = '/home/sbp29/R/Projects/methionine/figures/final/final/MET_CARS_STATS.csv')


##### JAKES MUTANTS
data.jm$id <- paste(data.jm$condition, data.jm$orf_name, sep = '_')
data.jm.pmt <- permutations(length(unique(data.jm$id)),2,unique(data.jm$id))

data.jm.stats <- NULL
for (i in seq(1,dim(data.jm.pmt)[1])) {
  
  temp <- data.jm[data.jm$id %in% data.jm.pmt[i,],] %>%
    filter(time == 't_final', attempt != 'pilot') %>%
    group_by(id, attempt_label, bio_rep) %>%
    summarise(relative_fitness = median(relative_cs, na.rm = T), .groups = 'keep') %>%
    data.frame()
  temp.p <- kruskal.test(relative_fitness ~ id, data = temp)
  
  data.jm.stats$orf_name_1[i] <- as.character(unique(data.jm$orf_name[data.jm$id == data.jm.pmt[i,1]]))
  data.jm.stats$condition_1[i] <- as.character(unique(data.jm$condition[data.jm$id == data.jm.pmt[i,1]]))
  data.jm.stats$orf_name_2[i] <- as.character(unique(data.jm$orf_name[data.jm$id == data.jm.pmt[i,2]]))
  data.jm.stats$condition_2[i] <- as.character(unique(data.jm$condition[data.jm$id == data.jm.pmt[i,2]]))

  
  data.jm.stats$method[i] <- temp.p$method
  data.jm.stats$statistic[i] <- temp.p$statistic
  data.jm.stats$p[i] <- temp.p$p.value
  
  temp <- data.jm[data.jm$id %in% data.jm.pmt[i,],] %>%
    filter(time == 't_final', attempt != 'pilot') %>%
    group_by(id) %>%
    summarise(relative_fitness = median(relative_cs, na.rm = T), .groups = 'keep') %>%
    data.frame()
  
  data.jm.stats$effsize[i] <- (temp$relative_fitness[temp$id == data.jm.pmt[i,2]] - temp$relative_fitness[temp$id == data.jm.pmt[i,1]])/
    temp$relative_fitness[temp$id == data.jm.pmt[i,1]]
}
data.jm.stats <- data.frame(data.jm.stats)
write.csv(data.jm.stats, file = '/home/sbp29/R/Projects/methionine/figures/final/final/MET_JM_STATS.csv')


##### RESPIRATION
data.res.gc.sum$id <- paste(data.res.gc.sum$condition, data.res.gc.sum$orf_name, sep = '_')
data.res.gc.sum.pmt <- permutations(length(unique(data.res.gc.sum$id)),2,unique(data.res.gc.sum$id))

data.res.gc.sum.stats <- NULL
for (i in seq(1,dim(data.res.gc.sum.pmt)[1])) {
  
  temp <- data.res.gc.sum[data.res.gc.sum$id %in% data.res.gc.sum.pmt[i,],] %>%
    group_by(id, expt_rep, bio_rep) %>%
    summarise(relative_fitness = median(rel_auc, na.rm = T), .groups = 'keep') %>%
    data.frame()
  temp.p <- kruskal.test(relative_fitness ~ id, data = temp)
  
  data.res.gc.sum.stats$orf_name_1[i] <- as.character(unique(data.res.gc.sum$orf_name[data.res.gc.sum$id == data.res.gc.sum.pmt[i,1]]))
  data.res.gc.sum.stats$condition_1[i] <- as.character(unique(data.res.gc.sum$condition[data.res.gc.sum$id == data.res.gc.sum.pmt[i,1]]))
  data.res.gc.sum.stats$base_1[i] <- as.character(unique(data.res.gc.sum$base[data.res.gc.sum$id == data.res.gc.sum.pmt[i,1]]))
  data.res.gc.sum.stats$methionine_1[i] <- as.character(unique(data.res.gc.sum$methionine[data.res.gc.sum$id == data.res.gc.sum.pmt[i,1]]))
  data.res.gc.sum.stats$cysteine_1[i] <- as.character(unique(data.res.gc.sum$cysteine[data.res.gc.sum$id == data.res.gc.sum.pmt[i,1]]))
  data.res.gc.sum.stats$carbon_1[i] <- as.character(unique(data.res.gc.sum$carbon[data.res.gc.sum$id == data.res.gc.sum.pmt[i,1]]))
  
  data.res.gc.sum.stats$orf_name_2[i] <- as.character(unique(data.res.gc.sum$orf_name[data.res.gc.sum$id == data.res.gc.sum.pmt[i,2]]))
  data.res.gc.sum.stats$condition_2[i] <- as.character(unique(data.res.gc.sum$condition[data.res.gc.sum$id == data.res.gc.sum.pmt[i,2]]))
  data.res.gc.sum.stats$base_2[i] <- as.character(unique(data.res.gc.sum$base[data.res.gc.sum$id == data.res.gc.sum.pmt[i,2]]))
  data.res.gc.sum.stats$methionine_2[i] <- as.character(unique(data.res.gc.sum$methionine[data.res.gc.sum$id == data.res.gc.sum.pmt[i,2]]))
  data.res.gc.sum.stats$cysteine_2[i] <- as.character(unique(data.res.gc.sum$cysteine[data.res.gc.sum$id == data.res.gc.sum.pmt[i,2]]))
  data.res.gc.sum.stats$carbon_2[i] <- as.character(unique(data.res.gc.sum$carbon[data.res.gc.sum$id == data.res.gc.sum.pmt[i,2]]))
  
  data.res.gc.sum.stats$method[i] <- temp.p$method
  data.res.gc.sum.stats$statistic[i] <- temp.p$statistic
  data.res.gc.sum.stats$p[i] <- temp.p$p.value
  
  temp <- data.res.gc.sum[data.res.gc.sum$id %in% data.res.gc.sum.pmt[i,],] %>%
    group_by(id) %>%
    summarise(relative_fitness = median(rel_auc, na.rm = T), .groups = 'keep') %>%
    data.frame()
  
  data.res.gc.sum.stats$effsize[i] <- (temp$relative_fitness[temp$id == data.res.gc.sum.pmt[i,2]] - temp$relative_fitness[temp$id == data.res.gc.sum.pmt[i,1]])/
    temp$relative_fitness[temp$id == data.res.gc.sum.pmt[i,1]]
}
data.res.gc.sum.stats <- data.frame(data.res.gc.sum.stats)
# write.csv(data.res.gc.sum.stats, file = '/home/sbp29/R/Projects/methionine/figures/final/final/MET_RESP_STATS.csv')

stats.res.supp <- data.res.gc.sum.stats[(data.res.gc.sum.stats$condition_1 == data.res.gc.sum.stats$condition_2 & data.res.gc.sum.stats$orf_name_1 == 'FY4') |
                        (data.res.gc.sum.stats$carbon_1 == data.res.gc.sum.stats$carbon_2 &
                           data.res.gc.sum.stats$methionine_1 == '+Met' & data.res.gc.sum.stats$methionine_2 == '-Met' &
                           data.res.gc.sum.stats$orf_name_1 == data.res.gc.sum.stats$orf_name_2),]
write.csv(stats.res.supp, file = '/home/sbp29/R/Projects/methionine/figures/final/final/RES_STATS_SUPP.csv')

##### PLASMID VALIDATION
head(data.sum.pv)
data.sum.pv <- merge(data.sum.pv, data.sum.pv %>%
        group_by(stage) %>%
        summarize(saturation = max(hours)) %>% data.frame(), by = 'stage')

data.sum.pv$id <- paste(data.sum.pv$stage, data.sum.pv$condition, data.sum.pv$orf_name, sep = '_')
data.sum.pv.pmt <- permutations(length(unique(data.sum.pv$id)),2,unique(data.sum.pv$id))

data.sum.pv.stats <- NULL
for (i in seq(1,dim(data.sum.pv.pmt)[1])) {
  
  temp <- data.sum.pv[data.sum.pv$id %in% data.sum.pv.pmt[i,],] %>%
    filter(hours == saturation) %>%
    group_by(id, expt_rep, bio_rep) %>%
    summarise(relative_fitness = median(relative_fitness, na.rm = T),
              cs = median(cs, na.rm = T), .groups = 'keep') %>%
    data.frame()
  temp$relative_fitness[is.na(temp$relative_fitness)] <- 0
  temp$cs[is.na(temp$cs)] <- 0
  
  temp.p <- kruskal.test(relative_fitness ~ id, data = temp)
  temp.p2 <- kruskal.test(cs ~ id, data = temp)
  
  data.sum.pv.stats$orf_name_1[i] <- as.character(unique(data.sum.pv$orf_name[data.sum.pv$id == data.sum.pv.pmt[i,1]]))
  data.sum.pv.stats$stage_1[i] <- as.character(unique(data.sum.pv$stage[data.sum.pv$id == data.sum.pv.pmt[i,1]]))
  data.sum.pv.stats$condition_1[i] <- as.character(unique(data.sum.pv$condition[data.sum.pv$id == data.sum.pv.pmt[i,1]]))
  data.sum.pv.stats$deletion1_1[i] <- as.character(unique(data.sum.pv$deletion1[data.sum.pv$id == data.sum.pv.pmt[i,1]]))
  data.sum.pv.stats$deletion2_1[i] <- as.character(unique(data.sum.pv$deletion2[data.sum.pv$id == data.sum.pv.pmt[i,1]]))
  data.sum.pv.stats$plasmid_backbone_1[i] <- as.character(unique(data.sum.pv$plasmid_backbone[data.sum.pv$id == data.sum.pv.pmt[i,1]]))
  data.sum.pv.stats$plasmid_orf_1[i] <- as.character(unique(data.sum.pv$plasmid_orf[data.sum.pv$id == data.sum.pv.pmt[i,1]]))
  
  data.sum.pv.stats$orf_name_2[i] <- as.character(unique(data.sum.pv$orf_name[data.sum.pv$id == data.sum.pv.pmt[i,2]]))
  data.sum.pv.stats$stage_2[i] <- as.character(unique(data.sum.pv$stage[data.sum.pv$id == data.sum.pv.pmt[i,2]]))
  data.sum.pv.stats$condition_2[i] <- as.character(unique(data.sum.pv$condition[data.sum.pv$id == data.sum.pv.pmt[i,2]]))
  data.sum.pv.stats$deletion1_2[i] <- as.character(unique(data.sum.pv$deletion1[data.sum.pv$id == data.sum.pv.pmt[i,2]]))
  data.sum.pv.stats$deletion2_2[i] <- as.character(unique(data.sum.pv$deletion2[data.sum.pv$id == data.sum.pv.pmt[i,2]]))
  data.sum.pv.stats$plasmid_backbone_2[i] <- as.character(unique(data.sum.pv$plasmid_backbone[data.sum.pv$id == data.sum.pv.pmt[i,2]]))
  data.sum.pv.stats$plasmid_orf_2[i] <- as.character(unique(data.sum.pv$plasmid_orf[data.sum.pv$id == data.sum.pv.pmt[i,2]]))
  
  data.sum.pv.stats$method[i] <- temp.p$method
  data.sum.pv.stats$statistic_fit[i] <- temp.p$statistic
  data.sum.pv.stats$p_fit[i] <- temp.p$p.value
  
  data.sum.pv.stats$statistic_cs[i] <- temp.p2$statistic
  data.sum.pv.stats$p_cs[i] <- temp.p2$p.value
  
  temp <- data.sum.pv[data.sum.pv$id %in% data.sum.pv.pmt[i,],] %>%
    group_by(id) %>%
    summarise(relative_fitness = median(relative_fitness, na.rm = T),
              cs = median(cs, na.rm = T), .groups = 'keep') %>%
    data.frame()
  
  data.sum.pv.stats$effsize_fit[i] <- (temp$relative_fitness[temp$id == data.sum.pv.pmt[i,2]] - temp$relative_fitness[temp$id == data.sum.pv.pmt[i,1]])/
    temp$relative_fitness[temp$id == data.sum.pv.pmt[i,1]]
  data.sum.pv.stats$effsize_cs[i] <- (temp$cs[temp$id == data.sum.pv.pmt[i,2]] - temp$cs[temp$id == data.sum.pv.pmt[i,1]])/
    temp$cs[temp$id == data.sum.pv.pmt[i,1]]
}
data.sum.pv.stats <- data.frame(data.sum.pv.stats)
write.csv(data.sum.pv.stats, file = '/home/sbp29/R/Projects/methionine/figures/final/final/MET_PV_STATS.csv')

##### PLASMID VALIDATION 2
head(data.sum.pv)
head(data.pv2)

data.pv %>%
  filter(stage == 'PS1', hours %in% c(48,115),
         orf_name %in% c('met15del_met12del_2M_empty','met15del_met12del_2M_met12')) %>%
  group_by(orf_name) %>%
  summarize(rel_f = median(relative_fitness, na.rm = T))

data.pv2 %>%
  filter(orf_name %in% c('Plasmid_5','Plasmid_11'),
         arm == 'PV_FY_MM', stage != 'PS2') %>%
  group_by(orf_name) %>%
  summarize(rel_f = median(relative_fitness, na.rm = T))

##### NO SULFATE
data.ns$id <- paste(data.ns$condition, data.ns$orf_name, sep = '_')
data.ns.pmt <- permutations(length(unique(data.ns$id)),2,unique(data.ns$id))

data.ns <- merge(data.ns, data.ns %>%
                   group_by(stage) %>%
                   summarize(.groups = 'keep', max_hrs = max(hours, na.rm = T)) %>%
                   data.frame(), by = 'stage')

data.ns.stats <- NULL
for (i in seq(1,dim(data.ns.pmt)[1])) {
  
  temp <- data.ns[data.ns$id %in% data.ns.pmt[i,],] %>%
    filter(hours == max_hrs, stage == 'S1', average != 0) %>%
    group_by(id, expt_rep, bio_rep) %>%
    summarise(cs = median(average, na.rm = T), .groups = 'keep') %>%
    data.frame()
  temp.p <- kruskal.test(cs ~ id, data = temp)
  
  data.ns.stats$orf_name_1[i] <- as.character(unique(data.ns$orf_name[data.ns$id == data.ns.pmt[i,1]]))
  data.ns.stats$condition_1[i] <- as.character(unique(data.ns$condition[data.ns$id == data.ns.pmt[i,1]]))
  data.ns.stats$media_1[i] <- as.character(unique(data.ns$media[data.ns$id == data.ns.pmt[i,1]]))
  data.ns.stats$base_1[i] <- as.character(unique(data.ns$base[data.ns$id == data.ns.pmt[i,1]]))
  data.ns.stats$ynb_type_1[i] <- as.character(unique(data.ns$ynb_type[data.ns$id == data.ns.pmt[i,1]]))
  data.ns.stats$methionine_1[i] <- as.character(unique(data.ns$methionine[data.ns$id == data.ns.pmt[i,1]]))
  data.ns.stats$sulfate_1[i] <- as.character(unique(data.ns$sulfate[data.ns$id == data.ns.pmt[i,1]]))
  
  data.ns.stats$orf_name_2[i] <- as.character(unique(data.ns$orf_name[data.ns$id == data.ns.pmt[i,2]]))
  data.ns.stats$condition_2[i] <- as.character(unique(data.ns$condition[data.ns$id == data.ns.pmt[i,2]]))
  data.ns.stats$media_2[i] <- as.character(unique(data.ns$media[data.ns$id == data.ns.pmt[i,2]]))
  data.ns.stats$base_2[i] <- as.character(unique(data.ns$base[data.ns$id == data.ns.pmt[i,2]]))
  data.ns.stats$ynb_type_2[i] <- as.character(unique(data.ns$ynb_type[data.ns$id == data.ns.pmt[i,2]]))
  data.ns.stats$methionine_2[i] <- as.character(unique(data.ns$methionine[data.ns$id == data.ns.pmt[i,2]]))
  data.ns.stats$sulfate_2[i] <- as.character(unique(data.ns$sulfate[data.ns$id == data.ns.pmt[i,2]]))
  
  data.ns.stats$method[i] <- temp.p$method
  data.ns.stats$statistic[i] <- temp.p$statistic[[1]]
  data.ns.stats$p[i] <- temp.p$p.value
  
  temp <- data.ns[data.ns$id %in% data.ns.pmt[i,],] %>%
    filter(hours == max_hrs, stage == 'S1', average != 0) %>%
    group_by(id) %>%
    summarise(cs = median(average, na.rm = T), .groups = 'keep') %>%
    data.frame()
  
  data.ns.stats$effsize[i] <- (temp$cs[temp$id == data.ns.pmt[i,2]] - temp$cs[temp$id == data.ns.pmt[i,1]])/
    temp$cs[temp$id == data.ns.pmt[i,1]]
}
data.ns.stats <- data.frame(data.ns.stats)
write.csv(data.ns.stats, file = '/home/sbp29/R/Projects/methionine/figures/final/final/MET_NS_STATS.csv')


##### REPEATED PINNING
data.rp$id <- paste(data.rp$table_name, data.rp$orf_name, sep = '_')
data.rp.pmt <- permutations(length(unique(data.rp$id)),2,unique(data.rp$id))

data.rp.stats <- NULL
for (i in seq(1,dim(data.rp.pmt)[1])) {
  temp <- data.rp[data.rp$id %in% data.rp.pmt[i,],] %>%
    filter(hours == max_hrs, average != 0) %>%
    group_by(id) %>%
    # summarise(cs = median(average, na.rm = T), .groups = 'keep') %>%
    data.frame()
  temp.p <- kruskal.test(average ~ id, data = temp)
  
  data.rp.stats$orf_name_1[i] <- as.character(unique(data.rp$orf_name[data.rp$id == data.rp.pmt[i,1]]))
  data.rp.stats$condition_1[i] <- as.character(unique(data.rp$condition[data.rp$id == data.rp.pmt[i,1]]))
  data.rp.stats$media_1[i] <- as.character(unique(data.rp$media[data.rp$id == data.rp.pmt[i,1]]))
  data.rp.stats$base_1[i] <- as.character(unique(data.rp$base[data.rp$id == data.rp.pmt[i,1]]))
  data.rp.stats$ynb_type_1[i] <- as.character(unique(data.rp$ynb_type[data.rp$id == data.rp.pmt[i,1]]))
  data.rp.stats$methionine_1[i] <- as.character(unique(data.rp$methionine[data.rp$id == data.rp.pmt[i,1]]))
  data.rp.stats$sulfate_1[i] <- as.character(unique(data.rp$sulfate[data.rp$id == data.rp.pmt[i,1]]))
  
  data.rp.stats$orf_name_2[i] <- as.character(unique(data.rp$orf_name[data.rp$id == data.rp.pmt[i,2]]))
  data.rp.stats$condition_2[i] <- as.character(unique(data.rp$condition[data.rp$id == data.rp.pmt[i,2]]))
  data.rp.stats$media_2[i] <- as.character(unique(data.rp$media[data.rp$id == data.rp.pmt[i,2]]))
  data.rp.stats$base_2[i] <- as.character(unique(data.rp$base[data.rp$id == data.rp.pmt[i,2]]))
  data.rp.stats$ynb_type_2[i] <- as.character(unique(data.rp$ynb_type[data.rp$id == data.rp.pmt[i,2]]))
  data.rp.stats$methionine_2[i] <- as.character(unique(data.rp$methionine[data.rp$id == data.rp.pmt[i,2]]))
  data.rp.stats$sulfate_2[i] <- as.character(unique(data.rp$sulfate[data.rp$id == data.rp.pmt[i,2]]))
  
  data.rp.stats$method[i] <- temp.p$method
  data.rp.stats$statistic[i] <- temp.p$statistic[[1]]
  data.rp.stats$p[i] <- temp.p$p.value
  
  temp <- data.rp[data.rp$id %in% data.rp.pmt[i,],] %>%
    filter(hours == max_hrs, average != 0) %>%
    group_by(id) %>%
    summarise(cs = median(average, na.rm = T), .groups = 'keep') %>%
    data.frame()
  
  data.rp.stats$effsize[i] <- (temp$cs[temp$id == data.rp.pmt[i,2]] - temp$cs[temp$id == data.rp.pmt[i,1]])/
    temp$cs[temp$id == data.rp.pmt[i,1]]
}
data.rp.stats <- data.frame(data.rp.stats)
write.csv(data.rp.stats, file = '/home/sbp29/R/Projects/methionine/figures/final/final/MET_RP_STATS.csv')
