

##### CARBON SOURCE

head(data.cbn)
data.cbn.sum <- rbind(data.cbn %>%
                        group_by(condition, orf_name, expt_rep, bio_rep) %>%
                        summarize(f = median(relative_fitness, na.rm = T), .groups = 'keep') %>%
                        data.frame(),
                      data.cbn.leu %>%
                        group_by(condition, orf_name, expt_rep, bio_rep) %>%
                        summarize(f = median(relative_fitness, na.rm = T), .groups = 'keep') %>%
                        data.frame())
unique(data.cbn.sum$condition)

data.cbn.sum %>%
  group_by(condition, orf_name) %>%
  summarize(f = median(f, na.rm = T), .groups = 'keep') %>%
  filter(condition == 'MiM_Glu', orf_name %in% c('BY4742','BY4741')) %>%
  data.frame()

compare_means(f ~ orf_name, data.cbn.leu %>%
                group_by(condition, orf_name, expt_rep, bio_rep) %>%
                summarize(f = median(relative_fitness, na.rm = T), .groups = 'keep') %>%
                filter(condition == 'SCmLeu', orf_name %in% c('BY4742','BY4741')) %>%
                data.frame(), 
              method = "kruskal.test", paired = FALSE)

data.cbn.stat <- NULL
for (c in unique(data.cbn.sum$condition)) {
  for (o1 in unique(data.cbn.sum$orf_name[data.cbn.sum$condition == c])) {
    for (o2 in unique(data.cbn.sum$orf_name[data.cbn.sum$condition == c])) {
      if (o1 != o2) {
        temp.fit <- data.cbn.sum %>%
          group_by(condition, orf_name) %>%
          summarize(f = median(f, na.rm = T), .groups = 'keep') %>%
          filter(condition == c, orf_name %in% c(o1,o2)) %>%
          data.frame()
        
        temp.stat <- compare_means(f ~ orf_name, data.cbn.sum %>%
                                     group_by(condition, orf_name, expt_rep, bio_rep) %>%
                                     summarize(f = median(f, na.rm = T), .groups = 'keep') %>%
                                     filter(condition == c, orf_name %in% c(o1,o2)) %>%
                                     data.frame(), 
                                   method = "kruskal.test", paired = FALSE)
        
        temp.es <- (temp.fit$f[temp.fit$orf_name == o1] - temp.fit$f[temp.fit$orf_name == o2])/temp.fit$f[temp.fit$orf_name == o2]
        data.cbn.stat <- rbind(data.cbn.stat,
                               data.frame(expt = 'carbon_source', condition = c, 
                   reference = o2, mutant = o1,
                   effect_size = temp.es, p = temp.stat$p))
      
        
      }
    }
  }
}
data.cbn.stat$p.adj <- adjust_pvalue(data.cbn.stat$p, method = 'BH')
write.csv(data.cbn.stat, file = 'results/final/stats/carbon_source_orfwise.csv',row.names = F)

data.cbn.stat <- NULL
for (o in unique(data.cbn.sum$orf_name)) {
  for (c1 in unique(data.cbn.sum$condition[data.cbn.sum$orf_name == o])) {
    for (c2 in unique(data.cbn.sum$condition[data.cbn.sum$orf_name == o])) {
      if (c1 != c2) {
        temp.fit <- data.cbn.sum %>%
          group_by(condition, orf_name) %>%
          summarize(f = median(f, na.rm = T), .groups = 'keep') %>%
          filter(orf_name == o, condition %in% c(c1, c2)) %>%
          data.frame()
        
        temp.stat <- compare_means(f ~ condition, data.cbn.sum %>%
                                     group_by(condition, orf_name, expt_rep, bio_rep) %>%
                                     summarize(f = median(f, na.rm = T), .groups = 'keep') %>%
                                     filter(orf_name == o, condition %in% c(c1,c2)) %>%
                                     data.frame(), 
                                   method = "kruskal.test", paired = FALSE)
        
        temp.es <- (temp.fit$f[temp.fit$condition == c1] - temp.fit$f[temp.fit$condition == c2])/temp.fit$f[temp.fit$condition == c2]
        data.cbn.stat <- rbind(data.cbn.stat,
                               data.frame(expt = 'carbon_source', orf_name = o, 
                                          reference = c2, stress = c1,
                                          effect_size = temp.es, p = temp.stat$p))
        
        
      }
    }
  }
}
data.cbn.stat$p.adj <- adjust_pvalue(data.cbn.stat$p, method = 'BH')
write.csv(data.cbn.stat, file = 'results/final/stats/carbon_source_conditionwise.csv',row.names = F)


##### JAKES MUTANTS
head(data.jm.2)
data.jm.sum <- data.jm.2[!is.na(data.jm.2$orf_name),] %>%
  group_by(condition, orf_name, attempt, bio_rep) %>%
  summarize(f = median(relative_fitness, na.rm = T), .groups = 'keep') %>%
  data.frame()
data.jm.sum$attempt[data.jm.sum$attempt == 'attempt 1'] <- 'copy1'
unique(data.jm.sum$condition)

data.jm.stat <- NULL
for (c in unique(data.jm.sum$condition)) {
  for (o1 in unique(data.jm.sum$orf_name[data.jm.sum$condition == c])) {
    for (o2 in unique(data.jm.sum$orf_name[data.jm.sum$condition == c])) {
      if (o1 != o2) {
        temp.fit <- data.jm.sum %>%
          group_by(condition, orf_name) %>%
          summarize(f = median(f, na.rm = T), .groups = 'keep') %>%
          filter(condition == c, orf_name %in% c(o1,o2)) %>%
          data.frame()
        
        temp.stat <- compare_means(f ~ orf_name, data.jm.sum %>%
                                     group_by(condition, orf_name, attempt, bio_rep) %>%
                                     summarize(f = median(f, na.rm = T), .groups = 'keep') %>%
                                     filter(condition == c, orf_name %in% c(o1,o2)) %>%
                                     data.frame(), 
                                   method = "kruskal.test", paired = FALSE)
        
        temp.es <- (temp.fit$f[temp.fit$orf_name == o1] - temp.fit$f[temp.fit$orf_name == o2])/temp.fit$f[temp.fit$orf_name == o2]
        data.jm.stat <- rbind(data.jm.stat,
                               data.frame(expt = 'jakes_mutants', condition = c, 
                                          reference = o2, mutant = o1,
                                          effect_size = temp.es, p = temp.stat$p))
        
        
      }
    }
  }
}
data.jm.stat$p.adj <- adjust_pvalue(data.jm.stat$p, method = 'BH')
write.csv(data.jm.stat, file = 'results/final/stats/jakes_mutants_orfwise.csv',row.names = F)


data.jm.stat <- NULL
for (o in unique(data.jm.sum$orf_name)) {
  for (c1 in unique(data.jm.sum$condition[data.jm.sum$orf_name == o])) {
    for (c2 in unique(data.jm.sum$condition[data.jm.sum$orf_name == o])) {
      if (c1 != c2) {
        temp.fit <- data.jm.sum %>%
          group_by(condition, orf_name) %>%
          summarize(f = median(f, na.rm = T), .groups = 'keep') %>%
          filter(orf_name == o, condition %in% c(c1, c2)) %>%
          data.frame()
        
        temp.stat <- compare_means(f ~ condition, data.jm.sum %>%
                                     group_by(condition, orf_name, attempt, bio_rep) %>%
                                     summarize(f = median(f, na.rm = T), .groups = 'keep') %>%
                                     filter(orf_name == o, condition %in% c(c1,c2)) %>%
                                     data.frame(), 
                                   method = "kruskal.test", paired = FALSE)
        
        temp.es <- (temp.fit$f[temp.fit$condition == c1] - temp.fit$f[temp.fit$condition == c2])/temp.fit$f[temp.fit$condition == c2]
        data.jm.stat <- rbind(data.jm.stat,
                               data.frame(expt = 'jakes_mutants', orf_name = o, 
                                          reference = c2, stress = c1,
                                          effect_size = temp.es, p = temp.stat$p))
        
        
      }
    }
  }
}
data.jm.stat$p.adj <- adjust_pvalue(data.jm.stat$p, method = 'BH')
write.csv(data.jm.stat, file = 'results/final/stats/jakes_mutants_conditionwise.csv',row.names = F)


##### H2S FLASK
merge(data.pred.h2s %>%
        melt(id.vars = 'Time', variable.name = 'Flask', value.name = 'OD'), 
      strain.labs.h2s, by.x = 'Flask', by.y = 'Falsk') %>%
  filter(Time == 3075, Expected.Initial.OD == 0.1, Replicate != 'Rep_1') %>%
  group_by(FeEDTA, Strain, Replicate) %>%
  summarize(OD = median(OD, na.rm = T), .groups = 'keep')

data.h2s.re[data.h2s.re$Time..Minutes. == 3075 &
              !(data.h2s.re$FeEDTA == 'Absent' & data.h2s.re$Replicate == 'R2') &
              !(data.h2s.re$FeEDTA == 'Present' & data.h2s.re$Replicate == 'R3'),-2] %>%
  group_by(FeEDTA, ORF, Replicate) %>%
  summarize(OD = median(OD, na.rm = T), .groups = 'keep') 


##### PV STATS
data.pv2.temp <- data.pv2 %>%
  filter(orf_name %in% c('Plasmid_1','Plasmid_3','Plasmid_5','Plasmid_7','Plasmid_11','Plasmid_13'),
         arm == 'PV_FY_MM', stage == 'PS1_2', hours == 336) %>%
  group_by(orf_name, expt_rep, bio_rep) %>%
  summarise(relative_fitness = median(relative_fitness, na.rm = T), .groups = 'keep') %>%
  data.frame()

data.pv2.stats <- NULL
for (o1 in unique(data.pv2.temp$orf_name)) {
  for (o2 in unique(data.pv2.temp$orf_name)) {
    if (o1 != o2) {
      temp.stat <- compare_means(relative_fitness ~ orf_name,
                                 data.pv2.temp %>%
                                   filter(orf_name %in% c(o1, o2)) %>%
                                   data.frame(),
                                 method = "kruskal.test", paired = FALSE)
      temp.es <- (median(data.pv2.temp$relative_fitness[data.pv2.temp$orf_name == o1], na.rm = T) - 
                    median(data.pv2.temp$relative_fitness[data.pv2.temp$orf_name == o2], na.rm = T))/
        median(data.pv2.temp$relative_fitness[data.pv2.temp$orf_name == o2], na.rm = T)
      data.pv2.stats <- rbind(data.pv2.stats,
                            data.frame(expt = 'plasmid_validation2', 
                                       ref_orf = o2,
                                       mut_orf = o1,
                                       effect_size = temp.es, p = temp.stat$p))
    }
  }
}
data.pv2.stats$p.adj <- adjust_pvalue(data.pv2.stats$p, method = 'BH')
write.csv(data.pv2.stats, file = 'results/final/stats/plasmid_validation2.csv',row.names = F)

