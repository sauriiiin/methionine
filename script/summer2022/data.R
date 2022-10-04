##### DATA FOR THE SUMMER SCREEN ANALYSIS
##### Screens done by Aaron Zhang, Alexis Berger and Brandon Garcia
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 07/28/2021 


##### INITIALIZE
source('/home/sbp29/R/Projects/methionine/paper/scripts/initialize.R')
fig_path <- '~/R/Projects/methionine/figures/summer2022/'

mysql_ids <- read.csv('/home/sbp29/R/Projects/methionine/data/summer2022/summer_screens_mysql.csv')

orf_cat_aaron <- read.csv(file = '/home/acwach/translatome2/oe_transient', stringsAsFactors = F)
orf_cat_aaron$orf_type <- "Unclassified"
orf_cat_aaron$orf_type[(orf_cat_aaron$is_transient + orf_cat_aaron$translated + orf_cat_aaron$is_candidate == 3)] <- "Transient"
orf_cat_aaron$orf_type[(orf_cat_aaron$is_conserved + orf_cat_aaron$translated + orf_cat_aaron$is_candidate == 3)] <- "Conserved"

p2c.table <- dbGetQuery(conn, 'select * from Brandon.bg_MET_OE2206_pos2coor order by density, plate_no, plate_col, plate_row')

conditions <- data.frame(condition_id = c('GA',
                                          'PM','MM',
                                          'PMEDTA','MMEDTA'),
                         condition = c('SC-Ura+Gal',
                                       'SD+Met-Cys-Ura+Gal',
                                       'SD-Met-Cys-Ura+Gal',
                                       'SD+Met-Cys-Ura+Gal+FeEDTA',
                                       'SD-Met-Cys-Ura+Gal+FeEDTA'))
cnd_limits <- c('SC-Ura+Gal',
                'SD+Met-Cys-Ura+Gal',
                'SD+Met-Cys-Ura+Gal+FeEDTA',
                'SD-Met-Cys-Ura+Gal',
                'SD-Met-Cys-Ura+Gal+FeEDTA')

##### GATHER DATA
data.fit <- NULL
data.stats <- NULL
for (i in 1:dim(mysql_ids)[1]) {
  temp.fit <- dbGetQuery(conn, sprintf('select * from %s.%s_%s_%s_%s_%s_%d_FITNESS',
                                       mysql_ids$database[i],
                                       mysql_ids$user_id[i],mysql_ids$expt_id[i],mysql_ids$stage[i],
                                       mysql_ids$condition_id[i],mysql_ids$bag[i],mysql_ids$density[i]))
  temp.fit$saturation <- max(unique(temp.fit$hours))
  temp.fit$stage <- mysql_ids$stage[i]
  temp.fit$condition_id <- mysql_ids$condition_id[i]
  temp.fit$bag <- mysql_ids$bag[i]
  data.fit <- rbind(data.fit, temp.fit)
  
  temp.fit.stats <- dbGetQuery(conn, sprintf('select a.*, b.stat, b.p
                                             from %s.%s_%s_%s_%s_%s_%d_FITNESS_STATS a, %s.%s_%s_%s_%s_%s_%d_PVALUE b
                                             where a.hours = b.hours and a.strain_id = b.strain_id',
                                       mysql_ids$database[i],
                                       mysql_ids$user_id[i],mysql_ids$expt_id[i],mysql_ids$stage[i],
                                       mysql_ids$condition_id[i],mysql_ids$bag[i],mysql_ids$density[i],
                                       mysql_ids$database[i],
                                       mysql_ids$user_id[i],mysql_ids$expt_id[i],mysql_ids$stage[i],
                                       mysql_ids$condition_id[i],mysql_ids$bag[i],mysql_ids$density[i]))
  temp.fit.stats$saturation <- max(unique(temp.fit.stats$hours))
  temp.fit.stats$stage <- mysql_ids$stage[i]
  temp.fit.stats$condition_id <- mysql_ids$condition_id[i]
  temp.fit.stats$bag <- mysql_ids$bag[i]
  data.stats <- rbind(data.stats, temp.fit.stats)
}
data.fit$rep <- as.numeric(str_trunc(as.character(data.fit$pos), 6, side = 'left', ellipsis = ''))
colnames(data.stats) <- str_replace_all(colnames(data.stats), 'cs', 'fitness')

data.fit <- merge(data.fit, conditions, by = 'condition_id')
data.stats <- merge(data.stats, conditions, by = 'condition_id')


head(data.fit)
head(data.stats)

sat.times <- data.stats %>%
  group_by(stage, condition, condition_id, bag) %>%
  summarise(saturation = unique(saturation), .groups = 'keep') %>%
  data.frame()

##### ASSIGN PHENOTYPES AND GET COUNTS
data.stats$phenotype[data.stats$p > 0.05] <- 'Neutral'
data.stats$phenotype[data.stats$p <= 0.05 & data.stats$stat > 0] <- 'Beneficial'
data.stats$phenotype[data.stats$p <= 0.05 & data.stats$stat < 0] <- 'Deleterious'

data.cnts.all <- merge(data.stats %>%
                         filter(hours == saturation) %>%
                         group_by(stage, condition, condition_id, bag, hours, phenotype) %>%
                         count() %>% data.frame(),
                       data.stats %>%
                         filter(hours == saturation) %>%
                         group_by(stage, condition, condition_id, bag, hours) %>%
                         count() %>% data.frame(),
                       by = c('stage','condition','condition_id','bag','hours'), suffixes = c('','_total'))
data.cnts.all$percentage <- data.cnts.all$n/data.cnts.all$n_total * 100

data.cnts.trn <- merge(data.stats %>%
                         filter(hours == saturation, orf_name %in% orf_cat_aaron$orf_name[orf_cat_aaron$orf_type == 'Transient']) %>%
                         group_by(stage, condition, condition_id, bag, hours, phenotype) %>%
                         count() %>% data.frame(),
                       data.stats %>%
                         filter(hours == saturation, orf_name %in% orf_cat_aaron$orf_name[orf_cat_aaron$orf_type == 'Transient']) %>%
                         group_by(stage, condition, condition_id, bag, hours) %>%
                         count() %>% data.frame(),
                       by = c('stage','condition','condition_id','bag','hours'), suffixes = c('','_total'))
data.cnts.trn$percentage <- data.cnts.trn$n/data.cnts.trn$n_total * 100

##### CLEAN DATA
data.sum <- data.fit %>%
  group_by(stage, condition_id, condition, bag, hours, saturation, strain_id, orf_name, rep) %>%
  summarise(fitness_median = median(fitness, na.rm = T), cs_median = median(average, na.rm = T),
            fitness_mad = mad(fitness, na.rm = T), cs_mad = mad(average, na.rm = T),
            .groups = 'keep') %>%
  data.frame()

data.fit.clean <- merge(data.fit, data.sum, by = c('stage','condition','condition_id','bag','hours','saturation','strain_id','orf_name','rep'))
data.fit.clean$fitness[data.fit.clean$fitness < (data.fit.clean$fitness_median - 2*data.fit.clean$fitness_mad) |
                         data.fit.clean$fitness > (data.fit.clean$fitness_median + 2*data.fit.clean$fitness_mad)] <- NA
data.fit.clean$average[data.fit.clean$average < (data.fit.clean$cs_median - 2*data.fit.clean$cs_mad) |
                         data.fit.clean$average > (data.fit.clean$cs_median + 2*data.fit.clean$cs_mad)] <- NA


data.sum <- data.fit.clean %>%
  group_by(stage, condition_id, condition, bag, hours, saturation, strain_id, orf_name, rep) %>%
  summarise(fitness_median = median(fitness, na.rm = T), cs_median = median(average, na.rm = T),
            .groups = 'keep') %>%
  data.frame()


##### GIANT MATRIX
q <- "select matrix.*,
FS1_GA_nobag_score + FS1_GA_inbag_score + FS1_PM_nobag_score + FS1_PM_inbag_score + FS1_MM_nobag_score + FS1_MM_inbag_score + FS1_PMEDTA_nobag_score + FS1_PMEDTA_inbag_score + FS1_MMEDTA_nobag_score + FS1_MMEDTA_inbag_score + 
FS3_GA_nobag_score + FS3_GA_inbag_score + FS3_PM_nobag_score + FS3_PM_inbag_score + FS3_MM_inbag_score + FS3_PMEDTA_nobag_score + FS3_PMEDTA_inbag_score + FS3_MMEDTA_inbag_score
as score
from

(select a.*,
b.FS1_GA_nobag_fitness, b.FS1_GA_nobag_phenotype, b.FS1_GA_nobag_score,
c.FS1_GA_inbag_fitness, c.FS1_GA_inbag_phenotype, c.FS1_GA_inbag_score,
d.FS1_PM_nobag_fitness, d.FS1_PM_nobag_phenotype, d.FS1_PM_nobag_score,
e.FS1_PM_inbag_fitness, e.FS1_PM_inbag_phenotype, e.FS1_PM_inbag_score,
f.FS1_MM_nobag_fitness, f.FS1_MM_nobag_phenotype, f.FS1_MM_nobag_score,
g.FS1_MM_inbag_fitness, g.FS1_MM_inbag_phenotype, g.FS1_MM_inbag_score,
h.FS1_PMEDTA_nobag_fitness, h.FS1_PMEDTA_nobag_phenotype, h.FS1_PMEDTA_nobag_score,
i.FS1_PMEDTA_inbag_fitness, i.FS1_PMEDTA_inbag_phenotype, i.FS1_PMEDTA_inbag_score,
j.FS1_MMEDTA_nobag_fitness, j.FS1_MMEDTA_nobag_phenotype, j.FS1_MMEDTA_nobag_score,
k.FS1_MMEDTA_inbag_fitness, k.FS1_MMEDTA_inbag_phenotype, k.FS1_MMEDTA_inbag_score,
l.FS3_GA_nobag_fitness, l.FS3_GA_nobag_phenotype, l.FS3_GA_nobag_score,
m.FS3_GA_inbag_fitness, m.FS3_GA_inbag_phenotype, m.FS3_GA_inbag_score,
n.FS3_PM_nobag_fitness, n.FS3_PM_nobag_phenotype, n.FS3_PM_nobag_score,
o.FS3_PM_inbag_fitness, o.FS3_PM_inbag_phenotype, o.FS3_PM_inbag_score,
p.FS3_MM_inbag_fitness, p.FS3_MM_inbag_phenotype, p.FS3_MM_inbag_score,
q.FS3_PMEDTA_nobag_fitness, q.FS3_PMEDTA_nobag_phenotype, q.FS3_PMEDTA_nobag_score,
r.FS3_PMEDTA_inbag_fitness, r.FS3_PMEDTA_inbag_phenotype, r.FS3_PMEDTA_inbag_score,
s.FS3_MMEDTA_inbag_fitness, s.FS3_MMEDTA_inbag_phenotype, s.FS3_MMEDTA_inbag_score

from
(select strain_id, orf_name from BARFLEX_SPACE_AGAR_180313
where strain_id > 0
union
select strain_id, orf_name from PROTOGENE_COLLECTION) a
left join

(select a.strain_id, a.orf_name, a.cs_median FS1_GA_nobag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS1_GA_nobag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS1_GA_nobag_score
from Brandon.bg_MET_OE2206_FS1_GA_nobag_6144_FITNESS_STATS a, Brandon.bg_MET_OE2206_FS1_GA_nobag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 24) b
on a.strain_id = b.strain_id and a.orf_name = b.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS1_GA_inbag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS1_GA_inbag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS1_GA_inbag_score
from alexis.anb_MET_OE2206_FS1_GA_inbag_6144_FITNESS_STATS a, alexis.anb_MET_OE2206_FS1_GA_inbag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and b.hours = 15) c
on a.strain_id = c.strain_id and a.orf_name = c.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS1_PM_nobag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS1_PM_nobag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS1_PM_nobag_score
from aaron_z.aaz_MET_OE2206_FS1_PM_nobag_6144_FITNESS_STATS a, aaron_z.aaz_MET_OE2206_FS1_PM_nobag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and b.hours = 24) d
on a.strain_id = d.strain_id and a.orf_name = d.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS1_PM_inbag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS1_PM_inbag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS1_PM_inbag_score
from alexis.anb_MET_OE2206_FS1_PM_inbag_6144_FITNESS_STATS a, alexis.anb_MET_OE2206_FS1_PM_inbag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and b.hours = 18) e
on a.strain_id = e.strain_id and a.orf_name = e.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS1_MM_nobag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS1_MM_nobag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS1_MM_nobag_score
from Brandon.bg_MET_OE2206_FS1_MM_nobag_6144_FITNESS_STATS a, Brandon.bg_MET_OE2206_FS1_MM_nobag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and b.hours = 147) f
on a.strain_id = f.strain_id and a.orf_name = f.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS1_MM_inbag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS1_MM_inbag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS1_MM_inbag_score
from Brandon.bg_MET_OE2206_FS1_MM_inbag_6144_FITNESS_STATS a, Brandon.bg_MET_OE2206_FS1_MM_inbag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and b.hours = 150) g
on a.strain_id = g.strain_id and a.orf_name = g.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS1_PMEDTA_nobag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS1_PMEDTA_nobag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS1_PMEDTA_nobag_score
from alexis.anb_MET_OE2206_FS1_PMEDTA_nobag_6144_FITNESS_STATS a, alexis.anb_MET_OE2206_FS1_PMEDTA_nobag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and b.hours = 24) h
on a.strain_id = h.strain_id and a.orf_name = h.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS1_PMEDTA_inbag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS1_PMEDTA_inbag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS1_PMEDTA_inbag_score
from alexis.anb_MET_OE2206_FS1_PMEDTA_inbag_6144_FITNESS_STATS a, alexis.anb_MET_OE2206_FS1_PMEDTA_inbag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 35) i
on a.strain_id = i.strain_id and a.orf_name = i.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS1_MMEDTA_nobag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS1_MMEDTA_nobag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS1_MMEDTA_nobag_score
from aaron_z.aaz_MET_OE2206_FS1_MMEDTA_nobag_6144_FITNESS_STATS a, aaron_z.aaz_MET_OE2206_FS1_MMEDTA_nobag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 161) j
on a.strain_id = j.strain_id and a.orf_name = j.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS1_MMEDTA_inbag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS1_MMEDTA_inbag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS1_MMEDTA_inbag_score
from Brandon.bg_MET_OE2206_FS1_MMEDTA_inbag_6144_FITNESS_STATS a, Brandon.bg_MET_OE2206_FS1_MMEDTA_inbag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 117) k
on a.strain_id = k.strain_id and a.orf_name = k.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS3_GA_nobag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS3_GA_nobag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS3_GA_nobag_score
from aaron_z.aaz_MET_OE2206_FS3_GA_nobag_6144_FITNESS_STATS a, aaron_z.aaz_MET_OE2206_FS3_GA_nobag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 34) l
on a.strain_id = l.strain_id and a.orf_name = l.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS3_GA_inbag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS3_GA_inbag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS3_GA_inbag_score
from alexis.anb_MET_OE2206_FS3_GA_inbag_6144_FITNESS_STATS a, alexis.anb_MET_OE2206_FS3_GA_inbag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 24) m
on a.strain_id = m.strain_id and a.orf_name = m.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS3_PM_nobag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS3_PM_nobag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS3_PM_nobag_score
from aaron_z.aaz_MET_OE2206_FS3_PM_nobag_6144_FITNESS_STATS a, aaron_z.aaz_MET_OE2206_FS3_PM_nobag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 26) n
on a.strain_id = n.strain_id and a.orf_name = n.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS3_PM_inbag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS3_PM_inbag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS3_PM_inbag_score
from alexis.anb_MET_OE2206_FS3_PM_inbag_6144_FITNESS_STATS a, alexis.anb_MET_OE2206_FS3_PM_inbag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 42) o
on a.strain_id = o.strain_id and a.orf_name = o.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS3_MM_inbag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS3_MM_inbag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS3_MM_inbag_score
from Brandon.bg_MET_OE2206_FS3_MM_inbag_6144_FITNESS_STATS a, Brandon.bg_MET_OE2206_FS3_MM_inbag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 163) p
on a.strain_id = p.strain_id and a.orf_name = p.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS3_PMEDTA_nobag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS3_PMEDTA_nobag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS3_PMEDTA_nobag_score
from Brandon.bg_MET_OE2206_FS3_PMEDTA_nobag_6144_FITNESS_STATS a, Brandon.bg_MET_OE2206_FS3_PMEDTA_nobag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 57) q
on a.strain_id = q.strain_id and a.orf_name = q.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS3_PMEDTA_inbag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS3_PMEDTA_inbag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS3_PMEDTA_inbag_score
from Brandon.bg_MET_OE2206_FS3_PMEDTA_inbag_6144_FITNESS_STATS a, Brandon.bg_MET_OE2206_FS3_PMEDTA_inbag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 154) r
on a.strain_id = r.strain_id and a.orf_name = r.orf_name
left join

(select a.strain_id, a.orf_name, a.cs_median FS3_MMEDTA_inbag_fitness,
(case when b.p > 0.05 then 'neutral' when b.p <= 0.05 and b.stat > 0 then 'beneficial' else 'deleterious' end) as FS3_MMEDTA_inbag_phenotype,
(case when b.p <= 0.05 and b.stat > 0 then 1 else 0 end) as FS3_MMEDTA_inbag_score
from alexis.anb_MET_OE2206_FS3_MMEDTA_inbag_6144_FITNESS_STATS a, alexis.anb_MET_OE2206_FS3_MMEDTA_inbag_6144_PVALUE b
where a.strain_id = b.strain_id and a.hours = b.hours and a.hours = 163) s
on a.strain_id = s.strain_id and a.orf_name = s.orf_name
) matrix"

fit.matrix <- dbGetQuery(conn, q)
