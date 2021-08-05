library(pzfx)
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(ggridges)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(plotly)
library(scales)
library(reshape2)
library(rstatix)
library(gtools)
library(effsize)
library(RMariaDB)

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/methionine/functions/colorstrip.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

fig_path <- "~/R/Projects/methionine/figures"
expt.name <- "carbon_source"

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

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
crbn_src_data <- NULL
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
      
      crbn_src_data <- rbind(crbn_src_data, temp)
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
    
    crbn_src_data <- rbind(crbn_src_data, temp)
  }
}
crbn_src_data$expt_rep[crbn_src_data$expt_rep == 'R1'] <- 'rep1'
crbn_src_data$expt_rep[crbn_src_data$expt_rep == 'R2'] <- 'rep2'
crbn_src_data$expt_rep[crbn_src_data$expt_rep == 'R3'] <- 'rep3'
crbn_src_data$orf_name[crbn_src_data$orf_name == 'FY4_met15del'] <- 'FY4-met15del'
crbn_src_data$orf_name[crbn_src_data$orf_name == 'FY4_met3del'] <- 'FY4-met3del'
crbn_src_data$bio_rep[crbn_src_data$row%%2==1 & crbn_src_data$col%%2==1] = '1'
crbn_src_data$bio_rep[crbn_src_data$row%%2==0 & crbn_src_data$col%%2==1] = '3'
crbn_src_data$bio_rep[crbn_src_data$row%%2==1 & crbn_src_data$col%%2==0] = '2'
crbn_src_data$bio_rep[crbn_src_data$row%%2==0 & crbn_src_data$col%%2==0] = '4'

crbn_src_data$condition <- factor(crbn_src_data$condition, levels = c('MiM_Glu','MiM_Gal','MiM_Et',
                                                                      'MiU_Glu','MiU_Gal','MiU_Et',
                                                                      'GLU','GAL','PlM_Glu'))
strain.labs <- c('FY4','FY4-*met3Δ*','FY4-*met15Δ*','BY4742','BY4741')
crbn_src_data <- crbn_src_data[crbn_src_data$expt_rep == 'rep1',c('condition','expt_rep','bio_rep','orf_name','relative_fitness')] 

anova.res <- NULL
for (c in unique(crbn_src_data$condition)) {
  for (s in unique(crbn_src_data$orf_name)) {
    if (s != 'FY4') {
      res.rcs.aov <- crbn_src_data[crbn_src_data$condition == c & crbn_src_data$orf_name %in% c(s, 'FY4'),] %>%
        data.frame() %>%
        anova_test(relative_fitness ~ orf_name * bio_rep)

      cs1 <- crbn_src_data$relative_fitness[crbn_src_data$condition == c & crbn_src_data$orf_name == s]
      cs2 <- crbn_src_data$relative_fitness[crbn_src_data$condition == c & crbn_src_data$orf_name == 'FY4']
      emp_effsize <- (median(cs1, na.rm = T) - median(cs2, na.rm = T))/median(cs2, na.rm = T)
      
      anova.res <- rbind(anova.res, cbind(c, s,
                                          res.rcs.aov$p[res.rcs.aov$Effect == 'orf_name'],
                                          res.rcs.aov$p[res.rcs.aov$Effect == 'orf_name:bio_rep'],
                                          emp_effsize))
    }
  }
}
colnames(anova.res) <- c('condition','strain','rcs_between','rcs_within','effect_size')
anova.res <- data.frame(anova.res, stringsAsFactors = F)
anova.res$rcs_between <- as.numeric(anova.res$rcs_between)
anova.res$rcs_within <- as.numeric(anova.res$rcs_within)
anova.res$effect_size <- as.numeric(anova.res$effect_size)
anova.res$rcs_between <- p.adjust(anova.res$rcs_between, method = 'BH')
anova.res$rcs_within <- p.adjust(anova.res$rcs_within, method = 'BH')
anova.res$label[anova.res$rcs_between > 0.05] <- 'ns'
anova.res$label[anova.res$rcs_between <= 0.05] <- '*'
anova.res$label[anova.res$rcs_between <= 0.01] <- '**'
anova.res$label[anova.res$rcs_between <= 0.001] <- '***'
anova.res$label[anova.res$rcs_between <= 0.0001] <- '****'
 
anova.res$condition <- factor(anova.res$condition, levels = c('MiM_Glu','MiM_Gal','MiM_Et',
                                                              'MiU_Glu','MiU_Gal','MiU_Et',
                                                              'GLU','GAL','PlM_Glu'))
# anova.res$condition <- factor(anova.res$condition, levels = c('MiM_Glu','MiU_Glu','MiM_Gal','MiU_Gal','MiM_Et','MiU_Et','PlM_Glu'))

plot.rcs.box.aov <- crbn_src_data[crbn_src_data$orf_name != "FY4-met3del",] %>%
  # summarise(relative_fitness = median(relative_fitness, na.rm = T)) %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  geom_boxplot(aes(fill = as.character(bio_rep)), size = .3, outlier.shape = NA) +
  geom_text(data = anova.res[anova.res$strain != "FY4-met3del",],
            aes(x = strain, y = 1.5, label = label), size = 2.2, col = 'red') +
  geom_text(data = anova.res[anova.res$strain != "FY4-met3del",],
            aes(x = strain, y = 1.4, label = sprintf('%0.2f%%',effect_size*100)),
            size = 2.2) +
  labs(x = 'Strains', y = 'Relative Colony Size') +
  scale_fill_manual(name = 'Biological Replicate',
                    values = c("1" = "#673AB7",
                               "2" = "#009688",
                               "3" = "#607D8B",
                               "4" = "#FFC107")) +
  scale_x_discrete(limits = c('FY4','FY4-met3del','FY4-met15del','BY4742','BY4741')[-2],
                   labels = strain.labs[-2]) +
  facet_wrap(.~condition, nrow = 3, 
             labeller = labeller(condition = 
                                   c('MiM_Glu'  = 'SD-Met-Cys+Glu',
                                     'MiM_Gal'  = 'SD-Met-Cys+Gal',
                                     'MiM_Et'   = 'SD-Met-Cys+EtOH',
                                     'PlM_Glu'  = 'SD+Met-Cys+Glu',
                                     'MiU_Glu'  = 'SD-Ura+Glu',
                                     'MiU_Gal'  = 'SD-Ura+Gal',
                                     'MiU_Et'   = 'SD-Ura+EtOH',
                                     'GLU'      = 'SC-Met-Cys+Glu',
                                     'GAL'      = 'SC-Met-Cys+Gal'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))
plot.rcs.box.aov <- colorstrip(plot.rcs.box.aov, c('#7C4DFF','#9C27B0','#303F9F',
                                                   '#D32F2F','#FF5722','#757575',
                                                   '#448AFF','#00BCD4','#8BC34A'))
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX_ANOVA.jpg",fig_path, expt.name),
       plot.rcs.box.aov,
       height = one.5c, width = one.5c, units = 'mm',
       dpi = 600)

##### REMOVE MET3 and MAKE SUPPLEMENTAL FIGURES WITH FY4 and MET3 ONLY
plot.rcs.box.aov <- crbn_src_data[crbn_src_data$orf_name %in% c("FY4", "FY4-met3del"),] %>%
  # summarise(relative_fitness = median(relative_fitness, na.rm = T)) %>%
  ggplot(aes(x = orf_name, y = relative_fitness)) +
  geom_boxplot(aes(fill = as.character(bio_rep)), size = 0.3, outlier.shape = NA) +
  geom_text(data = anova.res[anova.res$strain %in% c("FY4", "FY4-met3del"),],
            aes(x = strain, y = 1.5, label = label), size = 2.2, col = 'red') +
  geom_text(data = anova.res[anova.res$strain %in% c("FY4", "FY4-met3del"),],
            aes(x = strain, y = 1.4, label = sprintf('%0.2f%%',effect_size*100)),
            size = 2.2) +
  labs(x = 'Strains', y = 'Relative Colony Size') +
  scale_fill_manual(name = 'Biological Replicate',
                    values = c("1" = "#673AB7",
                               "2" = "#009688",
                               "3" = "#607D8B",
                               "4" = "#FFC107")) +
  scale_x_discrete(limits = c('BY4742','BY4741','FY4','FY4-met3del','FY4-met15del')[c(3,4)],
                   labels = strain.labs[c(1,2)]) +
  facet_wrap(.~condition, nrow = 3, 
             labeller = labeller(condition = 
                                   c('MiM_Glu' = 'SD-Met-Cys+Glu',
                                     'MiM_Gal' = 'SD-Met-Cys+Gal',
                                     'MiM_Et' = 'SD-Met-Cys+EtOH',
                                     'PlM_Glu' = 'SD+Met-Cys+Glu',
                                     'MiU_Glu' = 'SD-Ura+Glu',
                                     'MiU_Gal' = 'SD-Ura+Gal',
                                     'MiU_Et' = 'SD-Ura+EtOH',
                                     'GLU'      = 'SC-Met-Cys+Glu',
                                     'GAL'      = 'SC-Met-Cys+Gal'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, hjust = 1),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(0,1.5))
plot.rcs.box.aov <- colorstrip(plot.rcs.box.aov, c('#7C4DFF','#9C27B0','#303F9F',
                                                   '#D32F2F','#FF5722','#757575',
                                                   '#448AFF','#00BCD4','#8BC34A'))
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX_ANOVA2.jpg",fig_path, expt.name),
       plot.rcs.box.aov,
       height = one.5c, width = one.5c, units = 'mm',
       dpi = 600)

