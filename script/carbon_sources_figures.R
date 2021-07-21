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

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")

fig_path <- "~/R/Projects/methionine/figures"
expt.name <- "carbon_source"

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9

crbn_src_data <- NULL
for (t in pzfx_tables("/home/sbp29/R/Projects/methionine/data/carbon_source/box_plots.pzfx")) {
  temp <- read_pzfx("/home/sbp29/R/Projects/methionine/data/carbon_source/box_plots.pzfx", table = t)
  for(i in 1:dim(temp)[2]){
    temp[,i][isoutlier(temp[,i],2)] <- NA
  }
  temp <- temp/median(temp$FY4, na.rm = T)
  
  temp$condition <- t
  crbn_src_data <- rbind(crbn_src_data, temp)
}
crbn_src_data <- cbind(crbn_src_data, bio_rep = c(1,2,3,4))
crbn_src_data$condition <- factor(crbn_src_data$condition, levels = c('MiM_Glu','MiM_Gal','MiM_Et','PlM_Glu',
                                                                      'MiU_Glu','MiU_Gal','MiU_Et'))
strain.labs <- c('FY4','*met3Δ*', '*met15Δ*', 'BY4742', 'BY4741')

crbn_src_data <- crbn_src_data %>% 
  melt(id.vars = c('condition', 'bio_rep'), variable.name = 'strain', value.name = "relative_fitness") %>%
  group_by(condition, strain, bio_rep)

anova.res <- NULL
for (c in unique(crbn_src_data$condition)) {
  for (s in unique(crbn_src_data$strain)) {
    if (s != 'FY4') {
      res.rcs.aov <- crbn_src_data[crbn_src_data$condition == c & crbn_src_data$strain %in% c(s, 'FY4'),] %>%
        data.frame() %>%
        anova_test(relative_fitness ~ strain * bio_rep)

      cs1 <- crbn_src_data$relative_fitness[crbn_src_data$condition == c & crbn_src_data$strain == s]
      cs2 <- crbn_src_data$relative_fitness[crbn_src_data$condition == c & crbn_src_data$strain == 'FY4']
      emp_effsize <- (median(cs1, na.rm = T) - median(cs2, na.rm = T))/median(cs2, na.rm = T)
      
      anova.res <- rbind(anova.res, cbind(c, s,
                                          res.rcs.aov$p[res.rcs.aov$Effect == 'strain'],
                                          res.rcs.aov$p[res.rcs.aov$Effect == 'strain:bio_rep'],
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
 
anova.res$condition <- factor(anova.res$condition, levels = c('MiM_Glu','MiM_Gal','MiM_Et','PlM_Glu','MiU_Glu','MiU_Gal','MiU_Et'))

plot.rcs.box.aov <- crbn_src_data[crbn_src_data$strain != "met3D",] %>%
  # summarise(relative_fitness = median(relative_fitness, na.rm = T)) %>%
  ggplot(aes(x = strain, y = relative_fitness)) +
  geom_boxplot(aes(fill = as.character(bio_rep)), outlier.shape = NA) +
  geom_text(data = anova.res[anova.res$strain != "met3D",],
            aes(x = strain, y = 1.5, label = label), size = 2, col = 'red') +
  geom_text(data = anova.res[anova.res$strain != "met3D",],
            aes(x = strain, y = 1.4, label = sprintf('%0.2f%%',effect_size*100)),
            size = 1.5) +
  labs(x = 'Strains', y = 'Relative Colony Size') +
  scale_fill_discrete(name = 'Biological Replicate') +
  scale_x_discrete(limits = c('FY4','met3D','met15D','BY4742','BY4741')[-2],
                   labels = strain.labs[-2]) +
  # geom_jitter() +
  facet_wrap(.~condition, nrow = 2, 
             labeller = labeller(condition = 
                                   c('MiM_Glu' = 'SD-Met-Cys+Glu',
                                     'MiM_Gal' = 'SD-Met-Cys+Gal',
                                     'MiM_Et' = 'SD-Met-Cys+EtOH',
                                     'PlM_Glu' = 'SD+Met-Cys+Glu',
                                     'MiU_Glu' = 'SD-Ura+Glu',
                                     'MiU_Gal' = 'SD-Ura+Gal',
                                     'MiU_Et' = 'SD-Ura+EtOH'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.x = ggtext::element_markdown(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX_ANOVA.jpg",fig_path, expt.name),
       plot.rcs.box.aov,
       height = 140, width = two.c, units = 'mm',
       dpi = 600)

##### REMOVE MET3 and MAKE SUPPLEMENTAL FIGURES WITH FY4 and MET3 ONLY
plot.rcs.box.aov <- crbn_src_data[crbn_src_data$strain %in% c("FY4", "met3D"),] %>%
  # summarise(relative_fitness = median(relative_fitness, na.rm = T)) %>%
  ggplot(aes(x = strain, y = relative_fitness)) +
  geom_boxplot(aes(fill = as.character(bio_rep)), outlier.shape = NA) +
  geom_text(data = anova.res[anova.res$strain %in% c("FY4", "met3D"),],
            aes(x = strain, y = 1.5, label = label), size = 2, col = 'red') +
  geom_text(data = anova.res[anova.res$strain %in% c("FY4", "met3D"),],
            aes(x = strain, y = 1.4, label = sprintf('%0.2f%%',effect_size*100)),
            size = 1.5) +
  labs(x = 'Strains', y = 'Relative Colony Size') +
  scale_fill_discrete(name = 'Biological Replicate') +
  scale_x_discrete(limits = c('FY4','met3D','met15D','BY4742','BY4741')[1:2],
                   labels = strain.labs[1:2]) +
  # geom_jitter() +
  facet_wrap(.~condition, nrow = 1, 
             labeller = labeller(condition = 
                                   c('MiM_Glu' = 'SD-Met-Cys+Glu',
                                     'MiM_Gal' = 'SD-Met-Cys+Gal',
                                     'MiM_Et' = 'SD-Met-Cys+EtOH',
                                     'PlM_Glu' = 'SD+Met-Cys+Glu',
                                     'MiU_Glu' = 'SD-Ura+Glu',
                                     'MiU_Gal' = 'SD-Ura+Gal',
                                     'MiU_Et' = 'SD-Ura+GtOH'))) +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.x = ggtext::element_markdown(),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm")))
ggsave(sprintf("%s/%s/RELATIVE_COLONY_SIZE_BOX_ANOVA2.jpg",fig_path, expt.name),
       plot.rcs.box.aov,
       height = 70, width = two.c, units = 'mm',
       dpi = 600)

