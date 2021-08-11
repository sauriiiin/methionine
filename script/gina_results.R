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
expt.name <- "gina"

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9

pzfx_tables("/home/sbp29/R/Projects/methionine/data/gina/072321_gina_data.pzfx")
data.gina <- read_pzfx("/home/sbp29/R/Projects/methionine/data/gina/072321_gina_data.pzfx", table = "Data 2")
head(data.gina)

#####
data.gina <- melt(data.gina, id.vars = 'Time (min)', variable.name = 'ID', value.name = 'uM')
data.gina <- cbind(data.gina, str_split(data.gina$ID, '_', simplify = T))
colnames(data.gina) <- c(colnames(data.gina)[1:3], 'Sample', 'Replicate')
data.gina$Sample <- factor(data.gina$Sample, levels = c('Met15', 'Yll058w', 'None'))

data.gina$outlier <- NULL
for (t in unique(data.gina$`Time (min)`)) {
  for (s in unique(data.gina$Sample[data.gina$`Time (min)` == t])) {
    data.gina$outlier[data.gina$`Time (min)` == t & data.gina$Sample == s] <-
      isoutlier(data.gina$uM[data.gina$`Time (min)` == t & data.gina$Sample == s],3)
  }
}

#####
ttest.res <- NULL
for (s1 in unique(data.gina$Sample)) {
  for (s2 in unique(data.gina$Sample[data.gina$Sample != s1])) {
    temp.t <- t.test(uM ~ Sample, data = data.gina %>%
                       filter(Sample %in% c(s1, s2)) %>%
                       group_by(`Time (min)`, Sample) %>%
                       summarise(.groups = 'keep', uM = mean(uM, na.rm = T)) %>%
                       data.frame(), paired = TRUE, alternative = 'greater')
    ttest.res <- rbind(ttest.res, data.frame(reference = s1, query = s2, p = temp.t$p.value, stat = temp.t$statistic[[1]]))
  }
}
head(ttest.res)

ttest.res$label[ttest.res$p > 0.05] <- 'ns'
ttest.res$label[ttest.res$p <= 0.05] <- '*'
ttest.res$label[ttest.res$p <= 0.01] <- '**'
ttest.res$label[ttest.res$p <= 0.001] <- '***'
ttest.res$label[ttest.res$p <= 0.0001] <- '****'

ttest.res <- merge(ttest.res, data.gina %>%
                     filter(`Time (min)` == 15) %>%
                     group_by(Sample) %>%
                     summarise(.groups = 'keep', uM = mean(uM, na.rm = T)) %>%
                     data.frame(), by.x = 'query', by.y = 'Sample')


#####
for (t in unique(data.gina$`Time (min)`)) {
  for (s in unique(data.gina$Sample[data.gina$`Time (min)` == t & data.gina$Sample != 'None'])) {
    temp.kw <- compare_means(data = data.gina[data.gina$`Time (min)` == t & data.gina$Sample %in% c(s,'None') &
                                                data.gina$outlier == FALSE,],
                                method = 't.test', formula = uM ~ Sample) %>% data.frame()
    
  }
}

#####
plot.gina.res <- data.gina %>%
  # filter(outlier == FALSE) %>%
  ggplot(aes(x = `Time (min)`, y = uM)) +
  stat_summary(aes(group = Sample, col = Sample), fun=mean, geom="line", lwd = 0.7) +
  stat_summary(aes(group = Sample, col = Sample), fun.data = mean_se, geom = "errorbar", lwd = 0.7) +
  stat_summary(aes(group = Sample), fun=mean, geom="point", size =2) +
  stat_summary(aes(group = Sample, col = Sample), fun=mean, geom="point", size = 0.5) +
  labs(y = 'Homocysteine (log2(μM))') +
  scale_color_manual(name = 'Sample',
                     values = c('Met15' = "#607D8B",
                                'Yll058w' = "#9E9E9E",
                                'None' = "#212121")) +
  scale_y_continuous(trans = 'log2') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm"),
        legend.box.spacing = unit(0.5,"mm"),
        strip.text = element_text(size = txt,
                                  face = 'bold',
                                  margin = margin(0.1,0,0.1,0, "mm"))) +
  coord_cartesian(ylim = c(5,125))
fig6C <- plot.gina.res
save(fig6C, file = 'figures/final/fig6C.RData')
# ggsave(sprintf("%s/%s/GINA_RESULTS.jpg",fig_path, expt.name), plot.gina.res,
#        height = one.c, width = one.c, units = 'mm',
#        dpi = 600)

#####







