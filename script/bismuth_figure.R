
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
library(RColorBrewer)

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/methionine/functions/colorstrip.R")

fig_path <- "~/R/Projects/methionine/figures"
expt.name <- "bismuth"

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

strain.levels <- c('FY4', 'FY4-met12D', 'FY4-str3D', 'FY4-met3D', 'FY4-met15D', 'FY4-met2D',
                   'FY4-met6D', 'FY4-met13D', 'FY4-cys4D', 'BY4742', 'BY4741')
strain.labs <- c('FY4', 'FY4-*met12Δ*', 'FY4-*str3Δ*', 'FY4-*met3Δ*', 'FY4-*met15Δ*', 'FY4-*met2Δ*',
                 'FY4-*met6Δ*', 'FY4-*met13Δ*', 'FY4-*cys4Δ*', 'BY4742', 'BY4741')
strain.labs2 <- data.frame(orf_name = strain.levels, labels = strain.labs)
strain.labs2$auxotrophy[strain.labs2$orf_name %in% c('BY4742','FY4','FY4-met12D')] <- 'Prototroph'
strain.labs2$auxotrophy[strain.labs2$orf_name %in% c('BY4741','FY4-met15D','FY4-met3D','FY4-met2D',
                                                     'FY4-met6D','FY4-met13D','FY4-cys4D')] <- 'Presumed Auxotroph'
strain.labs2$auxotrophy[strain.labs2$orf_name %in% c('FY4-str3D')] <- 'Unknown'
strain.labs2$auxotrophy <- factor(strain.labs2$auxotrophy, levels = c('Prototroph', 'Presumed Auxotroph', 'Unknown'))

# color.levels <- colorRampPalette(c("#D7CCC8", "#5D4037"))
# color.levels <- color.levels(7)

##### GATHER DATA
data.bis <- NULL
for (t in pzfx_tables("/home/sbp29/R/Projects/methionine/data/bismuth/bismuth_graphs_2.pzfx")) {
  temp <- read_pzfx("/home/sbp29/R/Projects/methionine/data/bismuth/bismuth_graphs_2.pzfx", table = t)
  temp <- melt(temp, variable.name = 'Strain', value.name = 'HS')
  temp$Condition = t
  data.bis <- rbind(data.bis, temp)
}
data.bis <- data.frame(data.bis)
data.bis$Strain <- as.character(data.bis$Strain)
data.bis$Strain <- factor(data.bis$Strain, levels = strain.levels)
head(data.bis)

##### PLOTTING FIGURE
plot.bar.bis <- data.bis %>%
  # filter(Strain != 'FY4-cys4D') %>%
  ggplot(aes(x = Strain, y = HS)) +
  stat_summary(data = data.bis %>%
                 # filter(Strain != 'FY4-cys4D') %>% 
                 group_by(Condition, Strain) %>%
                 summarize(HS = mean(HS, na.rm = T), .groups = 'keep'),
               aes(fill = round(HS)), col = 'white', alpha = 0.9, size = 1,
               fun = mean, geom = "bar") +
  # stat_summary(fun.data = mean_se, geom = "errorbar") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", fun.args = list(mult = 1)) +
  scale_x_discrete(labels = strain.labs) +
  scale_y_continuous(breaks = seq(0,10,1)) +
  scale_fill_gradient(low = "#D7CCC8", high = "#5D4037", guide = F) +
  labs(y = 'Relative Hydrogen Sulfide',
       x = 'Strains') +
  facet_wrap(.~Condition, ncol = 2) +
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
  coord_cartesian(ylim = c(1,7))
# plot.bar.bis <- colorstrip(plot.bar.bis,c("#212121","#212121"))
fig3B <- plot.bar.bis
save(fig3B, file = 'figures/final/fig3B.RData')



