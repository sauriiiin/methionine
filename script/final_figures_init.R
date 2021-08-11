

library(ggplot2)
library(ggridges)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(scales)
library(egg)
library(zoo)
library(png)
library(ggrepel)
library(reshape2)

load("~/R/Projects/methionine/figures/final/fig6C.RData")
load("~/R/Projects/methionine/figures/final/fig6B.RData")
load("~/R/Projects/methionine/figures/final/fig5.RData")
load("~/R/Projects/methionine/figures/final/fig3B.RData")
load("~/R/Projects/methionine/figures/final/fig3A.RData")
load("~/R/Projects/methionine/figures/final/fig1C.RData")
load("~/R/Projects/methionine/figures/final/fig1B.RData")

fig_path <- "~/R/Projects/methionine/figures"

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 7
lbls <- 9


fig1A <- readPNG('figures/branden/Fig1A.png')
fig1A <- ggplot() + 
  background_image(fig1A) +
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
        plot.background = element_blank())

# fig3B <- readPNG('figures/branden/Fig3B.png')
# fig3B <- ggplot() + 
#   background_image(fig3B) +
#   theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
#         plot.background = element_blank())
# 
# fig3C <- readPNG('figures/branden/Fig3C.png')
# fig3C <- ggplot() + 
#   background_image(fig3C) +
#   theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
#         plot.background = element_blank())

fig6A <- rectGrob(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                  width = unit(1, "npc"), height = unit(1, "npc"),
                  just = "centre", hjust = NULL, vjust = NULL,
                  default.units = "npc", name = NULL,
                  gp=gpar(col=NA), vp = NULL)

fig6D <- readPNG('figures/branden/Fig6D.png')
fig6D <- ggplot() + 
  background_image(fig6D) +
  theme(plot.margin = margin(t=0, l=0, r=0, b=0, unit = "mm"),
        plot.background = element_blank())



fig1 <- cowplot::plot_grid(fig1A, cowplot::plot_grid(fig1B, fig1C, ncol = 2,
                                                     align = 'v', axis = 'l', rel_widths = c(1.5,1),
                                                     labels = c('B','C'),
                                                     label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                           ncol = 1, align = 'v', axis = 'l', rel_heights = c(1,1.5),
                           labels = c('A',''),
                           label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/final/Figure1.jpg",fig_path), fig1,
       height = 250, width = two.c, units = 'mm',
       dpi = 600)


fig3 <- cowplot::plot_grid(fig3A, fig3B,
                           labels = c('A','B'), ncol = 2, rel_widths = c(2,1),
                           align = 'hv', axis = 'tb',
                           label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/final/Figure3.jpg",fig_path), fig3,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)


ggsave(sprintf("%s/final/Figure5.jpg",fig_path), fig5,
       height = 140, width = two.c, units = 'mm',
       dpi = 600)


fig6 <- cowplot::plot_grid(fig6A, fig6B, fig6C, fig6D,
                           ncol = 2, align = 'v', axis = 'l',
                           labels = c('A','B','C','D'),
                           label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/final/Figure6.jpg",fig_path), fig6,
       height = one.5c, width = two.c, units = 'mm',
       dpi = 600)
