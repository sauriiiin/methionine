
library(cowplot)
##### FIGURE 1
fig1 <- plot_grid(fig1a,fig1b,
                  labels = c('A','B'), ncol = 1, rel_heights = c(1,1),
                  label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure1.jpg",fig_path), fig1,
       height = 170, width = two.c, units = 'mm',
       dpi = 600)


##### FIGURE 2
# fig2 <- plot_grid(plot_grid(fig1c, fig2a, fig2b.g,
#                             nrow = 1, rel_widths = c(1,1,1),
#                             labels = c('A','B',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
#                   NULL,
#                   fig2c,
#                   ncol = 1, rel_heights = c(1,0.1,1),
#                   labels = c('','','C'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
# ggsave(sprintf("%s/Figure2.jpg",fig_path), fig2,
#        height = two.c, width = two.c, units = 'mm',
#        dpi = 600)
# 
# with repeated pinning but no petite
fig2 <- plot_grid(plot_grid(fig1c, fig2a, fig2b.g,
                            nrow = 1, rel_widths = c(1,1,1),
                            labels = c('A','B',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  NULL,
                  NULL,
                  fig3a,
                  NULL,
                  ncol = 1, rel_heights = c(1,0.1,0.3,0.8,0.1),
                  labels = c('','','C','',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure2.jpg",fig_path), fig2,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)

##### FIGURE 3
fig3 <- plot_grid(fig3b.g, fig3cp, plot_grid(fig3c, fig3d.g,
                                             nrow = 1, rel_widths = c(1,1), align = 'h', axis = 'tb',
                                             labels = c('C','D'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  NULL,
          ncol = 1, rel_heights = c(1,1,1,0.1),
          labels = c('A','B','',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure3.jpg",fig_path), fig3,
       height = 220, width = two.c, units = 'mm',
       dpi = 600)

##### FIGURE 4
fig4 <- plot_grid(fig4a, plot_grid(fig4c.2, fig4c.3, fig4c,
                           ncol = 1, rel_heights = c(1,1,1),
                           labels = c('B','C','D'), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
          ncol = 2, rel_widths = c(1,2),
          labels = c('A',''), label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')
ggsave(sprintf("%s/Figure4.jpg",fig_path), fig4,
       height = two.c, width = two.c, units = 'mm',
       dpi = 600)

##### FIGURE 5
fig5 <- plot_grid(fig5a,
                  plot_grid(fig5b, fig5c,
                            nrow = 1, rel_widths = c(1.5,1),
                            labels = c('B','C'),
                            label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold'),
                  plot_grid(fig5d, fig5e,
                            ncol = 2, rel_widths = c(1,1),
                            labels = c('D','E'),
                            label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold',
                            align = 'hv', axis = 'tb'),
                  ncol = 1, rel_heights = c(1,1,1),
                  labels = c('A','',''),
                  label_size = lbls, label_fontfamily = 'sans', label_fontface = 'bold')

ggsave(sprintf("%s/Figure5.jpg",fig_path), fig5,
       height = 220, width = two.c, units = 'mm',
       dpi = 600)

##### FIGURE 6



