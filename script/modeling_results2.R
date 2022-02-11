
# figS9 <- data.mdl %>%
#   filter(label %in% c('C. B + Hypothesized YLL058W reaction',
#                       'D. C - All MET15 reactions')) %>%
#   ggplot(aes(x = log(Yll058w.Met15,10), y = Growth)) +
#   geom_line(aes(col = label), lwd = 1) +
#   geom_point(size = 2) +
#   geom_point(aes(col = label), size = 0.5) +
#   scale_x_continuous(trans = 'reverse',
#                      breaks = seq(10,-10,-1),
#                      labels = 10^seq(10,-10,-1)) +
#   scale_color_manual(name = 'Model with',
#                      values = c('C. B + Hypothesized YLL058W reaction' = '#333333',
#                                 'D. C - All MET15 reactions' = '#999999'),
#                      labels = c('Hypothesized *YLL058W* reaction<br/>w/ all *MET15* reactions',
#                                 'Hypothesized *YLL058W* reaction<br/>w/o all *MET15* reactions')) +
#   labs(x = 'Yll058w/Met15 Kcat Ratio',
#        y = 'Simulated Biomass Flux\n("Growth")') +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = ggtext::element_markdown(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm"))) +
#   guides(color = guide_legend(nrow = 2, byrow=F, order = 1,
#                               override.aes = list(shape = 15, size = 3, alpha = 1)))
# 
# ggsave(sprintf("%s/FigureS9.jpg",fig_path), figS9,
#        height = one.c, width = one.c, units = 'mm',
#        dpi = 300)


# figS9 <- melt(data.mdl.2, id.vars = c('Yll058w_Met15','Model'),
#               variable.name = 'Reaction', value.name = 'Flux') %>%
#   ggplot(aes(x = log(Yll058w_Met15,10), y = Flux)) +
#   geom_line(aes(col = Model), lwd = 1) +
#   geom_point(size = 2) +
#   geom_point(aes(col = Model), size = 0.5) +
#   scale_x_continuous(trans = 'reverse',
#                      breaks = seq(10,-10,-1),
#                      labels = 10^seq(10,-10,-1)) +
#   scale_color_manual(name = 'Model with',
#                      values = c('C' = '#333333',
#                                 'D' = '#999999'),
#                      labels = c('C'='Hypothesized *YLL058W* reaction w/ all *MET15* reactions',
#                                 'D'='Hypothesized *YLL058W* reaction w/o all *MET15* reactions')) +
#   labs(x = 'Yll058w/Met15 Kcat Ratio',
#        y = 'Simulated Flux') +
#   facet_wrap(.~Reaction, scales = 'free_y',
#              labeller = as_labeller(c('Biomass_Flux' = 'Overall biomass flux\n("Growth")',
#                                       'HS_Flux' = 'Flux through homocysteine\nproducing reactions'))) +
#   theme_linedraw() +
#   theme(plot.title = element_text(size = titles + 2, face = 'bold', hjust = 0.5),
#         axis.title = element_text(size = titles),
#         axis.text = element_text(size = txt),
#         legend.title = element_text(size = titles),
#         legend.text = ggtext::element_markdown(size = txt),
#         legend.position = 'bottom',
#         legend.key.size = unit(3, "mm"),
#         legend.box.spacing = unit(0.5,"mm"),
#         strip.text = element_text(size = txt,
#                                   face = 'bold',
#                                   margin = margin(0.1,0,0.1,0, "mm"))) +
#   guides(color = guide_legend(nrow = 2, byrow=F, order = 1,
#                               override.aes = list(shape = 15, size = 3, alpha = 1)))
# ggsave(sprintf("%s/FigureS9.jpg",fig_path), figS9,
#        height = one.c, width = two.c, units = 'mm',
#        dpi = 300)