hi <- dbGetQuery(conn, 'select * from brian_031918.`FITNESS_v2_KO1-MET_SD_screen` where exp_id = 20')
hi %>%
  ggplot(aes(x = average)) +
  geom_line(stat = 'density', aes(col = orf_type))

hi$orf_type[hi$orf_name %in% unique(hi2$orf_name)] <- 'SAP'
hi$orf_type[is.na(hi$orf_type)] <- 'not SAP'



hi2 <- dbGetQuery(conn, "select * from brian_031918.`FITNESS_v2_KO1-MET_SD_screen`
where orf_name in ('YAL012W', 'YER043C', 'YER091C', 'YKL001C', 'YFR030W', 'YDR502C', 'YGL125W', 'YGL184C', 'YGR155W', 'YJR010W', 'YJR130C', 'YJR137C', 'YLR180W', 'YLR303W', 'YNL277W', 'YPR167C', 'YLL058W') and exp_id = 20")

hi2 %>%
  ggplot(aes(x = average)) +
  geom_line(stat = 'density')


median(hi$average, na.rm = T)

median(hi2$average, na.rm = T)

"select * from `FITNESS_v2_KO1-MET_SD_screen` a#, KO_pos2coor_new2 b
where a.orf_name in ('YAL012W', 'YER043C', 'YER091C', 'YKL001C', 'YFR030W', 'YDR502C', 'YGL125W', 'YGL184C', 'YGR155W', 'YJR010W', 'YJR130C', 'YJR137C', 'YLR180W', 'YLR303W', 'YNL277W', 'YPR167C', 'YLL058W') and a.exp_id = 20 
order by orf_name #and a.pid = b.position"