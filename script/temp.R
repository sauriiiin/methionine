
data.jm.2 <- dbGetQuery(conn, 'select * from
                        (select a.pos, a.hours, a.average, c.density, c.plate, c.row, c.col, b.orf_name
                        from Branden.jakes_mutants_c1_MiM_Glu_new_384_CLEAN a, Branden.jakes_mutants_pos2orf_name b, Branden.jakes_mutants_pos2coor c
                        where a.pos = b.pos and b.pos = c.pos
                        union
                        select a.pos, a.hours, a.average, c.density, c.plate, c.row, c.col, b.orf_name
                        from Branden.jakes_mutants_c2_MiM_Glu_new_384_CLEAN a, Branden.jakes_mutants_pos2orf_name b, Branden.jakes_mutants_pos2coor c
                        where a.pos = b.pos and b.pos = c.pos) d
                        order by d.hours, d.plate, d.col, d.row')

data.jm.2 %>%
  # filter(orf_name == 'BY4741') %>%
  ggplot(aes(x = hours, y = average)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               aes(group = orf_name, fill = orf_name), geom="ribbon", alpha = 0.4) +
  stat_summary(aes(group = orf_name, col = orf_name), fun=mean, geom="line", lwd =0.7)



data <- dbGetQuery(conn, 'select * from
                        (select a.pos, a.hours, a.average, c.density, c.plate, c.row, c.col, b.orf_name
                        from Branden.carbon_rep1_MiM_Glu_384_CLEAN a, Branden.carbon_rep1_pos2orf_name b, Branden.carbon_rep1_pos2coor c
                        where a.pos = b.pos and b.pos = c.pos
                        union
                        select a.pos, a.hours, a.average, c.density, c.plate, c.row, c.col, b.orf_name
                        from Branden.carbon_rep2_MiM_Glu_384_CLEAN a, Branden.carbon_rep2_pos2orf_name b, Branden.carbon_rep2_pos2coor c
                        where a.pos = b.pos and b.pos = c.pos) d
                        order by d.hours, d.plate, d.col, d.row')

# OR

data <- dbGetQuery(conn, 'select a.pos, a.hours, a.average, c.density, c.plate, c.row, c.col, b.orf_name
                           from Branden.carbon_rep1_MiM_Glu_384_CLEAN a, Branden.carbon_rep1_pos2orf_name b, Branden.carbon_rep1_pos2coor c
                           where a.pos = b.pos and b.pos = c.pos
                           order by a.hours, c.plate, c.col, c.row')

data %>%
  # filter(orf_name != 'BY4742') %>%
  ggplot(aes(x = hours, y = average)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               aes(group = orf_name, fill = orf_name), geom="ribbon", alpha = 0.4) +
  stat_summary(aes(group = orf_name, col = orf_name), fun=mean, geom="line", lwd =0.7)





'select a.pos,
a. nours, a. average,
c. density,
c.plate, C.row,
C.coL, b.orf name
from Branden. respiration_exp1RMiM_Glu_Rep1_BVO_binary_384_CLEAN a, Branden.respiration_exp_1_JX_pos2orf_name b, Branden.respiration_exp_1_JX_pos2coor c
where a.pos = b.pos and b.pos = c.pos
order by a.hours, c.plate, c.col, c.row'

data <- dbGetQuery(conn, 'select a.pos, a.hours, a.average, c.density, c.plate, c.row, c.col, b.orf_name
                           from Branden.respiration_exp1R_MiM_Glu_Rep2_BVO_binary_384_CLEAN a, 
                          Branden.respiration_exp_1_JX_pos2orf_name b, Branden.respiration_exp_1_JX_pos2coor c
                   where a.pos = b.pos and b.pos = c.pos
                   order by a.hours, c.plate, c.col, c.row')

unique(data$hours)
data %>%
  filter(orf_name %in% c('BY4741','BY4742')) %>%
  ggplot(aes(x = hours, y = average)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               aes(group = orf_name, fill = orf_name), geom="ribbon", alpha = 0.4) +
  stat_summary(aes(group = orf_name, col = orf_name), fun=mean, geom="line", lwd =0.7)
