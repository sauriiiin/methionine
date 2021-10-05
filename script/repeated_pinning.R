

library(stringi)
library(dplyr)
library(RMariaDB)

# out_path <- "~/R/Projects/methionine/data"

source("~/R/Projects/adaptivefitness/R/functions/isoutlier.R")
source("~/R/Projects/adaptivefitness/R/functions/initialize.sql.R")
conn <- initialize.sql("saurin_test")

tables <- dbGetQuery(conn, 'SELECT table_name
                     FROM information_schema.tables
                     WHERE table_schema = "Branden"
                     and table_name like "Repeat_%_CLEAN" and table_rows < 10000')
temp <- str_split(tables$table_name, '_', simplify = T) %>% data.frame()

tables <- cbind(tables, temp[,c(2:4)])
colnames(tables) <- c('table_name','aux','carbon','pin')

head(tables)

data.rp <- NULL
for (t in unique(tables$table_name)) {
  temp <- dbGetQuery(conn, sprintf("select a.*, b.density, b.plate, b.row, b.col, c.orf_name
                                   from Branden.%s a, Repeat_11_pos2coor b, Repeat_11_pos2orf_name c
                                   where a.pos = b.pos and b.pos = c.pos
                                   order by a.hours, a.pos", t))
  temp$aux <- tables$aux[tables$table_name == t]
  temp$carbon <- tables$carbon[tables$table_name == t]
  temp$pin <- tables$pin[tables$table_name == t]
  temp$table_name <- tables$table_name[tables$table_name == t]
  
  for (h in unique(temp$hours)) {
    for (o in unique(temp$orf_name)) {
      temp$average[temp$orf_name == o & temp$hours == h][isoutlier(temp$average[temp$orf_name == o & temp$hours == h], 2) |
                                                           isoutlier(temp$fitness[temp$orf_name == o & temp$hours == h], 2)] <- NA
    }
    temp$fitness[temp$hours == h] <- temp$average[temp$hours == h]/
      median(temp$average[temp$hours == h & temp$orf_name == 'FY4'], na.rm = T)
  }
  data.rp <- rbind(data.rp, temp)
}


#####
temp <- data.rp %>%
  group_by(table_name) %>%
  summarise(.groups = 'keep', max_hrs = max(hours, na.rm = T)) %>%
  data.frame()
data.rp <- merge(data.rp, temp, by = 'table_name')

#####
data.rp[str_detect(data.rp$orf_name, 'FY'),] %>%
  filter(hours == max_hrs, orf_name != 'BOR') %>%
  ggplot(aes(x = orf_name, y = fitness)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(aux*carbon~pin)


data.rp %>%
  filter(orf_name != 'BOR') %>%
  ggplot(aes(x = hours, y = average)) +
  # geom_line(aes(col = orf_name)) +
  stat_summary(aes(col = orf_name), geom = 'line', fun = median) +
  facet_grid(aux*carbon~pin)




