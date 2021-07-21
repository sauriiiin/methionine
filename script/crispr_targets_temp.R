
data <- read.delim(file = '/home/sbp29/flashfry/scer/allorfs.output.scored', sep = '\t', header = T, stringsAsFactors = F)
head(data)

temp <- str_split(data$contig, pattern = ',', simplify = T)
temp <- str_split(temp[,1], pattern = '_', simplify = T)
temp <- temp[,c(1,2)]
colnames(temp) <- c('orf_name','standard_name')

data <- cbind(temp, data)
data <- data[,-3]
head(data)

data$candidates[data$Hsu2013 >= 95 & data$Doench2014OnTarget >= 0.9 & !is.na(data$Doench2014OnTarget)] <- 'Good' 
data$candidates[is.na(data$candidates)] <- 'Not Good'

data %>%
  filter(candidates == 'Good') %>%
  group_by(orf_name) %>% 
  data.frame() %>%
  plyr::count() 

library(seqinr)

smorfs <- dbGetQuery(conn, "select orf_name, seq_nt_atg from ORFS_SEQUENCES
                     where orf_name in
                     (select orf_name from PROTOGENES
                     where orf_name like 'sm%' and pg_2012 = 1)")
write.fasta(as.list(smorfs$seq_nt_atg), as.list(smorfs$orf_name), file.out = 'data/flashfry/smorfs.fasta', open = "w")


data %>%
  ggplot(aes(x = Doench2014OnTarget, y = Hsu2013)) +
  geom_point(size = 0.2) +
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
                                  margin = margin(0.1,0,0.1,0, "mm")))

data %>%
  ggplot() +
  geom_line(aes(y = Hsu2013), stat = 'density') +
  plot.theme + theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

ggpubr::ggarrange(plot.off, NULL, plot.onandoff, plot.on, ncol = 2, nrow = 2,  align = "hv")
