##### DELETION SCREEN SOM ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 08/01/2021 

##### INITIALIZE
library(ggplot2)
library(gridExtra)
library(grid)
library(tidyverse)
library(ggpubr)
library(stringr)
library(egg)
library(zoo)
library(ggrepel)
library(ggforce)
library(plotly)
library(scales)
library(reshape2)
library(locfit)
library(growthcurver)
library(rstatix)
library(gtools)
library(locfit)
library(growthrates)
library(RMariaDB)
library(genefilter)
library(apeglm)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(kohonen)
library(lattice)

out_path <- "~/R/Projects/methionine/data"
fig_path <- "~/R/Projects/methionine/figures"
res_path <- "~/R/Projects/methionine/results"
expt.name <- "deletion"

source("~/R/Projects/rnaseek/functions/gkEnrichSOM.R")
load(sprintf("%s/%s/210728_environment.RData",out_path,expt.name))

##### SOM FUNCTIONS
# load('~/R/Projects/rnaseek/functions/som_functions.RData')
convert_to_format_for_lattice <- function (  som, input_matrix, id_a, id_b, ncol, nrow) {
  # Check number of empty SOMs clusters
  num_empty_som_clusters <- 0
  for ( i in 1:(ncol*nrow) )  {   
    temp_num_genes  <- length(which(som$unit.classif == i))
    if ( temp_num_genes == 0 ) {  
      num_empty_som_clusters <- num_empty_som_clusters + 1
    }
  }
  data.dimension <- dim( input_matrix)
  # dimension 1 is the number of genes
  # dimension 2 is the number of time points
  data.total_length <- (data.dimension[1] + num_empty_som_clusters) * data.dimension[2]
  data.num_id_a <- data.dimension[1]
  
  som_data_pivot <- data.frame(  id_a= rep(0, data.total_length),  id_b= rep(0, data.total_length), 
                                 x = rep(0, data.total_length), y= rep(0, data.total_length), cluster=rep(0, data.total_length) )
  
  row_count <- 1
  for ( i in 1:(ncol*nrow) ) {    
    temp_input_matrix_data  <- input_matrix[ som$unit.classif == i, ]
    temp_id_a     		    <- id_a[ som$unit.classif == i] 
    temp_id_b     		    <- id_b[ som$unit.classif == i] 
    temp_num_genes  		<- length(which(som$unit.classif == i ) )
    if ( temp_num_genes  > 0 ) {     
      for ( j in 1:temp_num_genes ) {
        my_k_length <- 0
        if ( temp_num_genes > 1 ) {
          my_k_length <- length(temp_input_matrix_data[1,])
        } else {
          my_k_length <- length(temp_input_matrix_data)
        }       
        for ( k in 1:my_k_length ) {
          som_data_pivot[row_count,"id_a"] <- as.character(temp_id_a[j])
          som_data_pivot[row_count,"id_b"] <- as.character(temp_id_b[j])
          som_data_pivot[row_count,"x"]    <- k 
          som_data_pivot[row_count,"cluster"] <- i 
          if ( temp_num_genes > 1 ) {
            som_data_pivot[row_count,"y"] <- temp_input_matrix_data[j, k] 
          } else if( temp_num_genes == 1 ) {
            som_data_pivot[row_count,"y"] <- temp_input_matrix_data[k] 
          }  
          row_count <- row_count + 1
        }
      }
    } else {
      # treat empty clusters                  
      som_data_pivot <-  som_data_pivot 
      for ( num_time_points in 1:data.dimension[2] ) {     
        # id_a, id_b, x, y, cluster, values in that order in the vector
        som_data_pivot[row_count,] <-  c(NA, NA, num_time_points, NA, i )
        row_count <- row_count + 1
      }
    }
  }
  return ( som_data_pivot) 
}

run_soms_analysis <- function ( data_matrix, som_size_x, som_size_y, som_number_seed = 7 ) {
  data_matrix <- as.matrix(data_matrix)
  ### Initialize the X by Y SOM grid, we are using a rectangular grid here
  rectangular_x_by_y <- somgrid(xdim = som_size_x, ydim = som_size_y, topo = c("rectangular"), toroidal = FALSE)
  ### Need to set the seed so that every time you run the code, the results are the same
  set.seed(som_number_seed)
  all.som <- som(data_matrix, grid=rectangular_x_by_y , rlen = 100,   keep.data = TRUE)  
  ### Convert the data into a format that can be used by lattice to plot the SOM
  print("Perform data pivot")
  som_data_pivot <- convert_to_format_for_lattice(all.som, data_matrix,
                                                  rownames(data_matrix), rownames(data_matrix), 
                                                  ncol = som_size_x, nrow = som_size_y) 
  cluster_for_each_id_a <- unique(som_data_pivot[,c("id_a", "id_b", "cluster")])
  return( list(pivot=som_data_pivot, som=all.som) )
}

##### FIGURE SIZE
one.c <- 90 #single column
one.5c <- 140 #1.5 column
two.c <- 190 #full width

##### TEXT SIZE
titles <- 7
txt <- 5
lbls <- 9


##### REFERENCE LIMITS
data.lim2 <- data.lim[,c('arm','stage','hours','fitness_ll','fitness_ul')]
data.lim2$hours[data.lim2$arm == 'SD-Met-Cys+Gal' & data.lim2$stage == 'Pre-Screen #2'] <-
  data.lim2$hours[data.lim2$arm == 'SD-Met-Cys+Gal' & data.lim2$stage == 'Pre-Screen #2'] +
  max(data.lim2$hours[data.lim2$arm == 'SD-Met-Cys+Gal' & data.lim2$stage == 'Pre-Screen #1'])
data.lim2$hours[data.lim2$arm == 'SD-Met-Cys+Gal' & data.lim2$stage == 'Final Screen'] <-
  data.lim2$hours[data.lim2$arm == 'SD-Met-Cys+Gal' & data.lim2$stage == 'Final Screen'] +
  max(data.lim2$hours[data.lim2$arm == 'SD-Met-Cys+Gal' & data.lim2$stage == 'Pre-Screen #2'])
data.lim2$hours[data.lim2$arm == 'SD+Met-Cys+Gal' & data.lim2$stage == 'Pre-Screen #2'] <-
  data.lim2$hours[data.lim2$arm == 'SD+Met-Cys+Gal' & data.lim2$stage == 'Pre-Screen #2'] +
  max(data.lim2$hours[data.lim2$arm == 'SD+Met-Cys+Gal' & data.lim2$stage == 'Pre-Screen #1'])
data.lim2$hours[data.lim2$arm == 'SD+Met-Cys+Gal' & data.lim2$stage == 'Final Screen'] <-
  data.lim2$hours[data.lim2$arm == 'SD+Met-Cys+Gal' & data.lim2$stage == 'Final Screen'] +
  max(data.lim2$hours[data.lim2$arm == 'SD+Met-Cys+Gal' & data.lim2$stage == 'Pre-Screen #2'])

##### FITNESS THROUGH TIME
data.som <- data.sum %>%
  group_by(arm, stage, hours) %>%
  dplyr::count() %>% data.frame()
col.names <- NULL
for (o in unique(data.diff$orf_name)) {
  temp <- data.sum[data.sum$orf_name == o,] %>% group_by(stage, arm, hours) %>% summarize(fitness = mean(fitness, na.rm = T), .groups = 'keep')
  temp <- merge(data.som[,c('stage','arm','hours')], temp[,c('stage','arm','hours','fitness')], by = c('arm','stage','hours'), all = T)
  data.som <- cbind(data.som, temp$fitness)
  col.names <- cbind(col.names,o)
}
colnames(data.som) <- c(colnames(data.som)[1:4], col.names)
data.som$hours[data.som$arm == 'SD-Met-Cys+Gal' & data.som$stage == 'Pre-Screen #2'] <-
  data.som$hours[data.som$arm == 'SD-Met-Cys+Gal' & data.som$stage == 'Pre-Screen #2'] +
  max(data.som$hours[data.som$arm == 'SD-Met-Cys+Gal' & data.som$stage == 'Pre-Screen #1'])
data.som$hours[data.som$arm == 'SD-Met-Cys+Gal' & data.som$stage == 'Final Screen'] <-
  data.som$hours[data.som$arm == 'SD-Met-Cys+Gal' & data.som$stage == 'Final Screen'] +
  max(data.som$hours[data.som$arm == 'SD-Met-Cys+Gal' & data.som$stage == 'Pre-Screen #2'])
data.som$hours[data.som$arm == 'SD+Met-Cys+Gal' & data.som$stage == 'Pre-Screen #2'] <-
  data.som$hours[data.som$arm == 'SD+Met-Cys+Gal' & data.som$stage == 'Pre-Screen #2'] +
  max(data.som$hours[data.som$arm == 'SD+Met-Cys+Gal' & data.som$stage == 'Pre-Screen #1'])
data.som$hours[data.som$arm == 'SD+Met-Cys+Gal' & data.som$stage == 'Final Screen'] <-
  data.som$hours[data.som$arm == 'SD+Met-Cys+Gal' & data.som$stage == 'Final Screen'] +
  max(data.som$hours[data.som$arm == 'SD+Met-Cys+Gal' & data.som$stage == 'Pre-Screen #2'])
data.som[is.na(data.som)] <- 0

data.som2 <- NULL
hh <- NULL
col.names <- NULL
for (a in unique(data.som$arm)) {
  for (o in colnames(data.som)[5:dim(data.som)[[2]]]) {
    temp <- data.som[data.som$arm == a,c('hours',o)]
    temp <- temp[order(temp$hours),]

    lo <- loess.smooth(temp$hours, temp[,o],
                       span = 0.6, evaluation = 25, degree = 2,
                       family = 'gaussian')
    data.som2 <- cbind(data.som2,lo$y)
    col.names <- cbind(col.names,paste(a,o,sep = ','))
    hh <- c(hh, lo$x)
  }
}
data.som2 <- data.frame(data.som2)
colnames(data.som2) <- col.names
data.som2 <- melt(data.som2, variable.name = 'sample', value.name = 'fitness')
temp <- str_split(data.som2$sample, ',', simplify = T)
colnames(temp) <- c('arm','orf_name')
data.som2 <- cbind(temp, hours = hh, data.som2)

# data.som2 <- data.som[data.som$arm == 'SD-Met-Cys+Gal',]
# data.som2 <- data.som2[data.som2$hours %in% unique(data.som2$hours)[c(T,F)],]
# data.som2 <- melt(data.som2, id.vars = c('arm','stage','hours','n'),
#                   variable.name = 'orf_name', value.name = 'fitness')

head(data.som2)

## SOM ANALYSIS
temp <- data.som2 %>% filter(arm == 'SD-Met-Cys+Gal')
temp <- temp[,c('orf_name','fitness')]
data_matrix <- NULL
for (o in unique(temp$orf_name)) {
  data_matrix <- rbind(data_matrix, temp$fitness[temp$orf_name == o])
}
rownames(data_matrix) <- unique(temp$orf_name)
data_matrix <- data.frame(data_matrix)

som_data <- run_soms_analysis( data_matrix, 4, 3, 7)
som_data$pivot$hours <- data.som2$hours[data.som2$arm == 'SD-Met-Cys+Gal']

som_clusters <- som_data$pivot %>%
  group_by(cluster, id_a) %>%
  summarise(cluster = max(cluster), id_a = max(id_a)) %>%
  data.frame()
som_clusters <- merge(som_clusters, bitr(som_clusters$id_a, fromType = "ORF",
                                             toType = c("GENENAME","DESCRIPTION"),
                                             OrgDb = org.Sc.sgd.db), by.x = 'id_a', by.y = 'ORF', all = T)
write.csv(som_clusters, file = sprintf('%s/%s/SOM_results.csv', res_path, expt.name))

som_clusters %>%
  group_by(cluster) %>% count() %>% data.frame()

##### PLOTTING THE SOM RESULTS
all.som <- ggplot(som_data$pivot[som_data$pivot$cluster > 0,],
                  aes(x = hours,
                      y = y)) +
  geom_vline(xintercept = c(69,137), lwd = 0.2, linetype = 'dashed') +
  # geom_rect(mapping=aes(xmin=0, xmax=69, ymin=-20, ymax=25, fill= 'Pre-Screen #1', alpha= 'Pre-Screen #1')) +
  # geom_rect(mapping=aes(xmin=69, xmax=137, ymin=-20, ymax=25, fill= 'Pre-Screen #2', alpha= 'Pre-Screen #2')) +
  # geom_rect(mapping=aes(xmin=137, xmax=159, ymin=-20, ymax=25, fill= 'Final Screen', alpha= 'Final Screen')) +
  # geom_line(aes(group = id_a), col = '#BDBDBD', alpha = 0.3) +
  # geom_line(data = som_data$pivot[som_data$pivot$cluster > 0 &
  #                                   som_data$pivot$id_a %in% unique(data.sul$orf_name),],
  #           aes(x = data.som2$hours[data.som2$arm == 'SD-Met-Cys+Gal' &
  #                                     data.som2$orf_name %in% unique(data.sul$orf_name)],
  #               y = y, group = id_a), col = 'blue', alpha = 0.3) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               aes(group = 1), geom="ribbon", alpha = 0.4) +
  geom_smooth(data = data.lim2[data.lim2$arm == 'SD-Met-Cys+Gal',],
            aes(x = hours, y = fitness_ul), method = 'loess', se = F,
            col = 'red', linetype = 'dashed', lwd = 0.5) +
  geom_smooth(data = data.lim2[data.lim2$arm == 'SD-Met-Cys+Gal',],
              aes(x = hours, y = fitness_ll), method = 'loess', se = F,
              col = 'red', linetype = 'dashed', lwd = 0.5) +
  stat_summary(fun=mean, aes(group=1), geom="line", colour="blue") +
  scale_y_continuous(trans = 'pseudo_log') +
  labs(title = 'Fitness Dynamics',
       x = 'Time',
       y = 'Fitness') +
  # scale_fill_discrete(limits = c('Pre-Screen #1','Pre-Screen #2','Final Screen')) +
  # scale_alpha_manual(guide = F, values = c(0.3,0.3,0.3)) +
  facet_wrap(.~cluster, nrow = 4) +
  theme_linedraw() +
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        # axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position = 'bottom',
        legend.margin = margin(0.1,0.1,0.1,0.1, "mm"),
        strip.text = element_text(size = 5,
                                  margin = margin(0.1,0.1,0.1,0.1, "mm"))) +
  coord_cartesian(ylim = c(0,20))
ggsave(sprintf('%s/%s/SOM_ALL.png',fig_path, expt.name), all.som,
       width = two.c, height = one.5c, units = 'mm',
       dpi = 600)

##### GO ENRICHMENT OF CLUSTERS
som.enrich <- gkEnrichSOM(som_data$pivot)
write.csv(som.enrich, file = sprintf('%s/%s/SOM_enrichments.csv', res_path, expt.name))


##### REFERENCE DYNAMIC
data.ref <- data.lim
data.ref$hours[data.ref$arm == 'SD+Met-Cys+Gal' & data.ref$stage == 'Final Screen' & data.ref$hours > 0] <- 
  data.ref$hours[data.ref$arm == 'SD+Met-Cys+Gal' & data.ref$stage == 'Final Screen' & data.ref$hours > 0] - 1

data.ref <- merge(data.ref[data.ref$arm == 'SD-Met-Cys+Gal',c('stage','hours','cs_m')], 
      data.ref[data.ref$arm == 'SD+Met-Cys+Gal',c('stage','hours','cs_m')],
      by = c('stage','hours'), suffixes = c('_MM','_PM'))
data.ref$cs_m_rel <- data.ref$cs_m_MM/data.ref$cs_m_PM

plot.ref.dyn <- data.ref %>%
  ggplot(aes(x = hours, y = cs_m_rel)) +
  stat_summary(fun=mean, aes(group=1), geom="line", colour="blue") +
  geom_point() +
  labs(y = 'Relative Colony Size (-Met/+Met)',
       x = 'Time (hours)') +
  facet_wrap(.~stage) +
  theme_linedraw() +
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        # axis.text.x = element_text(angle = 90, size = 7, vjust = 0.5),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.position = 'bottom',
        legend.margin = margin(0.1,0.1,0.1,0.1, "mm"),
        strip.text = element_text(size = 7,
                                  margin = margin(0.1,0.1,0.1,0.1, "mm")))
ggsave(sprintf('%s/%s/REFERENCE_DYNAMICS.png',fig_path, expt.name), plot.ref.dyn,
       height = one.c, width = two.c, units = 'mm',
       dpi = 600)
