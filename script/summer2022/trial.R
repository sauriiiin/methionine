##### COLONY BASED FITNESS SCREEN ANALYSIS
##### Author : Saurin Parikh
##### Email  : dr.saurin.parikh@gmail.com
##### Date   : 07/07/2022

##### INITIALIZE
## you might need to install some of these libraries either from CRAN or from Bioconductor
## to install from Bioconductor install BiocManager and then used command BiocManager::install("package name")
library(RMariaDB)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(genefilter)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(stringr)

## Connecting to MYSQL
# enter your mysql database, username and password
conn <- dbConnect(MariaDB(), dbname = '',         #your mysql database name
                  usr = '',                       #your mysql username
                  password = '',                  #your mysql password
                  host = 'paris.csb.pitt.edu')

## Getting a list of all the proto-genes from MYSQL database
pgs <- dbGetQuery(conn, 'select orf_name from saurin_test.PROTOGENES where pg_2012 = 1')

##### SETTING FIGURE SIZE
one.c <- 90   #single column
one.5c <- 140 #1.5 column
two.c <- 190  #full width

##### SETIING TEXT SIZE
titles <- 8
txt <- 8
lbls <- 9

##### GATHER DATA
# database name, condition, stage, colony density, bag, fitness tablename
tables <- rbind(c('saurin_test','MM','FS1',6144,'No','MET_OE_FS_MM_6144_FITNESS_STATS', 'MET_OE_FS_MM_6144_PVALUE'),
            c('Brandon','MM','FS1',6144,'No','bg_MET_OE2206_FS1_MM_nobag_6144_FITNESS_STATS', 'bg_MET_OE2206_FS1_MM_nobag_6144_PVALUE'),
            c('aaron_z','PM','FS1',6144,'No','aaz_MET_OE2206_FS1_PM_nobag_6144_FITNESS_STATS','aaz_MET_OE2206_FS1_PM_nobag_6144_PVALUE')) %>%
  data.frame()
# update above table to add more conditions in the same format as others
colnames(tables) <- c('database','condition','stage','density','bag','fit_tablename','pval_tablename')
head(tables)

# loop through all the table information you provided above to now gather all the data from MYSQL
dim(tables)
fit_data <- NULL
for (i in seq(1,dim(tables)[1],1)) {
  temp <- dbGetQuery(conn, sprintf('select a.*, b.stat, b.p from %s.%s a, %s.%s b
                                   where a.hours = b.hours and a.strain_id = b.strain_id',
                                   tables$database[i], tables$fit_tablename[i],
                                   tables$database[i], tables$pval_tablename[i]))
  # the final time point in the table is set as the saturation time point
  temp$saturation <- max(unique(temp$hours))
  
  # entering all the information from the tables table from above
  temp$database <- tables$database[i]
  temp$condition <- tables$condition[i]
  temp$stage <- tables$stage[i]
  temp$density <- tables$density[i]
  temp$bag <- tables$bag[i]
  temp$tablename <- tables$fit_tablename[i]
  
  # making a master table with data from all the MYSQL tables
  fit_data <- rbind(fit_data, temp)
}
# filtering data to work with only the saturation time point data for now
fit_data <- fit_data %>% filter(hours == saturation)

# categorizing mutants into proto-gene, gene or reference using the pgs we downloaded earlier
fit_data$orf_type[fit_data$orf_name %in% pgs$orf_name] <- 'Proto-gene'
fit_data$orf_type[fit_data$orf_name == 'BF_control'] <- 'Reference'
fit_data$orf_type[is.na(fit_data$orf_type)] <- 'Gene'

# categorizing the phenotypes of the mutants
fit_data$phenotype[fit_data$p <= 0.05 & fit_data$stat > 0] <- 'Beneficial'
fit_data$phenotype[fit_data$p <= 0.05 & fit_data$stat < 0] <- 'Deleterious'
fit_data$phenotype[is.na(fit_data$phenotype)] <- 'Neutral'

##### CORRELATION WITH PREVIOUS DATA
old_data <- fit_data %>% filter(tablename == 'MET_OE_FS_MM_6144_FITNESS_STATS') # filtering old minus MET data from the master table
new_data <- fit_data %>% filter(tablename == 'bg_MET_OE2206_FS1_MM_nobag_6144_FITNESS_STATS') # filter new minus MET data from the master table

# creating a single table with the old and new data
old_vs_new <- merge(old_data[,c('strain_id','orf_name','cs_median','tablename')],
                    new_data[,c('strain_id','orf_name','cs_median','tablename')],
                    by = c('strain_id','orf_name'),
                    suffixes = c('_old','_new'))

# making a correlation plot between old and new data
# each point represents a single mutant
old_vs_new %>%
  ggplot(aes(x = cs_median_new, y = cs_median_old)) +
  geom_abline(linetype = 'dashed', col = 'red', lwd = 1) +
  geom_point() +
  stat_cor(method = 'pearson') +
  labs(x = 'Fitness in current screen',
       y = 'Fitness in previous screen') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm")) +
  coord_cartesian(xlim = c(0,6),
                  ylim = c(0,6))
# find a way to save this plot


##### PLUS MET AND MINUS MET DATA
# the plus met acts as a control for the minus met data here

pm_data <- fit_data %>% filter(tablename == 'aaz_MET_OE2206_FS1_PM_nobag_6144_FITNESS_STATS') # filter new plus MET data from the master table
mm_data <- fit_data %>% filter(tablename == 'bg_MET_OE2206_FS1_MM_nobag_6144_FITNESS_STATS') # filter new minus MET data from the master table

# creating a single table with the plus met and minus met data
pm_vs_mm <- merge(pm_data[,c('strain_id','orf_name','cs_median','tablename')],
                   mm_data[,c('strain_id','orf_name','cs_median','tablename')],
                   by = c('strain_id','orf_name'),
                   suffixes = c('_pm','_mm'))

# creating a new variable in the table that is the difference in the fitness between the minus met and plus met condition per mutant
pm_vs_mm <- pm_vs_mm %>%
  mutate(cs_diff = cs_median_mm - cs_median_pm) 
# this is a measure of differential fitness giving us the exact change in fitness that happened to a mutant due to the overexpression
# of the ORF in the minus MET condition as compared to the plus MET condition
# a neutral ORF in the minus MET condition might still be considered 'useful' if the mutant had a deleterious fitness in the plus MET condition
# why do you think that is?
# can you think of other such differential fitness assessments in your dataset

# plotting the plus met vs minus met results
pm_vs_mm %>%
  ggplot(aes(x = cs_median_pm, y = cs_median_mm)) +
  geom_abline(linetype = 'dashed', col = 'red', lwd = 1) +
  geom_point() +
  geom_text(data = pm_vs_mm %>% filter(abs(cs_diff) > 0.5),
             aes(x = cs_median_pm, y = cs_median_mm + 0.2, label = orf_name),
            size = 2.5) +
  labs(x = 'Fitness in +Methionine',
       y = 'Fitness in -Methionine') +
  theme_linedraw() +
  theme(plot.title = element_text(size = titles, face = 'bold', hjust = 0.5),
        axis.title = element_text(size = titles),
        axis.text = element_text(size = txt),
        legend.title = element_text(size = titles),
        legend.text = element_text(size = txt),
        legend.position = 'bottom',
        legend.key.size = unit(3, "mm")) +
  coord_cartesian(xlim = c(0,6),
                  ylim = c(0,6))
# what do you think about the interesting hits that have been annotated in this figure?
# could you figured them out using differential fitness as well?

##### PHENOTYPE COUNTS
# overall counts of beneficial, neutral and deleterious genes and proto-genes per condition
fit_data %>%
  group_by(tablename, orf_type, phenotype) %>%
  count() %>% data.frame()


##### ODDS RATIO OF A PROTO-GENE BEING BENEFICIAL WHEN OVEREXPRESSED VS GENE BEING BENEFICIAL WHEN OVEREXPRESSED
# Figure out using numbers from above


##### GO/KEGG ENRICHMENT IN VARIOUS PHENOTYPE CATEGORIES PER CONDITION
# try to figure out what the below code is doing on your own
goe <- data.frame()
kegg <- data.frame()
for (t in unique(fit_data$tablename)) {
  allgenes <- unique(fit_data$orf_name[fit_data$tablename == t])
  allgenes <- bitr(allgenes, fromType = "ORF",
                   toType = c("ENTREZID","GENENAME","ENSEMBL"),
                   OrgDb = org.Sc.sgd.db)
  allgenes <- allgenes[!is.na(allgenes$ENSEMBL),]
  for (p in unique(fit_data$phenotype[fit_data$tablename == t])) {
    temp.deg <- fit_data$orf_name[fit_data$phenotype == p & fit_data$tablename == t]
    temp.deg <- bitr(temp.deg, fromType = "ORF",
                     toType = c("ENTREZID","GENENAME","ENSEMBL"),
                     OrgDb = org.Sc.sgd.db)
    temp.deg <- temp.deg[!is.na(temp.deg$ENSEMBL),]
    
    temp.goe <- enrichGO(gene          = temp.deg$ENSEMBL,
                         universe      = allgenes$ENSEMBL,
                         OrgDb         = org.Sc.sgd.db,
                         keyType       = "ENSEMBL",
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05)
    if (dim(temp.goe)[1] == 0) {
      cat(sprintf('There are no GO term enrichment for %s ORFs in %s.\n',p,t))
    } else{
      cat(sprintf('GO term enrichment for %s ORFs in %s are:\n%s\n',
                  p,t,
                  paste(temp.goe$Description,collapse = ', ')))
      goe <- rbind(goe, data.frame(temp.goe, tablename = t, phenotype = p))
    }
    
    temp.kegg <- enrichKEGG(gene         = temp.deg$ENSEMBL,
                            universe     = allgenes$ENSEMBL,
                            organism     = 'sce',
                            pvalueCutoff = 0.05)
    if (dim(temp.kegg)[1] == 0) {
      cat(sprintf('There are no KEGG pathway enriched for %s ORFs in %s.\n',p,t))
    } else{
      cat(sprintf('KEGG pathway enrichment for %s ORFs in %s are:\n%s\n',
                  p,t,
                  paste(temp.kegg$Description,collapse = ', ')))
      kegg <- rbind(kegg, data.frame(temp.kegg, tablename = t, phenotype = p))
    }
  }
}
goe$GeneRatio <- as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$GeneRatio,'/',simplify = T)[,2])
goe$BgRatio <- as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(goe$BgRatio,'/',simplify = T)[,2])
goe$GO <- paste0(goe$ONTOLOGY, '_', goe$Description)

kegg$GeneRatio <- as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$GeneRatio,'/',simplify = T)[,2])
kegg$BgRatio <- as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,1])/
  as.numeric(str_split(kegg$BgRatio,'/',simplify = T)[,2])

goe <- goe[order(goe$tablename,goe$phenotype,-goe$GeneRatio,-goe$Count,goe$qvalue),]
kegg <- kegg[order(kegg$tablename,kegg$phenotype,-kegg$GeneRatio,-kegg$Count,kegg$qvalue),]

goe <- merge(goe, tables[,-7], by.x = 'tablename', by.y = 'fit_tablename')
kegg <- merge(kegg, tables[,-7], by.x = 'tablename', by.y = 'fit_tablename')
# Save the GO/KEGG results a .csv or .xlsx and view them locally on your computer
# See if the enriched categories make sense to you and why?



