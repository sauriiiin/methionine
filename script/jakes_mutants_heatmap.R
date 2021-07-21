library(Hmisc)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(PoiClaClu)
library(RColorBrewer)

head(data)

temp <- data %>%
  group_by(attempt, condition, time, orf_name, bio_rep) %>%
  summarize(relative_cs = median(relative_cs, na.rm = T)) %>%
  data.frame()

corr.dat <- temp[temp$attempt == 'pilot' & temp$condition == 'SD-MET+Glucose',c(3,4,5)]
col.names <- NULL
corr.temp <- NULL
for (a in unique(temp$attempt)) {
  for(c in unique(temp$condition[temp$attempt == a])) {
    corr.temp  <- cbind(corr.temp, temp$relative_cs[temp$attempt == a & temp$condition == c])
    col.names <- c(col.names, paste(a, c))
  }
}
colnames(corr.temp) <- col.names
corr.dat <- cbind(corr.dat, corr.temp)
head(corr.dat)

mat <- cor(corr.dat[,col.names], use = 'complete.obs', method = 'spearman')
mat  

col_fun = colorRamp2(c(0.3, 0.5, 1), c("white", "#FF5252", "#512DA8"))
ha <- HeatmapAnnotation(
  Replicate = str_split(col.names, ' ', simplify = T)[,1],
  Condition = str_split(col.names, ' ', simplify = T)[,2],
  col = list(Replicate = c("pilot" = "Navy", "copy1" = "#607D8B", "copy2" = "#FFA000"),
             Condition = c("YPDA" = "#B3E5FC", "SD-MET+Glucose" = "#448AFF")),
  gp = gpar(col = "black", lwd = 0.2)
)

ra <- rowAnnotation(
  Replicate = str_split(col.names, ' ', simplify = T)[,1],
  Condition = str_split(col.names, ' ', simplify = T)[,2],
  col = list(Replicate = c("pilot" = "Navy", "copy1" = "#607D8B", "copy2" = "#FFA000"),
             Condition = c("YPDA" = "#B3E5FC", "SD-MET+Glucose" = "#448AFF")),
  gp = gpar(col = "black", lwd = 0.2),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

hm <- Heatmap(mat,
              name = "Correlation",
              column_title = "Jake's Mutant Experiment",
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%0.2f", mat[i, j]), x, y, gp = gpar(fontsize = 4))
              },
              jitter = T,
              column_title_gp = gpar(fontsize = 12, fontface = "bold"),
              show_row_dend = F,
              rect_gp = gpar(col = "black", lwd = 0.2),
              show_row_names = F,
              show_column_names = F,
              border = T,
              top_annotation = ha)
draw(hm+ra)
