
library(growthcurver)

strain.labs.h2s <- read.csv(file = 'data/h2s/FeEDTA_samples.csv')
data.h2s <- read.csv(file = 'data/h2s/FeEDTA_readings.csv')


data.pred.h2s <- NULL
t <- data.h2s$Time
for (c in colnames(data.h2s)[2:dim(data.h2s)[2]]) {
  temp <- data.h2s[c]
  temp[temp <= 0] <- temp[temp <= 0] + 0.0001
  lo <- loess.smooth(t, log(temp),
                     span = 0.6, evaluation = 50, degree = 2,
                     family = 'gaussian')
  data.pred.h2s <- cbind(data.pred.h2s,exp(lo$y))
}
data.pred.h2s <- cbind(lo$x, data.pred.h2s)
colnames(data.pred.h2s) <- colnames(data.h2s)
data.pred.h2s <- data.frame(data.pred.h2s)
head(data.pred.h2s)

data.h2s.gc.res <- SummarizeGrowthByPlate(data.pred.h2s)

data.pred.h2s %>%
  melt(id.vars = 'Time', variable.name = 'Flask', value.name = 'OD') %>%
  ggplot(aes(x = Time, y = OD, col = Flask)) +
  geom_line()
