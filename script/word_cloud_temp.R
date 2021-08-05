library("tm")
library("SnowballC")
library("wordcloud")
library("RColorBrewer")

# Load the data as a corpus
docs <- Corpus(VectorSource(data.diff2$DESCRIPTION))

##### TEXT CLEANUP
# Convert the text to lower case
docs <- tm_map(docs, content_transformer(tolower))
# Remove numbers
docs <- tm_map(docs, removeNumbers)
# Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))
# Remove punctuations
docs <- tm_map(docs, removePunctuation)
# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)

##### GENERATE TERM-DOCUMENT MATRIX
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
head(d, 20)


docs <- Corpus(VectorSource(data.diff2$DESCRIPTION[data.diff2$fitness_diff >= 0.4 & data.diff2$stage == 'Final Screen']))

##### TEXT CLEANUP
# Convert the text to lower case
docs <- tm_map(docs, content_transformer(tolower))
# Remove numbers
docs <- tm_map(docs, removeNumbers)
# Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))
# Remove punctuations
docs <- tm_map(docs, removePunctuation)
# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)

##### GENERATE TERM-DOCUMENT MATRIX
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d2 <- data.frame(word = names(v),freq=v)
head(d2, 20)

head(d2[!(d2$word %in% d$word[d$freq >= 1500]),],20)

set.seed(1234)
wordcloud(words = d2$word[!(d2$word %in% d$word[d$freq >= 1500])], freq = d$freq, min.freq = 10,
          max.words=300, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))

d3 <- merge(d,d2,by = 'word')
d3$freq <- d3$freq.y/d3$freq.x

ggplot(d3, aes(x = freq)) +
  geom_line(stat = 'density')

d3[d3$freq > 0.1 & d3$freq.y > 50,]
