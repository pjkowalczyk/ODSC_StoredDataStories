library(rcdk)
library(tidyverse)
library(magrittr)
library(purrr)
library(stringr)
library(caret)
library(corrplot)
library(ggplot2)
library(ggthemes)

## read data

fileList <- dir('data')
for (j in c(fileList[3], fileList[11], fileList[16], fileList[24])) {
  
  inputFilename <- paste0('data/', j)
  mols <- load.molecules(inputFilename)
  
  for (i in 2:length(mols)) {
    q <- data.frame(t(unlist(get.properties(mols[[i]]))))
    Q <- rbind(Q, q)
  }
  
  descNames <-
    unique(unlist(sapply(get.desc.categories(), get.desc.names)))
  
  descs <- eval.desc(mols, descNames)
  
  alles <- cbind(Q, descs)
  
  outputFilename <-
    paste0('cache/', str_sub(j, 1, str_locate_all(j, '_')[[1]][2] - 1), '_desc.csv')
  write.csv(alles, outputFilename, row.names = FALSE)
  
}

## descriptor curation

### BioHL

#### training set
train_BioHL <- read.csv('cache/TR_BioHL_desc.csv', header = TRUE, stringsAsFactors = FALSE)
train_BioHL_meta <- train_BioHL[ , 1:13]
train_BioHL_y <- train_BioHL[ , 14]
train_BioHL_x <- train_BioHL[ , 15:ncol(train_BioHL)]

nzv <- nearZeroVar(train_BioHL_x, freqCut = 100/0)
train_BioHL_x <- train_BioHL_x[ , -nzv]

trans <- preProcess(train_BioHL_x, method = c('center', 'scale'))
train_BioHL_x <- predict(trans, train_BioHL_x)

correlations <- cor(train_BioHL_x)
corrplot(correlations, order = 'hclust')
highCorr <- findCorrelation(correlations, cutoff = 0.85)
train_BioHL_x <- train_BioHL_x[ , -highCorr]

#### test set

test_BioHL <- read.csv('cache/TST_BioHL_desc.csv', header = TRUE, stringsAsFactors = FALSE)
test_BioHL_meta <- test_BioHL[ , 1:13]
test_BioHL_y <- test_BioHL[ , 14]
test_BioHL_x <- test_BioHL[ , 15:ncol(test_BioHL)]

test_BioHL_x <- test_BioHL_x[ , -nzv]

test_BioHL_x <- predict(trans, test_BioHL_x)

test_BioHL_x <- test_BioHL_x[ , -highCorr]

### PCA

pca <- preProcess(train_BioHL_x, method = c('pca'))
train_BioHL_x_pca <- predict(pca, train_BioHL_x)
test_BioHL_x_pca <- predict(pca, test_BioHL_x)

train_pca <- train_BioHL_x_pca %>%
  select(PC1, PC2) %>%
  mutate(dataset = 'train')
test_pca <- test_BioHL_x_pca %>%
  select(PC1, PC2) %>%
  mutate(dataset = 'test')
pcaPts <- rbind(train_pca, test_pca)

p <-
  ggplot(pcaPts, aes(PC1, PC2)) +
  geom_point(aes(colour = factor(dataset), shape = factor(dataset))) +
  ggthemes::theme_tufte()
p

### k-nearest neighbor regression

fit <- knnreg(train_BioHL_x, train_BioHL_y, k = 5)
result <- data.frame(test_BioHL_y, predict(fit, test_BioHL_x))
colnames(result) <- c('observed', 'predicted')

p <-
  ggplot(result, aes(observed, predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  xlim(c(0, 3)) + ylim(c(0, 3)) +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0, slope = 1, colour = 'red')
p
