library(tidyverse)
library(magrittr)
library(caret)
library(ggplot2)
library(ggthemes)
library(jtools)
library(randomForest)

# read data
df <-
  read.csv('data/water_solubility.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  na.omit()

# view sample of data
head(df[sample(nrow(df), 10), ])

LogMolar <-
  ggplot(df, aes(LogMolar, stat(density))) +
  # geom_freqpoly(binwidth = 0.25, size = 1) +
  geom_histogram(binwidth = 0.25, color = 'white', fill = 'blue') +
  theme(legend.position = "none") +
  ggthemes::theme_hc()
LogMolar

ggsave('graphics/WS_LogMolar_Histogram.jpg', plot = LogMolar)

# train / text split
inTrain <- caret::createDataPartition(df$LogMolar, p = 0.8, list = FALSE)
train <- df[inTrain, ]
test <- df[-inTrain, ]

X_train <- train[ , 3:ncol(train)]
y_train <- train[ , 2] %>% data.frame()
colnames(y_train) <- 'LogMolar'
X_test <- test[ , 3:ncol(test)]
y_test <- test[ , 2] %>% data.frame()
colnames(y_test) <- 'LogMolar'

TRAIN <- train %>%
  mutate(set = 'train')
TEST <- test %>%
  mutate(set = 'test')

LogMolar <- rbind(TRAIN, TEST)

LogMolar_train_test <-
  ggplot(LogMolar, aes(LogMolar, stat(density), colour = set)) +
  geom_freqpoly(binwidth = 0.25, size = 1) +
  scale_color_manual(values = c('#EB6B4A', '#0B3087')) +
  theme(legend.position = "none") +
  ggthemes::theme_hc()
LogMolar_train_test

ggsave('graphics/WS_LogMolar_TrainTest.jpg', plot = LogMolar_train_test)

# curate data
dim(X_train)

## near-zero variance descriptors

nzv <- caret::nearZeroVar(X_train, freqCut = 100/0)
X_train <- X_train[ , -nzv]
### and
X_test <- X_test[ , -nzv]

dim(X_train)

## highly correlated descriptors

correlations <- cor(X_train)
corrplot::corrplot(correlations, order = 'hclust')
#
jpeg('graphics/WS_fullCorrelation.jpg')
corrplot::corrplot(correlations, order = 'hclust')
dev.off()
#
highCorr <- findCorrelation(correlations, cutoff = 0.85)
X_train <- X_train[ , -highCorr]
### and
X_test <- X_test[ , -highCorr]

correlations <- cor(X_train)
corrplot::corrplot(correlations, order = 'hclust')
#
jpeg('graphics/WS_reducedCorrelation.jpg')
corrplot::corrplot(correlations, order = 'hclust')
dev.off()
#

dim(X_train)

## linear combinations

comboInfo <- findLinearCombos(X_train)
X_train <- X_train[ , -comboInfo$remove]
### and
X_test <- X_test[ , -comboInfo$remove]

dim(X_train)

## center & scale descriptors

preProcValues <- preProcess(X_train, method = c("center", "scale"))

X_trainTransformed <- predict(preProcValues, X_train)
### and
X_testTransformed <- predict(preProcValues, X_test)

### PCA

pca <- preProcess(X_trainTransformed, method = c('pca'))
X_train_pca <- predict(pca, X_trainTransformed)
X_test_pca <- predict(pca, X_testTransformed)

train_pca <- X_train_pca %>%
  select(PC1, PC2) %>%
  mutate(dataset = 'train')
test_pca <- X_test_pca %>%
  select(PC1, PC2) %>%
  mutate(dataset = 'test')
pcaPts <- rbind(train_pca, test_pca)

LogMolar_PC <-
  ggplot(pcaPts, aes(PC1, PC2, colour = factor(dataset))) +
  geom_point(aes(shape = factor(dataset))) +
  labs(title = 'WS PCA') +
  scale_color_manual(values = c('#EB6B4A', '#0B3087')) +
  theme(legend.position="none") +
  ggthemes::theme_tufte()
LogMolar_PC

ggsave('graphics/WS_LogMolar_PC.jpg', plot = LogMolar_PC)

# models

fitControl <- trainControl(## 5-fold CV
  method = "repeatedcv",
  repeats = 5)

set.seed(42)

## multiple linear regression

trainSet <- cbind(y_train, X_trainTransformed)

mlr <- train(LogMolar ~ .,
             data = trainSet,
             method = 'lm',
             trControl = fitControl)

y_predict <- predict(mlr, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(LogMolar ~ Predicted, data = data2plot))
summ(lm(LogMolar ~ Predicted, data = data2plot))

LogMolar_mlr <-
  ggplot(data2plot, aes(Predicted, LogMolar)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Multiple Linear Regression\n test dat') +
  # subtitle = 'Multiple Linear Regression\n test data') +
  ggthemes::theme_tufte()
LogMolar_mlr <- LogMolar_mlr + geom_abline(intercept = 0,
                                           slope = 1,
                                           colour = 'red')
LogMolar_mlr

ggsave('graphics/WS_LogMolar_mlrTest.jpg', plot = LogMolar_mlr)

data2plot$res <- resid(lm(LogMolar ~ Predicted, data = data2plot))

LogMolar_mlr_res <-
  ggplot(data2plot, aes(LogMolar, res)) +
  geom_point(colour = "blue", size = 2) +
  # coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  # geom_smooth(method = 'lm') +
  labs(title = 'Residual Plot') +
  # subtitle = 'Multiple Linear Regression\n test data') +
  geom_hline(yintercept = 0,
             color = "red",
             size = 1.5) +
  ggthemes::theme_tufte()

LogMolar_mlr_res

ggsave('graphics/WS_LogMolar_mlrTestRes.jpg', plot = LogMolar_mlr_res)

y_predict <- predict(mlr, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

# mlrPR <- postResample(pred = y_predict, obs = X_trainTransformed)
# rmse_train = c(mlrPR[1])
# r2_train = c(mlrPR[2])

data2plot <- cbind(y_train, y_predict)

summary(lm(LogMolar ~ Predicted, data = data2plot))
summ(lm(LogMolar ~ Predicted, data = data2plot))

p <-
  ggplot(data2plot, aes(Predicted, LogMolar)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogMolar',
       subtitle = 'Multiple Linear Regression\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

ggsave('graphics/WS_LogMolar_mlrTrain.jpg', plot = p)

data2plot$res <- resid(lm(LogMolar ~ Predicted, data = data2plot))

LogMolar_mlr_res <-
  ggplot(data2plot, aes(LogMolar, res)) +
  geom_point(colour = "blue", size = 2) +
  # coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  # geom_smooth(method = 'lm') +
  labs(title = 'Residual Plot') +
  # subtitle = 'Multiple Linear Regression\n test data') +
  geom_hline(yintercept = 0,
             color = "red",
             size = 1.5) +
  ggthemes::theme_tufte()

LogMolar_mlr_res

ggsave('graphics/WS_LogMolar_mlrTrainRes.jpg', plot = LogMolar_mlr_res)

## k-nearest neighbors

knnModel <- train(
  LogMolar ~ .,
  data = trainSet,
  method = 'knn',
  tuneGrid = data.frame(.k = 1:20),
  trControl = fitControl
)

y_predict <- predict(knnModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(LogMolar ~ Predicted, data = data2plot))

LogMolar_kNN <-
  ggplot(data2plot, aes(Predicted, LogMolar)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'K-Nearest Neighbors\n test data') +
  # subtitle = 'K-Nearest Neighbors\n test data') +
  ggthemes::theme_tufte()
LogMolar_kNN <- LogMolar_kNN + geom_abline(intercept = 0,
                                           slope = 1,
                                           colour = 'red')
LogMolar_kNN

ggsave('graphics/WS_LogMolar_knnTest.jpg', plot = LogMolar_knn)

data2plot$res <- resid(lm(LogMolar ~ Predicted, data = data2plot))

LogMolar_knn_res <-
  ggplot(data2plot, aes(LogMolar, res)) +
  geom_point(colour = "blue", size = 2) +
  # coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  # geom_smooth(method = 'lm') +
  labs(title = 'Residual Plot') +
  # subtitle = 'Multiple Linear Regression\n test data') +
  geom_hline(yintercept = 0,
             color = "red",
             size = 1.5) +
  ggthemes::theme_tufte()

LogMolar_knn_res

ggsave('graphics/WS_LogMolar_knnTestRes.jpg', plot = LogMolar_mlr_res)

y_predict <- predict(knnModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(LogMolar ~ Predicted, data = data2plot))

p <-
  ggplot(data2plot, aes(Predicted, LogMolar)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogMolar',
       subtitle = 'K-Nearest Neighbors\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

ggsave('graphics/WS_LogMolar_knnTrain.jpg', plot = p)

data2plot$res <- resid(lm(LogMolar ~ Predicted, data = data2plot))

LogMolar_knn_res <-
  ggplot(data2plot, aes(LogMolar, res)) +
  geom_point(colour = "blue", size = 2) +
  # coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  # geom_smooth(method = 'lm') +
  labs(title = 'Residual Plot') +
  # subtitle = 'Multiple Linear Regression\n test data') +
  geom_hline(yintercept = 0,
             color = "red",
             size = 1.5) +
  ggthemes::theme_tufte()

LogMolar_knn_res

ggsave('graphics/WS_LogMolar_knnTrainRes.jpg', plot = LogMolar_knn_res)

## Random Forests

rf <- randomForest(LogMolar ~ ., data = trainSet)

y_predict <- predict(rf, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(LogMolar ~ Predicted, data = data2plot))

LogMolar_rf <-
  ggplot(data2plot, aes(Predicted, LogMolar)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Random Forest\n test data') +
  # subtitle = 'K-Nearest Neighbors\n test data') +
  ggthemes::theme_tufte()
LogMolar_rf <- LogMolar_rf + geom_abline(intercept = 0,
                                           slope = 1,
                                           colour = 'red')
LogMolar_rf

ggsave('graphics/WS_LogMolar_rfTest.jpg', plot = LogMolar_rf)

data2plot$res <- resid(lm(LogMolar ~ Predicted, data = data2plot))

LogMolar_rf_res <-
  ggplot(data2plot, aes(LogMolar, res)) +
  geom_point(colour = "blue", size = 2) +
  # coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  # geom_smooth(method = 'lm') +
  labs(title = 'Residual Plot') +
  # subtitle = 'Multiple Linear Regression\n test data') +
  geom_hline(yintercept = 0,
             color = "red",
             size = 1.5) +
  ggthemes::theme_tufte()

LogMolar_rf_res

ggsave('graphics/WS_LogMolar_rfTestRes.jpg', plot = LogMolar_rf_res)

y_predict <- predict(rf, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(LogMolar ~ Predicted, data = data2plot))

p <-
  ggplot(data2plot, aes(Predicted, LogMolar)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogMolar',
       subtitle = 'Random Forest\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

ggsave('graphics/WS_LogMolar_rfTrain.jpg', plot = p)

data2plot$res <- resid(lm(LogMolar ~ Predicted, data = data2plot))

LogMolar_rf_res <-
  ggplot(data2plot, aes(LogMolar, res)) +
  geom_point(colour = "blue", size = 2) +
  # coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  # geom_smooth(method = 'lm') +
  labs(title = 'Residual Plot') +
  # subtitle = 'Multiple Linear Regression\n test data') +
  geom_hline(yintercept = 0,
             color = "red",
             size = 1.5) +
  ggthemes::theme_tufte()

LogMolar_rf_res

ggsave('graphics/WS_LogMolar_rfTrainRes.jpg', plot = LogMolar_rf_res)
