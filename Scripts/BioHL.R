library(rcdk)
library(tidyverse)
library(magrittr)
library(purrr)
library(stringr)
library(caret)
library(corrplot)
library(ggplot2)
library(ggthemes)
library(egg)

# read data

## training data
train <-
  read.csv('cache/TR_BioHL_112_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,-CAS,-ROMol,-SMILES,-ID) %>%
  select(LogHalfLife, everything()) %>%
  na.omit()

X_train <- train %>%
  select(-LogHalfLife)
y_train <- train %>%
  select(LogHalfLife) %>%
  data.frame()

## test data
test <-
  read.csv('cache/TST_BioHL_38_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,-CAS,-ROMol,-SMILES,-ID) %>%
  select(LogHalfLife, everything()) %>%
  na.omit()

X_test <- test %>%
  select(-LogHalfLife)
y_test <- test %>%
  select(LogHalfLife) %>%
  data.frame()

# curate data

## near-zero variance descriptors

nzv <- nearZeroVar(X_train, freqCut = 100/0)
X_train <- X_train[ , -nzv]
### and
X_test <- X_test[ , -nzv]

## highly correlated descriptors

correlations <- cor(X_train)
corrplot::corrplot(correlations, order = 'hclust', tl.cex = 0.6)
highCorr <- findCorrelation(correlations, cutoff = 0.85)
X_train <- X_train[ , -highCorr]
### and
X_test <- X_test[ , -highCorr]

correlations <- cor(X_train)
corrplot::corrplot(correlations, order = 'hclust', tl.cex = 0.8)

## linear combinations

comboInfo <- findLinearCombos(X_train) # returns NULL
# X_train <- X_train[ , -comboInfo$remove]
# ### and
# X_test <- X_test[ , -nzv]

## center & scale descriptors

preProcValues <- preProcess(X_train, method = c("center", "scale"))

X_trainTransformed <- predict(preProcValues, X_train)
### and
X_testTransformed <- predict(preProcValues, X_test)

### PCA

# pca <- preProcess(X_trainTransformed, method = c('pca'))
# X_train_pca <- predict(pca, X_trainTransformed)
# X_test_pca <- predict(pca, X_testTransformed)
# 
# train_pca <- X_train_pca %>%
#   select(PC1, PC2) %>%
#   mutate(dataset = 'train')
# test_pca <- X_test_pca %>%
#   select(PC1, PC2) %>%
#   mutate(dataset = 'test')
# pcaPts <- rbind(train_pca, test_pca)
# 
# p <-
#   ggplot(pcaPts, aes(PC1, PC2)) +
#   geom_point(aes(colour = factor(dataset), shape = factor(dataset))) +
#   ggthemes::theme_tufte()
# p

# models

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  repeats = 5)

set.seed(350)

## multiple linear regression

trainSet <- cbind(y_train, X_trainTransformed)

mlr <- train(LogHalfLife ~ .,
             data = trainSet,
             method = 'lm',
             trControl = fitControl)

y_predict <- predict(mlr, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogHalfLife, data = data2plot))

pTest <-
  ggplot(data2plot, aes(LogHalfLife, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Biodegradability half-life',
       subtitle = 'Multiple Linear Regression\n test data') +
  ggthemes::theme_tufte()
pTest <- pTest + geom_abline(intercept = 0,
                             slope = 1,
                             colour = 'red')
pTest

y_predict <- predict(mlr, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

mlrPR <- postResample(pred = y_predict, obs = X_trainTransformed)
rmse_train = c(mlrPR[1])
r2_train = c(mlrPR[2])

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogHalfLife, data = data2plot))

pTrain <-
  ggplot(data2plot, aes(LogHalfLife, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Biodegradability half-life',
       subtitle = 'Multiple Linear Regression\n training data') +
  ggthemes::theme_tufte()
pTrain <- pTrain + geom_abline(intercept = 0,
                               slope = 1,
                               colour = 'red')
pTrain

mlr_figure <- egg::ggarrange(pTrain, pTest, nrow = 1)
mlr_figure

## partial least squares

plsModel <- train(
  LogHalfLife ~ .,
  data = trainSet,
  method = 'pls',
  tuneLength = 20,
  trControl = fitControl
)

y_predict <- predict(plsModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogHalfLife, data = data2plot))

pTest <-
  ggplot(data2plot, aes(LogHalfLife, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Biodegradability half-life',
       subtitle = 'Partial Least Squares\n test data') +
  ggthemes::theme_tufte()
pTest <- pTest + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
pTest

y_predict <- predict(plsModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogHalfLife, data = data2plot))

pTrain <-
  ggplot(data2plot, aes(LogHalfLife, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'Biodegradability half-life',
       subtitle = 'Partial Least Squares\n training data') +
  ggthemes::theme_tufte()
pTrain <- pTrain + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
pTrain

pls_figure <- egg::ggarrange(pTrain, pTest, nrow = 1)
pls_figure

## support vector machines

svmModel <- train(
  LogHalfLife ~ .,
  data = trainSet,
  method = 'svmRadial',
  # tuneLength = 14,
  trControl = fitControl
)

y_predict <- predict(svmModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogHalfLife, data = data2plot))

pTest <-
  ggplot(data2plot, aes(LogHalfLife, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Biodegradability half-life',
       subtitle = 'Support Vector Machines\n test data') +
  ggthemes::theme_tufte()
pTest <- pTest + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
pTest

y_predict <- predict(svmModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogHalfLife, data = data2plot))

pTrain <-
  ggplot(data2plot, aes(LogHalfLife, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'Biodegradability half-life',
       subtitle = 'Support Vector Machines\n training data') +
  ggthemes::theme_tufte()
pTrain <- pTrain + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
pTrain

svm_figure <- egg::ggarrange(pTrain, pTest, nrow = 1)
svm_figure

## k-nearest neighbors

knnModel <- train(
  LogHalfLife ~ .,
  data = trainSet,
  method = 'knn',
  tuneGrid = data.frame(.k = 1:20),
  trControl = fitControl
)

y_predict <- predict(knnModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogHalfLife, data = data2plot))

pTest <-
  ggplot(data2plot, aes(LogHalfLife, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Biodegradability half-life',
       subtitle = 'K-Nearest Neighbors\n test data') +
  ggthemes::theme_tufte()
pTest <- pTest + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
pTest

y_predict <- predict(knnModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogHalfLife, data = data2plot))

pTrain <-
  ggplot(data2plot, aes(LogHalfLife, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'Biodegradability half-life',
       subtitle = 'K-Nearest Neighbors\n training data') +
  ggthemes::theme_tufte()
pTrain <- pTrain + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
pTrain

knn_figure <- egg::ggarrange(pTrain, pTest, nrow = 1)
knn_figure

## boosted trees

gbmGrid <- expand.grid(
  .interaction.depth = seq(1, 7, by = 2),
  .n.trees = seq(100, 1000, by = 50),
  .shrinkage = c(0.001, 0.1),
  .n.minobsinnode = 3
)

treeModel <- train(
  LogHalfLife ~ .,
  data = trainSet,
  method = 'gbm',
  tuneGrid = gbmGrid,
  verbose = FALSE
)

y_predict <- predict(treeModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogHalfLife, data = data2plot))

pTest <-
  ggplot(data2plot, aes(LogHalfLife, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Biodegradability half-life',
       subtitle = 'Boosted Trees\n test data') +
  ggthemes::theme_tufte()
pTest <- pTest + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
pTest

y_predict <- predict(treeModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogHalfLife, data = data2plot))

pTrain <-
  ggplot(data2plot, aes(LogHalfLife, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'Biodegradability half-life',
       subtitle = 'Boosted Trees\n training data') +
  ggthemes::theme_tufte()
pTrain <- pTrain + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
pTrain

gbm_figure <- egg::ggarrange(pTrain, pTest, nrow = 1)
gbm_figure
