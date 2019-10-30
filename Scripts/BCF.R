library(tidyverse)
library(magrittr)
library(purrr)
library(stringr)
library(caret)
library(corrplot)
library(ggplot2)
library(ggthemes)

# read data

## training data
train <-
  read.csv('cache/TR_BCF_469_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X, -CAS, -ROMol, -SMILES, -ID) %>%
  select(LogBCF, everything()) %>%
  na.omit()

X_train <- train %>%
  select(-LogBCF)
y_train <- train %>%
  select(LogBCF) %>%
  data.frame()
# trainFP <-
#   read.csv('data/BCF/trainFP.csv',
#            header = TRUE,
#            stringsAsFactors = FALSE)
# TRAIN <- merge(train, trainFP, by.x = 'ChemID', by.y = 'Name')

## test data
test <-
  read.csv('cache/TST_BCF_157_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X, -CAS, -ROMol, -SMILES, -ID) %>%
  select(LogBCF, everything()) %>%
  na.omit()

X_test <- test %>%
  select(-LogBCF)
y_test <- test %>%
  select(LogBCF) %>%
  data.frame()
# testFP <-
#   read.csv('data/BCF/testFP.csv',
#            header = TRUE,
#            stringsAsFactors = FALSE)
# TEST <- merge(test, testFP, by.x = 'ChemID', by.y = 'Name')

#####
TRAIN <- train %>%
  mutate(set = 'train')
TEST <- test %>%
  mutate(set = 'test')
BCF <- rbind(TRAIN, TEST)
BCF_train_test <- ggplot(BCF, aes(LogBCF, stat(density), colour = set)) +
  geom_freqpoly(binwidth = 0.25, size = 1) +
  scale_color_manual(values = c('#EB6B4A', '#0B3087')) +
  theme(legend.position="none")
BCF_train_test

# alles <- rbind(train, test)
# 
# rm(TEST, test, testFP, TRAIN, train, trainFP)
# 
# include <- createDataPartition(alles$LogBCF, p = 0.8, list = FALSE, groups = 10)
# 
# train <- alles[include, ]
# X_train <- train %>%
#   select(-ChemID, LogBCF)
# y_train <- train %>%
#   select(LogBCF)
# test <- alles[-include, ]
# X_test <- test %>%
#   select(-ChemID, LogBCF)
# y_test <- test %>%
#   select(LogBCF)
# 
# summary(train$LogBCF)
# summary(test$LogBCF)

# curate data

## near-zero variance descriptors

nzv <- nearZeroVar(X_train, freqCut = 100/0)
X_train <- X_train[ , -nzv]
### and
X_test <- X_test[ , -nzv]

## highly correlated descriptors

correlations <- cor(X_train)
# corrplot::corrplot(correlations, order = 'hclust')
highCorr <- findCorrelation(correlations, cutoff = 0.85)
X_train <- X_train[ , -highCorr]
### and
X_test <- X_test[ , -highCorr]

# correlations <- cor(X_train)
# corrplot::corrplot(correlations, order = 'hclust')

## linear combinations

# comboInfo <- findLinearCombos(X_train) # returns NULL
# X_train <- X_train[ , -comboInfo$remove]
# ### and
# X_test <- X_test[ , -nzv]

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

BCF_PC <-
  ggplot(pcaPts, aes(PC1, PC2)) +
  geom_point(aes(colour = factor(dataset), shape = factor(dataset))) +
  labs(title = 'BCF') +
  theme(legend.position="none")
BCF_PC

# models

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  repeats = 5)

set.seed(350)

## multiple linear regression

trainSet <- cbind(y_train, X_trainTransformed)

mlr <- train(LogBCF ~ .,
             data = trainSet,
             method = 'lm',
             trControl = fitControl)

y_predict <- predict(mlr, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogBCF, data = data2plot))

BCF_mlr <-
  ggplot(data2plot, aes(LogBCF, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  # coord_equal() +
  # xlim(c(0, 3.5)) +
  # ylim(c(-1.0, 6.0)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Multiple Linear Regression') +
       # subtitle = 'Multiple Linear Regression\n test data') +
  ggthemes::theme_tufte()
BCF_mlr <- BCF_mlr + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
BCF_mlr

y_predict <- predict(mlr, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

mlrPR <- postResample(pred = y_predict, obs = X_trainTransformed)
rmse_train = c(mlrPR[1])
r2_train = c(mlrPR[2])

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogBCF, data = data2plot))

p <-
  ggplot(data2plot, aes(LogBCF, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogBCF',
       subtitle = 'Multiple Linear Regression\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

## partial least squares

plsModel <- train(
  LogBCF ~ .,
  data = trainSet,
  method = 'pls',
  tuneLength = 20,
  trControl = fitControl
)

y_predict <- predict(plsModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogBCF, data = data2plot))

BCF_pls <-
  ggplot(data2plot, aes(LogBCF, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  # coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Partial Least Squares') +
       # subtitle = 'Partial Least Squares\n test data') +
  ggthemes::theme_tufte()
BCF_pls <- BCF_pls + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
BCF_pls

y_predict <- predict(plsModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogBCF, data = data2plot))

p <-
  ggplot(data2plot, aes(LogBCF, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogBCF',
       subtitle = 'Partial Least Squares\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

## support vector machines

svmModel <- train(
  LogBCF ~ .,
  data = trainSet,
  method = 'svmRadial',
  # tuneLength = 14,
  trControl = fitControl
)

y_predict <- predict(svmModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogBCF, data = data2plot))

BCF_svm <-
  ggplot(data2plot, aes(LogBCF, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  # coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Support Vector Machines') +
       # subtitle = 'Support Vector Machines\n test data') +
  ggthemes::theme_tufte()
BCF_svm <- BCF_svm + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
BCF_svm

y_predict <- predict(svmModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogBCF, data = data2plot))

p <-
  ggplot(data2plot, aes(LogBCF, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogBCF',
       subtitle = 'Support Vector Machines\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

## k-nearest neighbors

knnModel <- train(
  LogBCF ~ .,
  data = trainSet,
  method = 'knn',
  tuneGrid = data.frame(.k = 1:20),
  trControl = fitControl
)

y_predict <- predict(knnModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogBCF, data = data2plot))

BCF_kNN <-
  ggplot(data2plot, aes(LogBCF, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  # coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'K-Nearest Neighbors') +
       # subtitle = 'K-Nearest Neighbors\n test data') +
  ggthemes::theme_tufte()
BCF_kNN <- BCF_kNN + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
BCF_kNN

y_predict <- predict(knnModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogBCF, data = data2plot))

p <-
  ggplot(data2plot, aes(LogBCF, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogBCF',
       subtitle = 'K-Nearest Neighbors\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

## boosted trees

gbmGrid <- expand.grid(
  .interaction.depth = seq(1, 7, by = 2),
  .n.trees = seq(100, 1000, by = 50),
  .shrinkage = c(0.001, 0.1),
  .n.minobsinnode = 3
)

treeModel <- train(
  LogBCF ~ .,
  data = trainSet,
  method = 'gbm',
  tuneGrid = gbmGrid,
  verbose = FALSE
)

y_predict <- predict(treeModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogBCF, data = data2plot))

BCF_gbm <-
  ggplot(data2plot, aes(LogBCF, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  # coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'Boosted Trees') +
       # subtitle = 'Boosted Trees\n test data') +
  ggthemes::theme_tufte()
BCF_gbm <- BCF_gbm + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
BCF_gbm

y_predict <- predict(treeModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogBCF, data = data2plot))

p <-
  ggplot(data2plot, aes(LogBCF, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogBCF',
       subtitle = 'Boosted Trees\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

alles_plot <-
  gridExtra::grid.arrange(BCF_train_test, BCF_mlr, BCF_pls, BCF_svm, BCF_kNN, BCF_gbm, nrow = 3)
ggsave('EcoTox/BCF_plots.png', alles_plot)
