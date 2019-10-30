library(rcdk)
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
  read.csv('cache/TR_AOH_516_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,-CAS,-ROMol,-SMILES,-ID) %>%
  select(LogOH, everything()) %>%
  na.omit()

X_train <- train %>%
  select(-LogOH)
y_train <- train %>%
  select(LogOH) %>%
  data.frame()

## test data
test <-
  read.csv('cache/TST_AOH_176_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,-CAS,-ROMol,-SMILES,-ID) %>%
  select(LogOH, everything()) %>%
  na.omit()

X_test <- test %>%
  select(-LogOH)
y_test <- test %>%
  select(LogOH) %>%
  data.frame()

# curate data

## near-zero variance descriptors

nzv <- nearZeroVar(X_train, freqCut = 100/0)
X_train <- X_train[ , -nzv]
### and
X_test <- X_test[ , -nzv]

## highly correlated descriptors

correlations <- cor(X_train)
corrplot::corrplot(correlations, order = 'hclust')
highCorr <- findCorrelation(correlations, cutoff = 0.85)
X_train <- X_train[ , -highCorr]
### and
X_test <- X_test[ , -highCorr]

correlations <- cor(X_train)
corrplot::corrplot(correlations, order = 'hclust')

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

mlr <- train(LogOH ~ .,
             data = trainSet,
             method = 'lm',
             trControl = fitControl)

y_predict <- predict(mlr, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogOH, data = data2plot))

p <-
  ggplot(data2plot, aes(LogOH, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'LogOH',
       subtitle = 'Multiple Linear Regression\n test data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

y_predict <- predict(mlr, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

mlrPR <- postResample(pred = y_predict, obs = X_trainTransformed)
rmse_train = c(mlrPR[1])
r2_train = c(mlrPR[2])

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogOH, data = data2plot))

p <-
  ggplot(data2plot, aes(LogOH, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogOH',
       subtitle = 'Multiple Linear Regression\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

## partial least squares

plsModel <- train(
  LogOH ~ .,
  data = trainSet,
  method = 'pls',
  tuneLength = 20,
  trControl = fitControl
)

y_predict <- predict(plsModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogOH, data = data2plot))

p <-
  ggplot(data2plot, aes(LogOH, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'LogOH',
       subtitle = 'Partial Least Squares\n test data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

y_predict <- predict(plsModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogOH, data = data2plot))

p <-
  ggplot(data2plot, aes(LogOH, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogOH',
       subtitle = 'Partial Least Squares\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

## support vector machines

svmModel <- train(
  LogOH ~ .,
  data = trainSet,
  method = 'svmRadial',
  # tuneLength = 14,
  trControl = fitControl
)

y_predict <- predict(svmModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogOH, data = data2plot))

p <-
  ggplot(data2plot, aes(LogOH, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'LogOH',
       subtitle = 'Support Vector Machines\n test data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

y_predict <- predict(svmModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogOH, data = data2plot))

p <-
  ggplot(data2plot, aes(LogOH, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogOH',
       subtitle = 'Support Vector Machines\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

## k-nearest neighbors

knnModel <- train(
  LogOH ~ .,
  data = trainSet,
  method = 'knn',
  tuneGrid = data.frame(.k = 1:20),
  trControl = fitControl
)

y_predict <- predict(knnModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogOH, data = data2plot))

p <-
  ggplot(data2plot, aes(LogOH, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'LogOH',
       subtitle = 'K-Nearest Neighbors\n test data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

y_predict <- predict(knnModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogOH, data = data2plot))

p <-
  ggplot(data2plot, aes(LogOH, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogOH',
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
  LogOH ~ .,
  data = trainSet,
  method = 'gbm',
  tuneGrid = gbmGrid,
  verbose = FALSE
)

y_predict <- predict(treeModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

summary(lm(Predicted ~ LogOH, data = data2plot))

p <-
  ggplot(data2plot, aes(LogOH, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method = 'lm') +
  labs(title = 'LogOH',
       subtitle = 'Boosted Trees\n test data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

y_predict <- predict(treeModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

summary(lm(Predicted ~ LogOH, data = data2plot))

p <-
  ggplot(data2plot, aes(LogOH, Predicted)) +
  geom_point(colour = "blue", size = 2) +
  coord_equal() +
  # xlim(c(0, 3.5)) + ylim(c(0, 3.5)) +
  geom_smooth(method='lm') +
  labs(title = 'LogOH',
       subtitle = 'Boosted Trees\n training data') +
  ggthemes::theme_tufte()
p <- p + geom_abline(intercept = 0,
                     slope = 1,
                     colour = 'red')
p

#####
##### 2018-10-08
#####
##### repeat the model building process, with descriptors calculated
##### using rcdk
#####

## rcdk descriptors
dc <- get.desc.categories()
dc
## hybrid descriptors
dn01 <- get.desc.names(dc[1])
dn01
## constitutional descriptors
dn02 <- get.desc.names(dc[2])
dn02
## topological descriptors
dn03 <- get.desc.names(dc[3])
dn03
## electronic descriptors
dn04 <- get.desc.names(dc[4])
dn04
## geometrical descriptors
dn05 <- get.desc.names(dc[5])
dn05

# read data

## training data
train <-
  read.csv('cache/TR_AOH_516_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(CAS,SMILES,LogOH) %>%
  na.omit()

mol <- parse.smiles(train$SMILES)

allDescs01 <- eval.desc(mol, dn01)
allDescs02 <- eval.desc(mol, dn02)
allDescs03 <- eval.desc(mol, dn03)
allDescs04 <- eval.desc(mol, dn04) %>%
  select(-tpsaEfficiency, -TopoPSA)
allDescs05 <- eval.desc(mol, dn05) %>%
  select(-topoShape, -geomShape, -PPSA.1, -PPSA.2, -PPSA.3, -PNSA.1, -PNSA.2, -PNSA.3) %>%
  select(-DPSA.1, -DPSA.2, -DPSA.3, -FPSA.1, -FPSA.2, -FPSA.3, -FNSA.1, -FNSA.2, -FNSA.3) %>%
  select(-WPSA.1, -WPSA.2, -WPSA.3, -WNSA.1, -WNSA.2, -WNSA.3) %>%
  select(-RPCG, -RNCG, -RPCS, -RNCS, -THSA, -TPSA, -RHSA, -RPSA)

allDescs <- cbind(allDescs01, allDescs02, allDescs03, allDescs04, allDescs05)

train <- cbind(train, allDescs)

train <- train %>%
  select(LogOH, everything())

X_train <- train %>%
  select(-LogOH, -CAS, -SMILES)
y_train <- train %>%
  select(LogOH) %>%
  data.frame()

## test data
test <-
  read.csv('cache/TST_AOH_176_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(CAS,SMILES,LogOH) %>%
  na.omit()

mol <- parse.smiles(test$SMILES)

allDescs01 <- eval.desc(mol, dn01)
allDescs02 <- eval.desc(mol, dn02)
allDescs03 <- eval.desc(mol, dn03)
allDescs04 <- eval.desc(mol, dn04) %>%
  select(-tpsaEfficiency, -TopoPSA)
allDescs05 <- eval.desc(mol, dn05) %>%
  select(-topoShape, -geomShape, -PPSA.1, -PPSA.2, -PPSA.3, -PNSA.1, -PNSA.2, -PNSA.3) %>%
  select(-DPSA.1, -DPSA.2, -DPSA.3, -FPSA.1, -FPSA.2, -FPSA.3, -FNSA.1, -FNSA.2, -FNSA.3) %>%
  select(-WPSA.1, -WPSA.2, -WPSA.3, -WNSA.1, -WNSA.2, -WNSA.3) %>%
  select(-RPCG, -RNCG, -RPCS, -RNCS, -THSA, -TPSA, -RHSA, -RPSA)

allDescs <- cbind(allDescs01, allDescs02, allDescs03, allDescs04, allDescs05)

test <- cbind(test, allDescs)

test <- test %>%
  select(LogOH, everything())

X_test <- test %>%
  select(-LogOH, -CAS, -SMILES)
y_test <- test %>%
  select(LogOH) %>%
  data.frame()

# curate data

## near-zero variance descriptors

nzv <- nearZeroVar(X_train, freqCut = 100/0)
X_train <- X_train[ , -nzv] %>%
  na.omit()
### and
X_test <- X_test[ , -nzv] %>%
  na.omit()

## highly correlated descriptors

correlations <- cor(X_train)
corrplot::corrplot(correlations, order = 'hclust')
highCorr <- findCorrelation(correlations, cutoff = 0.85)
X_train <- X_train[ , -highCorr]
### and
X_test <- X_test[ , -highCorr]

correlations <- cor(X_train)
corrplot::corrplot(correlations, order = 'hclust')

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


newdataset1<-dataset1[-which(is.na(dataset1$x)),]