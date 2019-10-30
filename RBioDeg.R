library(rcdk)
library(tidyverse)
library(magrittr)
library(purrr)
library(stringr)
library(caret)
library(corrplot)
library(ggplot2)
library(ggthemes)
library(pROC)
library(egg)

# read data

## training data
train <-
  read.csv('cache/TR_RBioDeg_1197_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,-CAS,-ROMol,-SMILES,-ID) %>%
  select(Ready_Biodeg, everything()) %>%
  na.omit()
train$Ready_Biodeg <- ifelse(train$Ready_Biodeg > 0.5, 'RB', 'NRB')

X_train <- train %>%
  select(-Ready_Biodeg)
y_train <- train %>%
  select(Ready_Biodeg) %>%
  data.frame()

## test data
test <-
  read.csv('cache/TST_RBioDeg_411_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,-CAS,-ROMol,-SMILES,-ID) %>%
  select(Ready_Biodeg, everything()) %>%
  na.omit()
test$Ready_Biodeg <- ifelse(test$Ready_Biodeg > 0.5, 'RB', 'NRB')

X_test <- test %>%
  select(-Ready_Biodeg)
y_test <- test %>%
  select(Ready_Biodeg) %>%
  data.frame()

# curate data

## near-zero variance descriptors

nzv <- nearZeroVar(X_train, freqCut = 100/0)
X_train <- X_train[ , -nzv]
### and
X_test <- X_test[ , -nzv]

## highly correlated descriptors

correlations <- cor(X_train)
alles_plot <- corrplot::corrplot(correlations, order = 'hclust', tl.cex = 0.6)
highCorr <- findCorrelation(correlations, cutoff = 0.85)
X_train <- X_train[ , -highCorr]
### and
X_test <- X_test[ , -highCorr]

correlations <- cor(X_train)
noCorr_plot <- corrplot::corrplot(correlations, order = 'hclust', tl.cex = 0.8)

## linear combinations

comboInfo <- findLinearCombos(X_train) # returns NULL
X_train <- X_train[ , -comboInfo$remove]
### and
X_test <- X_test[ , -comboInfo$remove]

## center & scale descriptors

preProcValues <- preProcess(X_train, method = c("center", "scale"))

X_trainTransformed <- predict(preProcValues, X_train)
### and
X_testTransformed <- predict(preProcValues, X_test)

# models

## support vector machines

set.seed(350)

ctrl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10,
  classProbs = TRUE)

library(kernlab)
sigmaRangeReduced <- sigest(as.matrix(X_trainTransformed))

svmGridReduced <- expand.grid(.sigma = sigmaRangeReduced[1],
                              .C = 2^(seq(-4, 4)))

trainSet <- cbind(y_train, X_trainTransformed)

svmModel <- train(Ready_Biodeg ~ .,
                  data = trainSet,
                  method = 'svmRadial',
                  metric = 'ROC',
                  tuneGrid = svmGridReduced,
                  fit = FALSE,
                  trControl = ctrl)

y_predict <- predict(svmModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

confusionMatrix(data = as.factor(y_predict$Predicted),
                reference = as.factor(y_test$Ready_Biodeg),
                positive = 'NRB')

library(pROC)

y_test$endpt <- ifelse(y_test$Ready_Biodeg == 'NRB', 0, 1)
y_predict$endpt <- ifelse(y_predict$Predicted == 'NRB', 0, 1)
rocCurve <- roc(response = y_test$endpt,
                predictor = y_predict$endpt)
auc(rocCurve)
plot(rocCurve, legacy.axes = TRUE)

#####...#####...#####...#####

library(randomForest)
library(caret)
library(e1071)

# read data

## training data
train <-
  read.csv('cache/TR_RBioDeg_1197_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,-CAS,-ROMol,-SMILES,-ID) %>%
  select(Ready_Biodeg, everything()) %>%
  na.omit()
train$Ready_Biodeg <- ifelse(train$Ready_Biodeg > 0.5, 'RB', 'NRB')
## test data
test <-
  read.csv('cache/TST_RBioDeg_411_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,-CAS,-ROMol,-SMILES,-ID) %>%
  select(Ready_Biodeg, everything()) %>%
  na.omit()
test$Ready_Biodeg <- ifelse(test$Ready_Biodeg > 0.5, 'RB', 'NRB')

## bind train and test, by row
alles <- rbind(train, test)

## data splitting
set.seed(350)
trainIndex <- createDataPartition(alles$Ready_Biodeg, p = .8, 
                                  list = FALSE, 
                                  times = 1)

train <- alles[trainIndex, ]
test <- alles[-trainIndex, ]

X_train <- train %>%
  select(-Ready_Biodeg)
y_train <- train %>%
  select(Ready_Biodeg) %>%
  data.frame() %>%
  mutate(Ready_Biodeg = as.factor(Ready_Biodeg))

X_test <- test %>%
  select(-Ready_Biodeg)
y_test <- test %>%
  select(Ready_Biodeg) %>%
  data.frame() %>%
  mutate(Ready_Biodeg = as.factor(Ready_Biodeg))

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
X_train <- X_train[ , -comboInfo$remove]
### and
X_test <- X_test[ , -comboInfo$remove]

# 10 fold; repeat 3 times
control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3)

# metric: Accuracy
metric <- "Accuracy"

mtry <- sqrt(ncol(X_train))

tunegrid <- expand.grid(.mtry=mtry)

data2model <- cbind(y_train, X_train)

rf_default <- train(
  Ready_Biodeg ~ .,
  data = data2model,
  method = 'rf',
  metric = 'Accuracy',
  tuneGrid = tunegrid,
  trControl = control
)

print(rf_default)

y_predict <- predict(rf_default, newdata = X_test) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

confusionMatrix(data = as.factor(y_predict$Predicted),
                reference = as.factor(y_test$Ready_Biodeg),
                positive = 'NRB')

library(doParallel)
cores <- 3
registerDoParallel(cores = cores)

mtry <- sqrt(ncol(X_train))

#ntree: Number of trees to grow.
ntree <- 3

control <- trainControl(
  method = 'repeatedcv',
  number = 10,
  repeats = 3,
  search = 'random'
)

#
rf_random <- train(Ready_Biodeg ~ .,
                   data = data2model,
                   method = 'rf',
                   metric = 'Accuracy',
                   tuneLength  = 15, 
                   trControl = control)

print(rf_random)

plot(rf_random)

#
control <- trainControl(method='repeatedcv', 
                        number=10, 
                        repeats=3, 
                        search='grid')
#
tunegrid <- expand.grid(.mtry = (1:15)) 

rf_gridsearch <- train(Ready_Biodeg ~ ., 
                       data = data2model,
                       method = 'rf',
                       metric = 'Accuracy',
                       tuneGrid = tunegrid)

 