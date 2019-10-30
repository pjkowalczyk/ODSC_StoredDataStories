library(tidyverse)
library(magrittr)
library(caret)
library(ggplot2)
library(ggthemes)
library(jtools)
library(randomForest)
library(ISLR)

# read data
df <-
  read.csv('data/ready_biodeg.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  na.omit()

# view sample of data
head(df[sample(nrow(df), 10), ])

df0 <- df %>%
  group_by(Ready_Biodeg) %>%
  summarise(count = n())

ReadyBiodeg <- ggplot(df0, aes(x = Ready_Biodeg, y = count)) +
  theme_bw() +
  geom_bar(stat = "identity", fill = 'blue')
ReadyBiodeg

ggsave('graphics/RB_Ready_Biodeg_Histogram.jpg', plot = ReadyBiodeg)

# train / text split
inTrain <- caret::createDataPartition(df$Ready_Biodeg, p = 0.8, list = FALSE)
train <- df[inTrain, ]
test <- df[-inTrain, ]

X_train <- train[ , 3:ncol(train)]
y_train <- train[ , 2] %>% data.frame()
colnames(y_train) <- 'Ready_Biodeg'
X_test <- test[ , 3:ncol(test)]
y_test <- test[ , 2] %>% data.frame()
colnames(y_test) <- 'Ready_Biodeg'

TRAIN <- train %>%
  mutate(set = 'train')
TEST <- test %>%
  mutate(set = 'test')

RB <- rbind(TRAIN, TEST)

df1 <- RB %>%
  group_by(Ready_Biodeg, set) %>%
  summarise(count = n()) %>%
  data.frame()
df1$total <- ifelse(df1$set == 'train', nrow(train), nrow(test))
df1$pct <- df1$count / df1$total

p <- ggplot(df1, aes(x=Ready_Biodeg, y=pct, fill=set)) +
  scale_color_manual(values = c('#EB6B4A', '#0B3087')) +
  geom_bar(stat="identity", position=position_dodge())
p

ggsave('graphics/RB_Ready_Biodeg_Histogram_TrainTest.jpg', plot = p)

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
jpeg('graphics/RB_fullCorrelation.jpg')
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
jpeg('graphics/RB_reducedCorrelation.jpg')
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

RB_PC <-
  ggplot(pcaPts, aes(PC1, PC2, colour = factor(dataset))) +
  geom_point(aes(shape = factor(dataset))) +
  labs(title = 'Ready Biodeg PCA') +
  scale_color_manual(values = c('#EB6B4A', '#0B3087')) +
  theme(legend.position="none") +
  ggthemes::theme_tufte()
RB_PC

ggsave('graphics/RB_Ready_Biodeg_PC.jpg', plot = RB_PC)

# models

fitControl <- trainControl(## 5-fold CV
  method = "repeatedcv",
  repeats = 5)

set.seed(42)

trainSet <- cbind(y_train, X_trainTransformed)
trestSet <- cbind(y_test, X_testTransformed)

## Logistic Regression

lrModel <- glm(Ready_Biodeg ~ .,
          data = trainSet,
          family = binomial(link="logit"))

y_predict <- predict(lrModel, newdata = X_testTransformed, type = 'response') %>%
  data.frame()
colnames(y_predict) <- c('Probability')
y_predict$Predicted <- as_factor(ifelse(y_predict$Probability >= 0.5, 'RB', 'NRB'))

data2plot <- cbind(y_test, y_predict)

confusionMatrix(data2plot$Predicted, data2plot$Ready_Biodeg, positive = 'RB', dnn = c("Prediction", "Reference"))

cm <- confusionMatrix(data2plot$Predicted, data2plot$Ready_Biodeg, positive = 'RB', dnn = c("Prediction", "Reference"))
cm
fourfoldplot(cm$table)

y_predict <- predict(lrModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Probability')
y_predict$Predicted <- as_factor(ifelse(y_predict$Probability >= 0.5, 'RB', 'NRB'))

data2plot <- cbind(y_train, y_predict)

confusionMatrix(data2plot$Predicted, data2plot$Ready_Biodeg, positive = 'RB', dnn = c("Prediction", "Reference"))

cm <- confusionMatrix(data2plot$Predicted, data2plot$Ready_Biodeg, positive = 'RB', dnn = c("Prediction", "Reference"))
cm
fourfoldplot(cm$table)

## k-nearest neighbors

knnModel <- train(
  Ready_Biodeg ~ .,
  data = trainSet,
  method = 'knn',
  tuneGrid = data.frame(.k = 1:20),
  trControl = fitControl
)

y_predict <- predict(knnModel, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

confusionMatrix(data2plot$Predicted, data2plot$Ready_Biodeg, positive = 'RB', dnn = c("Prediction", "Reference"))

cm <- confusionMatrix(data2plot$Predicted, data2plot$Ready_Biodeg, positive = 'RB', dnn = c("Prediction", "Reference"))
cm
fourfoldplot(cm$table)

y_predict <- predict(knnModel, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

confusionMatrix(data2plot$Predicted, data2plot$Ready_Biodeg, positive = 'RB', dnn = c("Prediction", "Reference"))

cm <- confusionMatrix(data2plot$Predicted, data2plot$Ready_Biodeg, positive = 'RB', dnn = c("Prediction", "Reference"))
cm
fourfoldplot(cm$table)

## Random Forests

rf <- randomForest(Ready_Biodeg ~ ., data = trainSet, ntree = 250, nodesize = 7)

y_predict <- predict(rf, newdata = X_testTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_test, y_predict)

confusionMatrix(data2plot$Predicted, data2plot$Ready_Biodeg, positive = 'RB', dnn = c("Prediction", "Reference"))

cm <- confusionMatrix(data2plot$Predicted, data2plot$Ready_Biodeg, positive = 'RB', dnn = c("Prediction", "Reference"))
cm
fourfoldplot(cm$table)

y_predict <- predict(rf, newdata = X_trainTransformed) %>%
  data.frame()
colnames(y_predict) <- c('Predicted')

data2plot <- cbind(y_train, y_predict)

cm <- confusionMatrix(data2plot$Predicted, data2plot$Ready_Biodeg, positive = 'RB', dnn = c("Prediction", "Reference"))
cm
fourfoldplot(cm$table)
