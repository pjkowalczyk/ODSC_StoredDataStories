# load libraries
library(dplyr)
library(magrittr)
library(caret)
library(ggplot2)

# retrieve data
df <- read.csv('data/EGFR_data.csv')

# review column names
names(df)

#####
# inStudy <-  createDataPartition(df$Target, p = 0.2, list = FALSE, groups = min(5, nrow(df)))
# df <- df[inStudy, ]
#####

# Histogram overlaid with kernel density curve
ggplot(df, aes(x = Target)) +
  geom_histogram(
    aes(y = ..density..),
    # Histogram with density instead of count on y-axis
    binwidth = 0.25,
    colour = "black",
    fill = "white"
  ) +
  geom_density(alpha = .2, fill = "#FF6666")  # Overlay with transparent density plot

X <- df %>%
  dplyr::select(-c('Identifier', 'Target')) %>%
  data.frame()
y <- df %>%
  dplyr::select(c('Target')) %>%
  data.frame()

inTrain = createDataPartition(df$Target, p = 0.8, list = FALSE, groups = min(5, length(y)))
X_train <- X[inTrain, ]
X_test <- X[-inTrain, ]
y_train <- y[inTrain]
y_test <- y[-inTrain]

y01 <- y_train %>%
  data.frame() %>%
  rename('Target' = '.') %>%
  mutate('tranche' = 'train')
y02 <- y_test %>%
  data.frame() %>%
  rename('Target' = '.') %>%
  mutate('tranche' = 'test')
y03 <- rbind(y01, y02)
# Density plots with semi-transparent fill
ggplot(y03, aes(x = Target, fill = tranche)) + geom_density(alpha = .3)

# a bit of housekeeping
rm(y01, y02, y03)

nzv <- nearZeroVar(X_train, freqCut = 95/5)
filtered_X_train <- X_train[ , -nzv]
filtered_X_test <- X_test[ , -nzv] 

descrCor <- cor(filtered_X_train)
highCorr <- findCorrelation(descrCor, cutoff = 0.9)
filtered_X_train <- filtered_X_train[ , -highCorr]
filtered_X_test <- filtered_X_test[ , -highCorr]

comboInfo <- findLinearCombos(filtered_X_train)
# comboInfo
if( !is.null(comboInfo$remove) ) {
  filtered_X_train <- filtered_X_train[ , -comboInfo$remove]
  filtered_X_test <- filtered_X_test[ , -comboInfo$remove]
}

training <- filtered_X_train
test <- filtered_X_test

preProcValues <- preProcess(training, method = c("center", "scale"))

X_TRAIN <- predict(preProcValues, training)
X_TEST <- predict(preProcValues, test)

TRAIN <- cbind(y_train, X_TRAIN) %>%
  rename(Target = y_train)
TEST <- cbind(y_test, X_TEST) %>%
  rename(Target = y_test)

save(df, TRAIN, TEST, X_TRAIN, X_TEST, y_train, y_test, file = "data/data.RData")

load("data/data.RData")
