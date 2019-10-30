library(tidyverse)
library(magrittr)

# read data

## training data
train <-
  read.csv('cache/TR_RBioDeg_1197_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,-ROMol,-SMILES,-ID) %>%
  select(CAS, Ready_Biodeg, everything()) %>%
  na.omit()
train$Ready_Biodeg <- ifelse(train$Ready_Biodeg > 0.5, 'RB', 'NRB')

## test data
test <-
  read.csv('cache/TST_RBioDeg_411_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,-ROMol,-SMILES,-ID) %>%
  select(CAS, Ready_Biodeg, everything()) %>%
  na.omit()
test$Ready_Biodeg <- ifelse(test$Ready_Biodeg > 0.5, 'RB', 'NRB')

df <- rbind(train, test) %>%
  data.frame()

write.csv(df, 'ready_biodeg.csv', row.names = FALSE)
