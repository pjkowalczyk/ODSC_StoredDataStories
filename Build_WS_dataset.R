library(tidyverse)
library(magrittr)

# read data

## training data
train <-
  read.csv('cache/TR_WS_3158_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,CAS,-ROMol,-SMILES,-ID) %>%
  select(CAS, LogMolar, everything()) %>%
  na.omit()

## test data
test <-
  read.csv('cache/TST_WS_1066_descrs.csv',
           header = TRUE,
           stringsAsFactors = FALSE) %>%
  select(-X,CAS,-ROMol,-SMILES,-ID) %>%
  select(CAS, LogMolar, everything()) %>%
  na.omit()

ws <- rbind(train, test)

write.csv(ws, file = "data/water_solubility.csv", row.names = FALSE)