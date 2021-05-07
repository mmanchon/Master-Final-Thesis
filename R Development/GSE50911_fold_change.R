library(dplyr)
library(tidyverse)

args <- commandArgs(TRUE)

TOP_N <- strtoi(args[1])

df <- read.csv("../data/geo/GSE50911_series_matrix.txt/GSE50911_expression.csv")
df <- unique(df)

df$fold_increase <- df$GSM1232317 / df$GSM1232315

df <- df[c("entrez","fold_increase")]

df <- df %>% arrange(desc(fold_increase)) 

df <- df[!duplicated(df$entrez),]

df <- df %>% top_n(TOP_N) %>% select("entrez")

write.table(df,"../data/geo/fold_increase/GSE50911_fold_increase.csv", row.names = FALSE,sep = ',', col.names = FALSE)
