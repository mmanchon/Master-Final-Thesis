library(dplyr)
library(tidyverse)

args <- commandArgs(TRUE)

TOP_N <- strtoi(args[1])

df <- read.csv("../data/geo/GSE12059_series_matrix.txt/GSE12059_expression.csv")

df$control <- rowMeans(df[c("GSM304508","GSM304510","GSM304512")])

df$mutated <- rowMeans(df[c("GSM304514","GSM304516", "GSM304518")])

df$fold_increase <- df$mutated / df$control

df <- df[c("entrez","fold_increase")]

df <- df %>% arrange(desc(fold_increase)) 

df <- df[!duplicated(df$entrez),]

df <- df %>% top_n(TOP_N) %>% select("entrez")

write.table(df,"../data/geo/fold_increase/GSE12059_fold_increase_over.csv", row.names = FALSE,sep = ",", col.names = FALSE)

