library(dplyr)
library(tidyverse)


args <- commandArgs(TRUE)

TOP_N <- strtoi(args[1])

df <- read.csv("../data/geo/GSE55942_series_matrix.txt/GSE55942_expression.csv")


df$fold_1 <- df$GSM1348905 / df$GSM1348907
df$fold_2 <- df$GSM1348906 / df$GSM1348908

df$fold_increase <- rowMeans(df[c("fold_1","fold_2")])

df <- df[c("entrez","fold_increase")]

df <- df[!is.infinite(df$fold_increase),]

df <- df %>% arrange(desc(fold_increase)) 

df <- df[!duplicated(df$entrez),]

df <- df %>% top_n(TOP_N) %>% select("entrez")

write.table(df,"../data/geo/fold_increase/GSE55942_fold_increase_down.csv", row.names = FALSE,sep = ',', col.names = FALSE)
