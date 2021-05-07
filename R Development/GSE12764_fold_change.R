library(dplyr)
library(tidyverse)


args <- commandArgs(TRUE)

TOP_N <- strtoi(args[1])

df <- read.csv("../data/geo/GSE12764_series_matrix.txt/GSE12764_expression.csv")

df <- unique(df)

df$fold_1 <- df$GSM320249 / df$GSM320243 
df$fold_2 <- df$GSM320250 / df$GSM320244 
df$fold_3 <- df$GSM320251 / df$GSM320245 
df$fold_4 <- df$GSM320252 / df$GSM320246 
df$fold_5 <- df$GSM320253 / df$GSM320247 

df$fold_increase <- rowMeans(df[c("fold_1","fold_2","fold_3","fold_4","fold_5")])

df <- df[c("entrez","fold_increase")]

df <- df %>% arrange(desc(fold_increase)) 

df <- df[!duplicated(df$entrez),]

df <- df %>% top_n(TOP_N) %>% select("entrez")

write.table(df,"../data/geo/fold_increase/GSE12764_fold_increase.csv", row.names = FALSE,sep=",", col.names = FALSE)

