library(dplyr)
library(tidyverse)


args <- commandArgs(TRUE)

TOP_N <- strtoi(args[1])

df <- read.csv("../data/geo/GSE115406_series_matrix.txt/GSE115406_expression.csv")


df$fold_1 <- df$GSM3177685 / df$GSM3177682
df$fold_2 <- df$GSM3177686 / df$GSM3177683
df$fold_3 <- df$GSM3177687 / df$GSM3177684

df$fold_increase <- rowMeans(df[c("fold_1","fold_2","fold_3")])

df <- df[c("entrez","fold_increase")]

df <- df[!is.infinite(df$fold_increase),]

df <- df %>% arrange(desc(fold_increase)) 

df <- df[!duplicated(df$entrez),]

df <- df %>% top_n(TOP_N) %>% select("entrez")

write.table(df,"../data/geo/fold_increase/GSE115406_fold_increase_down.csv", row.names = FALSE,sep = ',', col.names = FALSE)
