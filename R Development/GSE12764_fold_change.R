library(dplyr)

df <- read.csv("../data/geo/GSE12764_series_matrix.txt/GSE12764_expression.csv")

df$fold_1 <- df$GSM320243 / df$GSM320249
df$fold_2 <- df$GSM320244 / df$GSM320250
df$fold_3 <- df$GSM320245 / df$GSM320251
df$fold_4 <- df$GSM320246 / df$GSM320252
df$fold_5 <- df$GSM320247 / df$GSM320253

df$fold_increase <- rowMeans(df[c("fold_1","fold_2","fold_3","fold_4","fold_5")])

df <- df[c("entrez","fold_increase")]

df <- df %>% arrange(desc(fold_increase)) %>% top_n(50)

write.table(df,"../data/geo/fold_increase/GSE12764_fold_increase.csv", row.names = FALSE,sep=",", col.names = FALSE)

