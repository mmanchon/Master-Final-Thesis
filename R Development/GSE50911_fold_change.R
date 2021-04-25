df <- read.csv("../data/geo/GSE50911_series_matrix.txt/GSE50911_expression.csv")

df2 <- read.table(file = '../data/geo/GSE50911_series_matrix.txt/GSE50911.top.table.tsv', sep = '\t', header = TRUE)

df$control_2 <- rowMeans(df[c("GSM1232298","GSM1232300")])

df$mutated_2 <- rowMeans(df[c("GSM1232302","GSM1232305", "GSM1232307", "GSM1232308", "GSM1232309")])

df$control_4 <- rowMeans(df[c("GSM1232314","GSM1232315")])

df$mutated_4 <- rowMeans(df[c("GSM1232316", "GSM1232317", "GSM1232318", "GSM1232319")])

df <- df[c("entrez","control_2","mutated_2","control_4","mutated_4")]

df$fold_2 <- df$mutated_2 / df$control_2

df$fold_4 <- df$mutated_4 /df$control_4

df$fold_increase <- rowMeans(df[c("fold_2","fold_4")])

df <- df[c("entrez","fold_increase")]

df <- df %>% arrange(desc(fold_increase)) %>% top_n(50)

write.table(df,"../data/geo/fold_increase/GSE50911_fold_increase.csv", row.names = FALSE,sep = ',', col.names = FALSE)
