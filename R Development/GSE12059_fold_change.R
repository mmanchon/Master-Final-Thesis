library(dplyr)

df <- read.csv("../data/geo/GSE12059_series_matrix.txt/GSE12059_expression.csv")

df$fold_1 <- df$GSM304508 / df$GSM304514
df$fold_2 <- df$GSM304509 / df$GSM304515
df$fold_3 <- df$GSM304510 / df$GSM304516
df$fold_4 <- df$GSM304511 / df$GSM304517
df$fold_5 <- df$GSM304512 / df$GSM304518
df$fold_6 <- df$GSM304513 / df$GSM304519

df$fold_increase <- rowMeans(df[c("fold_1","fold_2","fold_3","fold_4","fold_5","fold_6")])

df <- df[c("entrez","fold_increase")]

df <- df %>% arrange(desc(fold_increase)) %>% top_n(50)

write.csv(df,"../data/geo/fold_increase/GSE12059_fold_increase.csv", row.names = FALSE, col.names = FALSE)

