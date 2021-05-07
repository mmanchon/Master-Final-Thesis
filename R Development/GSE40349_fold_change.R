# library(dplyr)
# 
# df <- read.csv("../data/geo/GSE40349_series_matrix.txt/GSE40349_expression.csv")
# 
# df$fold_1 <- df$GSM991787 / df$GSM991789
# df$fold_1 <- df$GSM991787 / df$GSM991789
# df$fold_1 <- df$GSM991787 / df$GSM991789
# df$fold_1 <- df$GSM991787 / df$GSM991789
# 
# df$fold_1 <- df$GSM991787 / df$GSM991789
# df$fold_1 <- df$GSM991787 / df$GSM991789
# df$fold_1 <- df$GSM991787 / df$GSM991789
# df$fold_1 <- df$GSM991787 / df$GSM991789
# 
# 
# df$fold_increase <- rowMeans(df[c("fold_1","fold_2","fold_3","fold_4","fold_5","fold_6")])
# 
# df <- df[c("entrez","fold_increase")]
# 
# df <- df %>% arrange(desc(fold_increase)) %>% top_n(50)
# 
# write.csv(df,"../data/geo/fold_increase/GSE40349_fold_increase.csv", row.names = FALSE)
# 
