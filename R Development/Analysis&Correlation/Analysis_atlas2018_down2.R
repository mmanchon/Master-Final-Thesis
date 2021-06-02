library(utils)
library(openxlsx)
library(data.table)  # transpose
library(stringr)
library(corrplot)
library(dplyr)
library(randomForest)
library(caret)
library(fastDummies)
library(VIM)
library(corrplot)

datos_exp <- read.csv(file = "ClinicalOutcomesDS/CO_atlas2018_down2.csv")


# Analysis of missing values
nrow(datos_exp)
# Variables clínicas
colnames(datos_exp)

sapply(datos_exp[,c("AGE", "SEX", "RACE","ETHNICITY", "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT",
                    "AJCC_PATHOLOGIC_TUMOR_STAGE" ,"WEIGHT","OS_STATUS","OS_MONTHS" ,"DSS_STATUS",
                    "DSS_MONTHS", "DFS_STATUS" , "DFS_MONTHS", "PFS_STATUS", "PFS_MONTHS", "Diagnosis.Age",
                    "Mutation.Count", "Prior.Diagnosis", "Fraction.Genome.Altered", "DSS_MONTHS",
                    "DFS_STATUS","DFS_MONTHS","PFS_STATUS","PFS_MONTHS")], 
       function(x) sum(is.na(x)))

sapply(datos_exp[,c("OS_STATUS","OS_MONTHS","DSS_STATUS","DSS_MONTHS","DFS_STATUS","DFS_MONTHS","PFS_STATUS","PFS_MONTHS")], function(x) sum(is.na(x)))

aggr(datos_exp[,c("OS_STATUS","OS_MONTHS","DSS_STATUS","DSS_MONTHS","DFS_STATUS","DFS_MONTHS","PFS_STATUS","PFS_MONTHS")], col = c("skyblue", "red"),prop=FALSE, numbers=TRUE,
     cex.numbers=1, cex.axis=0.6, sortVars=TRUE,
     ylab = c("NA histogram", "Patron"), gap=2, oma = c(7,2,1,1))


table(datos_exp$OS_STATUS)
table(datos_exp$DSS_STATUS)
table(datos_exp$PFS_STATUS)

# Nulos en los genes
colnames(datos_exp)

sapply(datos_exp[,c("ARL4C" ,
                    "CDH15"  ,"TADA3" ,
                    "PDLIM5" ,"GAS2L1",
                    "OSBPL10","EGFR"  ,
                    "EFEMP1" ,"FRMD4B",
                    "ZC3H7B" ,"POFUT2",
                    "SASH1"  ,"DOCK9" ,
                    "KANK2"  ,"PCLO"  ,
                    "HFE"    ,"IL1A"  ,
                    "KCNMA1" ,"MOG"   ,
                    "PCDHGC3","PDE4D" ,
                    "NAT8B"  ,"ZNF44" ,
                    "PML"    ,"EPB41L4B" ,                              
                    "ELOVL2" ,"DOCK10",
                    "CHD7"   ,"PCDHGA3" ,                               
                    "SIM1"   ,"VDR"   ,
                    "CA12"   ,"PAX8"  ,
                    "FZD5"   ,"KCNAB1",
                    "ANTXR1" ,"LGR5"  ,
                    "RUNX1"  ,"TNFRSF25",                               
                    "NRP1"   ,"SKAP2" ,
                    "EIF1AY" ,"CD44")],
       function(x) sum(is.na(x)))

# Correlation analysis
datos_exp$SNORD116.22 <- NULL

temp<-datos_exp[,c("ARL4C" ,
                   "CDH15"  ,"TADA3" ,
                   "PDLIM5" ,"GAS2L1",
                   "OSBPL10","EGFR"  ,
                   "EFEMP1" ,"FRMD4B",
                   "ZC3H7B" ,"POFUT2",
                   "SASH1"  ,"DOCK9" ,
                   "KANK2"  ,"PCLO"  ,
                   "HFE"    ,"IL1A"  ,
                   "KCNMA1" ,"MOG"   ,
                   "PCDHGC3","PDE4D" ,
                   "NAT8B"  ,"ZNF44" ,
                   "PML"    ,"EPB41L4B" ,                              
                   "ELOVL2" ,"DOCK10",
                   "CHD7"   ,"PCDHGA3" ,                               
                   "SIM1"   ,"VDR"   ,
                   "CA12"   ,"PAX8"  ,
                   "FZD5"   ,"KCNAB1",
                   "ANTXR1" ,"LGR5"  ,
                   "RUNX1"  ,"TNFRSF25",                               
                   "NRP1"   ,"SKAP2" ,
                   "EIF1AY" ,"CD44",
                   "OS_STATUS","OS_MONTHS","DSS_MONTHS","PFS_MONTHS")]

vars<-c("ARL4C" ,
        "CDH15"  ,"TADA3" ,
        "PDLIM5" ,"GAS2L1",
        "OSBPL10","EGFR"  ,
        "EFEMP1" ,"FRMD4B",
        "ZC3H7B" ,"POFUT2",
        "SASH1"  ,"DOCK9" ,
        "KANK2"  ,"PCLO"  ,
        "HFE"    ,"IL1A"  ,
        "KCNMA1" ,"MOG"   ,
        "PCDHGC3","PDE4D" ,
        "NAT8B"  ,"ZNF44" ,
        "PML"    ,"EPB41L4B" ,                              
        "ELOVL2" ,"DOCK10",
        "CHD7"   ,"PCDHGA3" ,                               
        "SIM1"   ,"VDR"   ,
        "CA12"   ,"PAX8"  ,
        "FZD5"   ,"KCNAB1",
        "ANTXR1" ,"LGR5"  ,
        "RUNX1"  ,"TNFRSF25",                               
        "NRP1"   ,"SKAP2" ,
        "EIF1AY" ,"CD44","OS_STATUS")

vars_1<-c(
        "PCDHGC3","PDE4D" ,
        "NAT8B"  ,"ZNF44" ,
        "PML"    ,"EPB41L4B" ,                              
        "ELOVL2" ,"DOCK10",
        "CHD7"   ,"PCDHGA3" ,                               
        "SIM1"   ,"VDR"   ,
        "CA12"   ,"PAX8"  ,
        "FZD5"   ,"KCNAB1",
        "ANTXR1" ,"LGR5"  ,
        "RUNX1"  ,"TNFRSF25",                               
        "NRP1"   ,"SKAP2" ,
        "EIF1AY" ,"CD44","OS_STATUS")

vars_2<-c("ARL4C" ,
        "CDH15"  ,"TADA3" ,
        "PDLIM5" ,"GAS2L1",
        "OSBPL10","EGFR"  ,
        "EFEMP1" ,"FRMD4B",
        "ZC3H7B" ,"POFUT2",
        "SASH1"  ,"DOCK9" ,
        "KANK2"  ,"PCLO"  ,
        "HFE"    ,"IL1A"  ,
        "KCNMA1" ,"MOG"   ,
        "OS_STATUS")



# Convierto las variables de los genes a tipo numérico

temp[vars] <- sapply(temp[vars],as.numeric)
str(temp)

corrplot(cor(temp[,vars], method= "pearson"), method="color",type="upper",tl.col="black", tl.cex=0.6, tl.srt=45,addCoef.col="black",number.cex=0.7)
corrplot(cor(temp[,vars_1], method= "pearson"), method="color",type="upper",tl.col="black", tl.cex=0.6, tl.srt=45,addCoef.col="black",number.cex=0.7)
corrplot(cor(temp[,vars_2], method= "pearson"), method="color",type="upper",tl.col="black", tl.cex=0.6, tl.srt=45,addCoef.col="black",number.cex=0.7)

aux <- temp[, which(names(temp) %in% c("PDLIM5", "OSBPL10", "EGFR", "FRMD4B", "POFUT2", "HFE", "KCNMA1", "MOG", "ZNF44", "DOCK10", "VDR", "CA12", "RUNX1", "TNFRSF25", "SKAP2","OS_STATUS"))]

#aux <- aux[, -which(names(aux) %in% c("OS_MONTHS","DSS_MONTHS","PFS_MONTHS"))]

aux$OS_STATUS <- as.factor(aux$OS_STATUS)

# Hago el split entre train y test
train_index <- createDataPartition(aux$OS_STATUS, p = .6, 
                                   list = FALSE, 
                                   times = 1)

train <- aux[train_index,]
test  <- aux[-train_index,]
table(train$OS_STATUS)
table(test$OS_STATUS)
# Entreno el modelo
trControl <- trainControl(method = "repeatedcv",
                          number = 3,
                          repeats = 0)
rf<- train(OS_STATUS ~ .,
           data = train,
           method = "rf",
           metric = "Accuracy",
           #tuneGrid = tuneGrid,
           trControl = trControl,
           importance = TRUE,
           #nodesize = 5,
           #maxnodes = 8,
           #ntree = 300,
           #strata = train_dummy$Claudicacion, 
           #sampsize = c('1'=2381,'0'=2381),
           na.action=na.omit)
rf
# Hago la predicción sobre test3_dummy para confirmar las métricas de modelo
pred <-predict(rf, test, na.action=na.omit)
#pred <-predict(rf, test, na.action=na.omit)
# Extraigo la matriz de confusión
confMat<-confusionMatrix(pred, test$OS_STATUS, positive="1")
#confMat3<-confusionMatrix(pred3, test3$Claudicacion_umbral3,positive="1")
confMat

draw_confusion_matrix <- function(cm) {
        total <- sum(cm$table)
        res <- as.numeric(cm$table)
        # Generate color gradients. Palettes come from RColorBrewer.
        greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
        redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
        getColor <- function (greenOrRed = "green", amount = 0) {
                if (amount == 0)
                        return("#FFFFFF")
                palette <- greenPalette
                if (greenOrRed == "red")
                        palette <- redPalette
                colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
        }
        # set the basic layout
        layout(matrix(c(1,1,2)))
        par(mar=c(2,2,2,2))
        plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
        title('CONFUSION MATRIX', cex.main=2)
        # create the matrix 
        classes = colnames(cm$table)
        rect(150, 430, 240, 370, col=getColor("green", res[1]))
        text(195, 435, classes[1], cex=1.2)
        rect(250, 430, 340, 370, col=getColor("red", res[3]))
        text(295, 435, classes[2], cex=1.2)
        text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
        text(245, 450, 'Actual', cex=1.3, font=2)
        rect(150, 305, 240, 365, col=getColor("red", res[2]))
        rect(250, 305, 340, 365, col=getColor("green", res[4]))
        text(140, 400, classes[1], cex=1.2, srt=90)
        text(140, 335, classes[2], cex=1.2, srt=90)
        # add in the cm results
        text(195, 400, res[1], cex=1.6, font=2, col='white')
        text(195, 335, res[2], cex=1.6, font=2, col='white')
        text(295, 400, res[3], cex=1.6, font=2, col='white')
        text(295, 335, res[4], cex=1.6, font=2, col='white')
        # add in the specifics 
        plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
        text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
        text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
        text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
        text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
        text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
        text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
        text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
        text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
        text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
        text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
        # add in the accuracy information 
        text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
        text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
        text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
        text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}

draw_confusion_matrix(confMat)
# Importancia de las variables
importance <- varImp(rf,scale=TRUE)
print(importance)
plot(importance)

