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
library(neuralnet)
library(mxnet)

datos_exp <- read.csv(file = "dataset.csv")


# Analysis of missing values



datos_exp <- datos_exp[,c("ANO1" ,"C1orf86","CD44","CRYBA2","DCT","EIF4E","FLT4" ,"GNAO1","GULP1","ITGB1","NPAS3"  ,"PRLR","ROBO4" ,"SPAG6" ,"SRC"  ,"TFRC","OS_STATUS")]

datos_exp$OS_STATUS<-as.character(datos_exp$OS_STATUS)
datos_exp$OS_STATUS[which(datos_exp$OS_STATUS=="0:LIVING")]<-"LIVING"
datos_exp$OS_STATUS[which(datos_exp$OS_STATUS=='1:DECEASED')]<-"DECEASED"

# Hago el split entre train y test
train_index <- createDataPartition(datos_exp$OS_STATUS, p = .6, 
                                   list = FALSE, 
                                   times = 1)

train <- datos_exp[train_index,]
test  <- datos_exp[-train_index,]

nn=neuralnet(OS_STATUS ~ .,data=train, hidden=c(5,3),
             linear.output = FALSE, threshold=0.01)

plot(nn)

temp_test <- subset(test, select = c("ANO1" ,"C1orf86","CD44","CRYBA2","DCT","EIF4E","FLT4" ,"GNAO1","GULP1","ITGB1","NPAS3"  ,"PRLR","ROBO4" ,"SPAG6" ,"SRC"  ,"TFRC"))

Predict=compute(nn,temp_test)
Predict$net.result

prob <- Predict$net.result
pred <- ifelse(prob>0.5, 1, 0)
pred


predicts=data.frame("OS_STATUS"=ifelse(max.col(Predict$net.result[ ,1:2])==1, "LIVING",
                              ifelse(max.col(Predict$net.result[ ,1:2])==2, "DECEASED", "NA")))

test$OS_STATUS <- as.character(test$OS_STATUS)
test$OS_STATUS[which(test$OS_STATUS=="0:LIVING")]<-"LIVING"
test$OS_STATUS[which(test$OS_STATUS=='1:DECEASED')]<-"DECEASED"

cm = confusionMatrix(as.factor(test$OS_STATUS), predicts$OS_STATUS)

print(cm)

results <- data.frame(actual = as.numeric(test$OS_STATUS), prediction = Predict$net.result)

roundedresults<-sapply(results,round,digits=0)
roundedresultsdf=data.frame(roundedresults)
attach(roundedresultsdf)
table(actual,prediction)
