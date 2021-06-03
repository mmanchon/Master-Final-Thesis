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

datos_exp <- read.csv(file = "ClinicalOutcomesDS/CO_atlas2018_over1.csv")


colnames(datos_exp)[which(names(datos_exp) == "OS_STATUS")] <- "OsStatus"

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
maxmindf <- as.data.frame(lapply(datos_exp[,c("ANO1" ,"C1orf86","CD44","CRYBA2","DCT","EIF4E","FLT4" ,"GNAO1","GULP1","ITGB1","NPAS3"  ,"PRLR","ROBO4" ,"SPAG6" ,"SRC"  ,"TFRC")], normalize))
datos_exp[,c("ANO1" ,"C1orf86","CD44","CRYBA2","DCT","EIF4E","FLT4" ,"GNAO1","GULP1","ITGB1","NPAS3"  ,"PRLR","ROBO4" ,"SPAG6" ,"SRC"  ,"TFRC")] <- maxmindf
data <- datos_exp[,c("ANO1" ,"C1orf86","CD44","CRYBA2","DCT","EIF4E","FLT4" ,"GNAO1","GULP1","ITGB1","NPAS3"  ,"PRLR","ROBO4" ,"SPAG6" ,"SRC"  ,"TFRC","OsStatus")]



data_matrix <- model.matrix(~ANO1+ C1orf86+CD44+CRYBA2+DCT + EIF4E + FLT4 +GNAO1+GULP1+ITGB1+NPAS3+PRLR+ROBO4
                              +SPAG6+SRC+TFRC+OsStatus, data=data)

colnames(data_matrix)
colnames(data_matrix)[18] <- "OsStatusDeceased"

col_list <- paste(c(colnames(data_matrix[,-c(1,18)])),collapse="+")
col_list <- paste(c("OsStatusDeceased~",col_list),collapse="")
f <- formula(col_list)

train <- as.data.frame(data_matrix[1:106,-c(1)])
test <- as.data.frame(data_matrix[107:173,-c(1)])

nmodel <- neuralnet(f,data=train,hidden=1,
                    threshold = 0.01,
                    learningrate.limit = NULL,
                    learningrate.factor =
                      list(minus = 0.5, plus = 1.2),
                    algorithm = "rprop+")


plot(nmodel)


temp_test <- test
temp_test$OsStatusDeceased <- NULL

nn.results <- compute(nmodel, temp_test)


results <- data.frame(actual = test$OsStatusDeceased, prediction = nn.results$net.result)


results


roundedresults<-sapply(results,round,digits=0)
roundedresultsdf=data.frame(roundedresults)



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

confMat<-confusionMatrix(as.factor(roundedresultsdf$prediction), as.factor(roundedresultsdf$actual))

draw_confusion_matrix(confMat)


nn_backprop <- neuralnet(f, data=train,
                         learningrate = 0.0001,
                         hidden=c(15,10,5),
                         algorithm = "backprop",
                         act.fct = "logistic",
                         linear.output = FALSE,
                         stepmax = 1000000000000)

plot(nn_backprop)

temp_test <- test
temp_test$OsStatusDeceased <- NULL

nn.results <- compute(nn_backprop, temp_test)


results <- data.frame(actual = test$OsStatusDeceased, prediction = nn.results$net.result)


results


roundedresults<-sapply(results,round,digits=0)
roundedresultsdf=data.frame(roundedresults)
confMat<-confusionMatrix(as.factor(roundedresultsdf$prediction), as.factor(roundedresultsdf$actual))

draw_confusion_matrix(confMat)
