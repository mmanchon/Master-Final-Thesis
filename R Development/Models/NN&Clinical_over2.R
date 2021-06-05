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

datos_exp <- read.csv(file = "ClinicalOutcomesDS/CO_atlas2018_over2.csv")


categorical_data <- datos_exp[,c("X","AGE", "SEX", "RACE", "ETHNICITY", "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT",
                                 "AJCC_PATHOLOGIC_TUMOR_STAGE", "WEIGHT", "Diagnosis.Age", "Mutation.Count",
                                 "Prior.Diagnosis", "Fraction.Genome.Altered")]

categorical_data$AgeRange <- ifelse(categorical_data$AGE < 57,"X57",ifelse(categorical_data$AGE > 73, "X73","57x73"))
categorical_data$DiagnosisRange <- ifelse(categorical_data$Diagnosis.Age < 57,"X57",ifelse(categorical_data$Diagnosis.Age > 73, "X73","57x73"))
categorical_data$MutationRange <- ifelse(categorical_data$Mutation.Count < 24,"X24",ifelse(categorical_data$Mutation.Count > 47, "X47","24x47"))
categorical_data$AlteredRange <- ifelse(categorical_data$Fraction.Genome.Altered < 0.001325,"X0.001325",ifelse(categorical_data$Fraction.Genome.Altered > 0.229775, "X0.229775","0.001325x0.229775"))
categorical_data$WEIGHT <- NULL

colnames(categorical_data)[which(names(categorical_data) == "AJCC_PATHOLOGIC_TUMOR_STAGE")] <- "STAGE"
colnames(categorical_data)[which(names(categorical_data) == "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT")] <- "RepeteadTumor"
colnames(datos_exp)[which(names(datos_exp) == "OS_STATUS")] <- "OsStatus"

categorical_data$STAGE <- as.character(categorical_data$STAGE)
categorical_data$STAGE[is.na(categorical_data$STAGE)] <- "Stage I"
categorical_data$STAGE <- as.factor(categorical_data$STAGE)

lista<-levels(categorical_data$STAGE)
lista_early<-lista[-c(6:7)]
lista_late<-lista[-c(1:5)]

levels(categorical_data$STAGE)[levels(categorical_data$STAGE) %in% lista_early]<-"Stage_I_II"
levels(categorical_data$STAGE)[levels(categorical_data$STAGE) %in% lista_late]<-"Stage_III_IV"

categorical_data <- categorical_data[,-which(colnames(categorical_data) %in% c("AGE","Diagnosis.Age","Mutation.Count","Fraction.Genome.Altered"))]

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
maxmindf <- as.data.frame(lapply(datos_exp[,c("EGLN2","C19orf26","POM121L8P", "MCM4", "MLXIPL", "ROBO4", "LGR4","TFRC")], normalize))
datos_exp[,c("EGLN2","C19orf26","POM121L8P", "MCM4", "MLXIPL", "ROBO4", "LGR4","TFRC")] <- maxmindf
data <- merge(x = categorical_data, y = datos_exp[,c("X","EGLN2","C19orf26","POM121L8P", "MCM4", "MLXIPL", "ROBO4", "LGR4","TFRC","OsStatus")], by = "X", all.x = TRUE)
data$STAGE <- as.factor(categorical_data$STAGE)
data[is.na(data)] <- -1
data$ETHNICITY <- as.character(data$ETHNICITY)
data$NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT <- as.character(data$NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT)
data[is.na(data)] <- ''

data_matrix <- model.matrix(~AgeRange+SEX+RACE+ETHNICITY+RepeteadTumor+STAGE+
                              DiagnosisRange+MutationRange+AlteredRange+
                              EGLN2+C19orf26+POM121L8P+MCM4+MLXIPL+ROBO4+LGR4+TFRC+OsStatus, data=data)
colnames(data_matrix)
colnames(data_matrix)[24] <- "OsStatusDeceased"
colnames(data_matrix)[5] <- "RACEBlackOrAfricanAmerican"
colnames(data_matrix)[7] <- "ETHNICITYHispanicOrLatino"
colnames(data_matrix)[8] <- "ETHNICITYNotHispanicOrLatino"
col_list <- paste(c(colnames(data_matrix[,-c(1,24)])),collapse="+")
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
