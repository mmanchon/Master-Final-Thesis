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

categorical_data <- datos_exp[,c("X","AGE", "SEX", "RACE", "ETHNICITY", "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT",
                              "AJCC_PATHOLOGIC_TUMOR_STAGE", "WEIGHT", "Diagnosis.Age", "Mutation.Count",
                              "Prior.Diagnosis", "Fraction.Genome.Altered")]

categorical_data$AGE_RANGE <- ifelse(categorical_data$AGE < 57,"X<57",ifelse(categorical_data$AGE > 73, "X>73","57<=x<=73"))
categorical_data$DIAGNOSIS_RANGE <- ifelse(categorical_data$Diagnosis.Age < 57,"X<57",ifelse(categorical_data$Diagnosis.Age > 73, "X>73","57<=x<=73"))
categorical_data$MUTATION_RANGE <- ifelse(categorical_data$Mutation.Count < 24,"X<24",ifelse(categorical_data$Mutation.Count > 47, "X>47","24<=x<=47"))
categorical_data$ALTERED_RANGE <- ifelse(categorical_data$Fraction.Genome.Altered < 0.001325,"X<0.001325",ifelse(categorical_data$Fraction.Genome.Altered > 0.229775, "X>0.229775","0.001325<=x<=0.229775"))
categorical_data$WEIGHT <- NULL

colnames(categorical_data)[which(names(categorical_data) == "AJCC_PATHOLOGIC_TUMOR_STAGE")] <- "STAGE"

categorical_data$STAGE <- as.character(categorical_data$STAGE)
categorical_data$STAGE[is.na(categorical_data$STAGE)] <- "Stage I"
categorical_data$STAGE <- as.factor(categorical_data$STAGE)

lista<-levels(categorical_data$STAGE)
lista_early<-lista[-c(6:7)]
lista_late<-lista[-c(1:5)]

levels(categorical_data$STAGE)[levels(categorical_data$STAGE) %in% lista_early]<-"Stage_I_II"
levels(categorical_data$STAGE)[levels(categorical_data$STAGE) %in% lista_late]<-"Stage_III_IV"
categorical_data <- categorical_data[,-which(colnames(categorical_data) %in% c("AGE","Diagnosis.Age","Mutation.Count","Fraction.Genome.Altered"))]

#dummy <- dummyVars(" ~ .", data=categorical_data)
#newdata <- data.frame(predict(dummy, newdata = categorical_data))
data <- merge(x = categorical_data, y = datos_exp[,c("X","PDLIM5", "OSBPL10", "EGFR", "FRMD4B", "POFUT2", "HFE", "KCNMA1", "MOG", "ZNF44", "DOCK10", "VDR", "CA12", "RUNX1", "TNFRSF25", "SKAP2","OS_STATUS")], by = "X", all.x = TRUE)
categorical_data$STAGE <- as.factor(categorical_data$STAGE)

data$X <- NULL
data[is.na(data)] <- -1
data$ETHNICITY <- as.character(data$ETHNICITY)
data$NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT <- as.character(data$NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT)
data[is.na(data)] <- ''

train_index <- createDataPartition(data$OS_STATUS, p = .6, 
                                   list = FALSE, 
                                   times = 1)
train <- data [train_index,]
test  <- data [-train_index,]

trControl <- trainControl(method = "repeatedcv",
                          number = 3,
                          repeats = 0)
#repeats = 3)
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
# Utilizando las variables enteras, sin dummys 
rf


pred <-predict(rf, test, na.action = na.omit)
#pred <-predict(rf, test, na.action=na.omit)
# Extraigo la matriz de confusión

confMat<-confusionMatrix(pred, na.omit(test)$OS_STATUS)
#confMat3<-confusionMatrix(pred3, test3$Claudicacion_umbral3,positive="1")
confMat
draw_confusion_matrix(confMat)
# Importancia de las variables
importance <- varImp(rf,scale=TRUE)
print(importance)
plot(importance)


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
