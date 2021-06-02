library(openxlsx)
library(utils)
library(stringr)
library(writexl)

HOME <- "../data/clinical_outcome/paad_tcga_pan_can_atlas_2018/"
RNA <- "data_RNA_Seq_v2_mRNA_median_all_sample_Zscores.txt"
SAMPLES <- "data_clinical_sample.txt"
PATIENTS <- "data_clinical_patient.txt"
FIRMA <- "../gene signature/firma_over_1.xlsx"
CLINICAL <- "../data/clinical_outcome/paad_tcga_pan_can_atlas_2018/paad_tcga_pan_can_atlas_2018_clinical_data.tsv"

datos_RNA<-read.delim(paste(HOME,RNA,sep = ""), sep="\t", dec=".",na.strings=c("","NA"))

datos_samples<-read.delim(paste(HOME,SAMPLES, sep = ""), sep="\t",na.strings=c("","NA"))

datos_patients<-read.delim(paste(HOME,PATIENTS, sep = ""), sep="\t",na.strings=c("","NA"))

datos_geneSign<-read.xlsx(FIRMA, sheet=1)

clinical_data <- read.table(file = CLINICAL, sep = '\t', header = TRUE)


datos_geneSign$Gene<-trimws(datos_geneSign$Gene, which = c("both"))
datos_geneSign$Gene[which(datos_geneSign$Gene=='MIEN1')]<-"C17orf37"

RNA_geneSign<-merge(x=datos_geneSign, y=datos_RNA, by.x="Gene", by.y="Hugo_Symbol", all.x=TRUE)

t1<-data.frame(t(RNA_geneSign))

rownames(t1)<-colnames(RNA_geneSign)
colnames(t1)<-RNA_geneSign$Gene

t1$Samples<- colnames(RNA_geneSign)
datos_exp<-t1[-c(1,2),]
rownames(datos_exp)<-NULL
datos_exp$Samples <- str_replace_all(datos_exp$Samples,'\\.', '-')

datos_exp<-merge(x=datos_exp, y=datos_samples[,c("SAMPLE_ID","PATIENT_ID")], by.x="Samples", by.y="SAMPLE_ID", all.x=TRUE)

datos_exp<-merge(x=datos_exp, y=datos_patients[,c("PATIENT_ID","AGE","SEX","RACE","ETHNICITY","NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT","AJCC_PATHOLOGIC_TUMOR_STAGE","WEIGHT", "OS_STATUS","OS_MONTHS","DSS_STATUS","DSS_MONTHS","DFS_STATUS","DFS_MONTHS","PFS_STATUS","PFS_MONTHS")], by.x="PATIENT_ID", by.y="PATIENT_ID", all.x=TRUE)

datos_exp<-merge(x = datos_exp, y = clinical_data[,c("Patient.ID","Diagnosis.Age","Mutation.Count","Prior.Diagnosis","Fraction.Genome.Altered")],by.x = "PATIENT_ID", by.y = "Patient.ID", all.x = TRUE)

write.csv(x = datos_exp, file = "ClinicalOutcomesDS/CO_atlas2018_over1.csv")
write_xlsx(datos_exp, "ClinicalOutcomesDS/CO_atlas2018_over1.xlsx")
