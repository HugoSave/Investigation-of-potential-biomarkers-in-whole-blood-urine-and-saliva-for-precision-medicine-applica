library(VennDiagram)
library(PerformanceAnalytics)
library(ggplot2)
library(limma)
library(reshape2)
library(readxl)
library(tidyr)
library(raster)
library(openxlsx)
library(zipR)
#chart.correlation

setwd("/")

protein_path = "proteinGroups_saliva-urine-blood.txt"
FDA_approved_biomarkers_path = "manualy_formated_FDA_approved_biomarkers.xls"
non_FDA_approved_biomarkers_path = "manualy_formated_non_FDA_approved_biomarkers.xls"

Read.MaxQuant <- function(file){
  
  return(read.table(file, header = T,sep="\t", fill=T))
}


filter_CON_REV <- function(prot_data){
  
  filtered_data <- prot_data[-which(substr(prot_data$Majority.protein.IDs, 0, 3)== "CON"|substr(prot_data$Majority.protein.IDs, 0, 3)== "REV"),]
  
  return (filtered_data)
}

amount_protein_samples <- function(prot_data){
  
  
  blood_DF <- prot_data[, grepl("Intensity.b", names(prot_data))]
  saliva_DF <- prot_data[, grepl("Intensity.s", names(prot_data))]
  urine_DF <- prot_data[, grepl("Intensity.u", names(prot_data))]
  
  all_samples_intensity <- prot_data[, grepl("Intensity.", names(prot_data))]

  
}

binary_func <- function(protdataList, firstSampleColNr, onlyFDA){
  #View(protdataList[[1]])
  
  if (onlyFDA == TRUE){
    
    protdataList[[1]] <- protdataList[[1]][protdataList[[1]]$approved_by_FDA == TRUE | protdataList[[1]]$not_approved_by_FDA == TRUE,]
    protdataList[[2]] <- protdataList[[2]][protdataList[[2]]$approved_by_FDA == TRUE | protdataList[[2]]$not_approved_by_FDA == TRUE,]
    protdataList[[3]] <- protdataList[[3]][protdataList[[3]]$approved_by_FDA == TRUE | protdataList[[3]]$not_approved_by_FDA == TRUE,]
    
     }
  blood_sum <- protdataList[[1]][,c("intensity", "gene_name", "uniprot_name")]
  blood_sum$intensity[blood_sum$intensity != 0] <- 1
  blood_sum$intensity[is.na(blood_sum$intensity)] <- 0
  
  saliva_sum <- protdataList[[2]][,c("intensity", "gene_name", "uniprot_name")]
  saliva_sum$intensity[saliva_sum$intensity != 0] <- 1
  saliva_sum$intensity[is.na(saliva_sum$intensity)] <- 0

  urine_sum <- protdataList[[3]][,c("intensity", "gene_name", "uniprot_name")]
  urine_sum$intensity[urine_sum$intensity != 0] <- 1
  urine_sum$intensity[is.na(urine_sum$intensity)] <- 0
  
  #View(blood_sum)
  
  
  #blood_sum <- apply(protdata[firstSampleColNr:firstSampleColNr+8],1,sum, na.rm=TRUE)
  #saliva_sum <- apply(protdata[salive_start: salive_start +8], 1, sum,na.rm=TRUE)
  #urine_sum <- apply(protdata[urine_start: urine_start +8], 1, sum,na.rm=TRUE)
  return (list(blood_sum,saliva_sum,urine_sum))
  B_S_U_combinded <- cbind(blood_sum[, "intensity"],saliva_sum[, "intensity"], urine_sum[, "intensity"])
  #View(B_S_U_combinded)
  B_S_U_TrueFalse <- B_S_U_combinded != 0
  
  #Ombandlar TRUE till 1 och FALSE till 0
  #B_S_U_0101 <- B_S_U_0101*1
  
  return (B_S_U_TrueFalse)
  }
#ean_of_patients(intensity,1,c("Blood101","Blood102","Blood103"),9)

#send in all samples with their intensity
function_on_patients_bad_code <- function(protein_data_local, Func){
  
  protein_data_local[protein_data_local == 0] <- NA
  blood_column_101 <- apply(protein_data_local[1:3], 1, Func, na.rm =TRUE)
  blood_column_102 = apply(protein_data_local[4:6], 1, Func, na.rm =TRUE)
  blood_column_103 = apply(protein_data_local[7:9], 1, Func, na.rm =TRUE) 
  saliva_column_101 = apply(protein_data_local[10:12], 1, Func, na.rm =TRUE)
  saliva_column_102 = apply(protein_data_local[13:15], 1, Func, na.rm =TRUE)
  saliva_column_103 = apply(protein_data_local[16:18], 1, Func, na.rm =TRUE)
  urine_column_101 = apply(protein_data_local[19:21], 1, Func, na.rm =TRUE)
  urine_column_102 = apply(protein_data_local[22:24], 1, Func, na.rm =TRUE)
  urine_column_103 = apply(protein_data_local[25:27], 1, Func, na.rm =TRUE)
  
  combinded_tests= cbind(blood_column_101,blood_column_102,blood_column_103,
                         saliva_column_101,saliva_column_102,saliva_column_103,
                         urine_column_101,urine_column_102,urine_column_103
                         )
  combinded_tests[is.na(combinded_tests)] <- NA
  return(as.data.frame(combinded_tests))
}

#send in the patient catagoriezed intensities
patien_to_BSU_categories_median <- function(data, Func){
  
  data <- as.data.frame(data)
  
  

  #print(assign(paste("blood", funcname,sep=""), apply(data[1:3], 1, Func, na.rm = TRUE)))
  blood_median <- apply(data[1:3], 1, Func, na.rm = TRUE)
  saliva_median <- apply(data[4:6], 1, Func, na.rm = TRUE)
  urine_median <- apply(data[7:9], 1, Func, na.rm = TRUE)
  
  rest_columns <- data[10:length(colnames(data))]
  
  
  combined_BSU_test <- cbind(blood_median,saliva_median,urine_median, rest_columns)
  
  #capture.output(Func) == capture.output(cv)
  if (identical(Func,cv)){
    colnames(combined_BSU_test) <- c("blood_cv", "saliva_cv", "urine_cv", colnames(rest_columns))
  }
  else if(identical(Func,sd)){
    colnames(combined_BSU_test) <- c("blood_sd", "saliva_sd", "urine_sd", colnames(rest_columns))
  }
  else if(identical(Func,mean)){
    colnames(combined_BSU_test) <- c("blood_mean", "saliva_mean", "urine_mean", colnames(rest_columns))
  }
  else if (identical(Func,median)){
    colnames(combined_BSU_test) <- c("blood_median", "saliva_median", "urine_median", colnames(rest_columns))
  }
  else {
    colnames(combined_BSU_test) <- c("blood_unknown_func", "saliva_unknown_func", "urine_unknown_func", colnames(rest_columns))
  }
  
  combined_BSU_test[is.na(combined_BSU_test)] <- NA
  return (combined_BSU_test)
}


mean_of_patients_func <- function(protdata,firstSampleColNr, column_names, amountOfSamples){
  
  samplesPerPatient = amountOfSamples/length(column_names)
  iterations = amountOfSamples/samplesPerPatient
  colNr = firstSampleColNr
  
  Samples_mean <- data.frame()
  #Samples_mean[,ncol(Samples_mean)+1] <- NA
  #Samples_mean[, 4] <- c("s","d","a")
  print(dim(Samples_mean))
  
  sample_list = list()
  
  print(length(Samples_mean))
  for (i in range(iterations)){
    sample_mean <- apply(protdata[colNr:colNr+samplesPerPatient], 1, mean)
    colNr <- colNr + samplesPerPatient
    df <- data.frame(sample_mean)
    Samples_mean <- cbind(Samples_mean, df)
  }
  
  return (Samples_mean)
  
  
}


plot_venndiagram <- function(protdataList, firstColNummer, onlyFDA){
  
  binary_data <- binary_func(protdataList, firstColNummer, onlyFDA)
  
  
  blood_sum <- binary_data[[1]]$intensity
  saliva_sum <- binary_data[[2]]$intensity
  urine_sum <- binary_data[[3]]$intensity
  B_S_U_combinded <- cbind(blood_sum,saliva_sum, urine_sum)
  View(B_S_U_combinded)
  B_S_U_TrueFalse <- B_S_U_combinded != 0
  
  print(vennCounts(B_S_U_TrueFalse))
  vennDiagram(vennCounts(B_S_U_TrueFalse))
}


#cv_from_patient <- function_on_patients_bad_code(intensity, cv)


create_melted_B_S_U_list_with_sd <- function(median_BSU, sd_BSU){
  BSU_melted <- melt(median_BSU, id = c("gene_name", "uniprot_name"), value.name = "intensity")
  BSU_sd_melted <- melt(sd_BSU, id=c("gene_name", "uniprot_name"), value.name = "sd" )
  
  B_melted <- BSU_melted[grep("blood", BSU_melted$variable),]
  S_melted <- BSU_melted[grep("saliva", BSU_melted$variable),]
  U_melted <- BSU_melted[grep("urine", BSU_melted$variable),]
  #View(B_melted)
  B_sd <- BSU_sd_melted[grep("blood", BSU_sd_melted$variable),]
  S_sd <- BSU_sd_melted[grep("saliva", BSU_sd_melted$variable),]
  U_sd <- BSU_sd_melted[grep("urine", BSU_sd_melted$variable),]
  
  B_melted_sd <- cbind(B_melted, B_sd[4])
  S_melted_sd <- cbind(S_melted, S_sd[4])
  U_melted_sd <- cbind(U_melted, U_sd[4])
  
  #names(B_melted_sd) <- c("gene_name", "sample_type", "intensity", "sd")
  #names(S_melted_sd) <- c("gene_name", "sample_type", "intensity", "sd")
  #names(U_melted_sd) <- c("gene_name", "sample_type", "intensity", "sd")
  
  
  return(list(B_melted_sd, S_melted_sd, U_melted_sd))
  
}

plot_intensity_to_gene_name <- function(sample_with_sd){
  
  
  sample_with_sd <- sample_with_sd[complete.cases(sample_with_sd$intensity),]
  
  sample_with_sd$color_group <- "not known biomarker"
  
  sample_with_sd[ sample_with_sd$not_approved_by_FDA ,  "color_group"]  <- "potential biomarker"
  
  sample_with_sd[ sample_with_sd$approved_by_FDA ,  "color_group"]  <- "approved biomarker"
  
  sample_with_sd$color_group <- factor(sample_with_sd$color_group, levels = c("approved biomarker", "potential biomarker", "not known biomarker"))
  #View(sample_with_sd)
  #sort uniqe????+
  
  View(sample_with_sd)
  
  approved_genes <- sample_with_sd$gene_name[sample_with_sd$approved_by_FDA == TRUE]
  
  approved_genes <- approved_genes[sample(NROW(approved_genes),NROW(approved_genes)*(0.3))]
  
  View(approved_genes)
  theme_1 <- theme(axis.title.y = element_text(size = 8, colour = "black"),
                   axis.title.x = element_blank(), 
                   axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0, face = "bold", colour = "black"),
                   axis.text.y = element_text(size = 14, colour = "black"),
                   plot.title = element_text(size=12, face="bold")) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(plot.margin = unit(c(1,0.5,1,0.5), "cm"))
  
  ggplot(sample_with_sd, aes(x= reorder(gene_name, -intensity), y= log10(intensity))) + 
    geom_jitter(aes(color=color_group))  +
    geom_errorbar(aes(ymin=log10(intensity-sd), ymax=log10(intensity+sd), color=color_group)) + 
    scale_color_manual(values = c( "green","blue", "grey")) +
    scale_x_discrete(breaks =approved_genes) +
    theme_1 +
    labs(x = "Protein gene names(only randomly chosen \"approved biomarkers\" are shown)", y = "log10(Intensity)", color = "FDA classification")
  
  #levels = gene_name[sort(unique(-intensity, na.last = NA))]
  #ggplot(sample_with_sd, aes(x= gene_name, y= log10(intensity))) + geom_point() + geom_errorbar(aes(ymin=log10(intensity-sd), ymax=log10(intensity+sd)))
  
   }

separate_uniprot_names <- function(data) {
  #replace or and x with ,
  data$uniprot_name <- gsub("or",",",gsub("[+]", ",", data$uniprot_name))
  #remove whitespaces
  data$uniprot_name <- gsub( " ", "", data$uniprot_name)
  
  #delete rows with ? and ~
  data <- data[data$uniprot_name != "?",]
  data <- data[data$uniprot_name != "~",]
  
  data <- separate_rows(data, uniprot_name)
  
  return (data)
  
}

human_readable_table <- function(each_sample_intensity_w_gene, approved_biomarkers, not_approved_biomarkers){
  
  #Vilka gener vill vi ha? J?mf?r med biomarkers
  
  combinded_biomarkers_names <- rbind(approved_biomarkers[,1:2], not_approved_biomarkers[,1:2])
  data <- each_sample_intensity_w_gene[each_sample_intensity_w_gene$uniprot_name %in% combinded_biomarkers_names$uniprot_name,]
  #data <- combinded_biomarkers_names[each_sample_intensity_w_gene$uniprot_name %in% combinded_biomarkers_names$uniprot_name,"uniprot_name"]
  
  
  data_true_false <- data[, grep("Intensity.", colnames(data))] > 0
  data[, grep("Intensity.", colnames(data))] <- data_true_false[]
  data[data[] == TRUE] <- 1
  
  #View(data)
  
  number_in_blood_samples <- apply(data[, grep("Intensity.b", colnames(data))] , 1, sum, na.rm = TRUE)
  
  number_in_saliva_samples <- apply(data[, grep("Intensity.s", colnames(data))] , 1, sum, na.rm = TRUE)
  
  number_in_urine_samples <- apply(data[, grep("Intensity.u", colnames(data))] , 1, sum, na.rm = TRUE)
  
  blood_protein_in_patient_101 <- apply(data[,grep("Intensity.b_101", colnames(data))], 1 ,sum, na.rm=TRUE)
  blood_protein_in_patient_101[blood_protein_in_patient_101 > 0] <- 1
  blood_protein_in_patient_102 <- apply(data[,grep("Intensity.b_102", colnames(data))], 1 ,sum, na.rm=TRUE)
  blood_protein_in_patient_102[blood_protein_in_patient_102 > 0] <- 1
  blood_protein_in_patient_103 <- apply(data[,grep("Intensity.b_103", colnames(data))], 1 ,sum, na.rm=TRUE)
  blood_protein_in_patient_103[blood_protein_in_patient_103 > 0] <- 1
  
  saliva_protein_in_patient_101 <- apply(data[,grep("Intensity.s_101", colnames(data))], 1 ,sum, na.rm=TRUE)
  saliva_protein_in_patient_101[saliva_protein_in_patient_101 > 0] <- 1
  saliva_protein_in_patient_102 <- apply(data[,grep("Intensity.s_102", colnames(data))], 1 ,sum, na.rm=TRUE)
  saliva_protein_in_patient_102[saliva_protein_in_patient_102 > 0] <- 1
  saliva_protein_in_patient_103 <- apply(data[,grep("Intensity.s_103", colnames(data))], 1 ,sum, na.rm=TRUE)
  saliva_protein_in_patient_103[saliva_protein_in_patient_103 > 0] <- 1
  
  urine_protein_in_patient_101 <- apply(data[,grep("Intensity.u_101", colnames(data))], 1 ,sum, na.rm=TRUE)
  urine_protein_in_patient_101[urine_protein_in_patient_101 > 0] <- 1
  urine_protein_in_patient_102 <- apply(data[,grep("Intensity.u_102", colnames(data))], 1 ,sum, na.rm=TRUE)
  urine_protein_in_patient_102[urine_protein_in_patient_102 > 0] <- 1
  urine_protein_in_patient_103 <- apply(data[,grep("Intensity.u_103", colnames(data))], 1 ,sum, na.rm=TRUE)
  urine_protein_in_patient_103[urine_protein_in_patient_103 > 0] <- 1
  
  #blood_protein_in_nr_patients = cbind(blood_protein_in_patient_101,blood_protein_in_patient_102,blood_protein_in_patient_103)
  blood_protein_in_nr_patients = apply(cbind(blood_protein_in_patient_101,blood_protein_in_patient_102,blood_protein_in_patient_103), 1, sum)
  
  #blood_protein_in_nr_patients = cbind(blood_protein_in_patient_101,blood_protein_in_patient_102,blood_protein_in_patient_103)
  saliva_protein_in_nr_patients = apply(cbind(saliva_protein_in_patient_101,saliva_protein_in_patient_102,saliva_protein_in_patient_103), 1, sum)
  #saliva_protein_in_nr_patients = apply(blood_protein_in_nr_patients, 1, sum)
  
  urine_protein_in_nr_patients = apply(cbind(urine_protein_in_patient_101,urine_protein_in_patient_102,urine_protein_in_patient_103), 1, sum)
  
  combined_DF <- as.data.frame(cbind(blood_protein_in_nr_patients, number_in_blood_samples, saliva_protein_in_nr_patients, number_in_saliva_samples, urine_protein_in_nr_patients, number_in_urine_samples))
  print(class(combined_DF))
  print(class(blood_protein_in_nr_patients))
  combined_DF <- combined_DF %>% unite(col = "Blood, protein accurance in different patients / in number of samples", c(blood_protein_in_nr_patients,number_in_blood_samples), sep = "/", remove = TRUE) 
  combined_DF <- combined_DF %>% unite(col = "Saliva, protein accurance in different patients / in number of samples", c(saliva_protein_in_nr_patients,number_in_saliva_samples), sep = "/" ,remove = TRUE)
  combined_DF <- combined_DF %>% unite(col = "Urine, protein accurance in different patients / in number of samples", c(urine_protein_in_nr_patients,number_in_urine_samples), sep = "/", remove = TRUE)
  
  combinded_biomarkers_names <- combinded_biomarkers_names[!duplicated(combinded_biomarkers_names$uniprot_name), ]
  
  #View(combinded_biomarkers_names)
  #View(combinded_biomarkers_names)
  #combined_DF$protein_name_FDA <- combinded_biomarkers_names$Protein_name_FDA
  combined_DF$uniprot_name <- data$uniprot_name
  temp <- combined_DF
  #View(temp)
  combined_DF$Protein_name_FDA <- NA
  
  #Denna rad sorterar inte
  common_names <- combinded_biomarkers_names[combinded_biomarkers_names$uniprot_name %in% combined_DF$uniprot_name ,c("Protein_name_FDA", "uniprot_name")]
  #Sort them so that they can be combined
  common_names <- common_names[order(common_names$uniprot_name),]
  combined_DF <- combined_DF[order(combined_DF$uniprot_name),]
 # View(common_names)
  combined_DF$Protein_name_FDA <- common_names$Protein_name_FDA
  #combined_DF$Protein_name_FDA <- combinded_biomarkers_names[combined_DF$uniprot_name %in% combinded_biomarkers_names$uniprot_name,"Protein_name_FDA"]
  #View(combined_DF)
  #View(combined_DF)
  #apply(combinded_biomarkers_names, MARGIN = 2, if(combinded_biomarkers_names$uniprot_name %in% combined_DF$uniprot_name) combinded_biomarkers_names$protein_name_fda)
    
 #   combined_DF$uniprot_name %in% combinded_biomarkers_names$Protein_name_FDA
  return(combined_DF)
  #combined_DF <- cbind(combined_DF %>% unite(col = "Saliva, protein accurance in different patients / in number of samples", c(saliva_protein_in_nr_patients,number_in_saliva_samples), sep = "_" ))
 # combined_DF <- cbind(combined_DF %>% unite(col = "Urine, protein accurance in different patients / in number of samples"), c(urine_protein_in_nr_patients,number_in_urine_samples), sep = "_")
  
  
  #View(blood_1010)
  
}


cv_table_under_percent <- function(data, limit, colNR1, colNr2){
  
  data[,colNR1:colNr2][data[,colNR1:colNr2] > limit] <- NA
  #View(data)
  Na_rows <- apply(a[,colNR1:colNr2], 1, function(x) all(is.na(x)))
  data <- data[!Na_rows,]
  
  return(data)

}


correlatiation_graph <- function(data){
  
  chart.Correlation(log10(data), histogram = TRUE)
  
}

trim_and_add_FDA_protein_names <- function(data, FDA_names){
  
  FDA_names <- FDA_names[!duplicated(FDA_names$Protein_name_FDA),]
  data <- data[data$uniprot_name %in% FDA_names$uniprot_name,]
  common_names <- FDA_names[FDA_names$uniprot_name %in% data$uniprot_name,]
  common_names <- common_names[order(common_names$uniprot_name),]
  View(common_names)
  data <- data[order(data$uniprot_name),]
  data$protein_name_fda <- NA
  data$protein_name_fda <- common_names$Protein_name_FDA
  
  return(data)
  
  
}




FDA_approved_data <- read_excel(FDA_approved_biomarkers_path)
View(FDA_approved_data)
Non_FDA_approved_data <- read_excel(non_FDA_approved_biomarkers_path)

extended_FDA_approved_data <- separate_uniprot_names(FDA_approved_data)
extended_non_FDA_approved_data <- separate_uniprot_names(Non_FDA_approved_data)
combinded_biomarker_names <- rbind(extended_FDA_approved_data[,1:2], extended_non_FDA_approved_data[,1:2])
combinded_biomarker_names <- combinded_biomarker_names[!duplicated(combinded_biomarker_names$uniprot_name), ]


protein_data <- Read.MaxQuant(protein_path)
filtered_data <- filter_CON_REV(protein_data)
filtered_data$gene <- gsub("PE=.*", "", gsub(".*GN=", "", filtered_data$Majority.protein.IDs))
#Get gene name between | and | 
filtered_data$uniprot_name <- substr(filtered_data$Majority.protein.IDs, 4, 9)

intensity <- filtered_data[128:154]
#Remove two outlierproteins who have very low intensity < 1000. All other intensities are around 10^6 or over
intensity$Intensity.u_103_3[intensity$Intensity.u_103_3 < 1000] <- 0

intensity$gene_name <- filtered_data$gene
intensity$uniprot_name <- filtered_data$uniprot_name

#intensity <- filtered_data$gene




#sample_groups_binary <- binary_func(intensity, 1)

mean_from_patients <- function_on_patients_bad_code(intensity, mean)
mean_from_patients$gene_name <- filtered_data$gene
mean_from_patients$uniprot_name <-filtered_data$uniprot_name

sd_from_patients <- function_on_patients_bad_code(intensity, sd)
sd_from_patients$gene_name <- filtered_data$gene
sd_from_patients$uniprot_name <- filtered_data$uniprot_name

cv_from_patients <- function_on_patients_bad_code(intensity, cv)
cv_from_patients$gene_name <- filtered_data$gene
cv_from_patients$uniprot_name <- filtered_data$uniprot_name
cv_patients_biomarkers <-trim_and_add_FDA_protein_names(cv_from_patients,combinded_biomarker_names)


BSU_median <- patien_to_BSU_categories_median(mean_from_patients, median)
BSU_sd <- patien_to_BSU_categories_median(mean_from_patients, sd)
BSU_cv <- patien_to_BSU_categories_median(mean_from_patients, cv)




B_S_U_list <- create_melted_B_S_U_list_with_sd(BSU_median,BSU_sd)

B_S_U_list[[1]]$approved_by_FDA <- FALSE
B_S_U_list[[2]]$approved_by_FDA <- FALSE
B_S_U_list[[3]]$approved_by_FDA <- FALSE

B_S_U_list[[1]]$not_approved_by_FDA <- FALSE
B_S_U_list[[2]]$not_approved_by_FDA <- FALSE
B_S_U_list[[3]]$not_approved_by_FDA <- FALSE



#Replaces False with True on rows where uniprot names match from FDA list
B_S_U_list[[1]][B_S_U_list[[1]]$uniprot_name %in% extended_FDA_approved_data$uniprot_name, "approved_by_FDA"] <- TRUE
B_S_U_list[[2]][B_S_U_list[[2]]$uniprot_name %in% extended_FDA_approved_data$uniprot_name, "approved_by_FDA"] <- TRUE
B_S_U_list[[3]][B_S_U_list[[3]]$uniprot_name %in% extended_FDA_approved_data$uniprot_name, "approved_by_FDA"] <- TRUE

B_S_U_list[[1]][B_S_U_list[[1]]$uniprot_name %in% extended_non_FDA_approved_data$uniprot_name, "not_approved_by_FDA"] <- TRUE
B_S_U_list[[2]][B_S_U_list[[2]]$uniprot_name %in% extended_non_FDA_approved_data$uniprot_name, "not_approved_by_FDA"] <- TRUE
B_S_U_list[[3]][B_S_U_list[[3]]$uniprot_name %in% extended_non_FDA_approved_data$uniprot_name, "not_approved_by_FDA"] <- TRUE


#N?gon skilnad p? B_S_U_list och B_melted_with_sd. hitta den

plot_venndiagram(B_S_U_list,1, onlyFDA = FALSE)
plot_venndiagram(B_S_U_list,1, onlyFDA = TRUE)
plot_intensity_to_gene_name(B_S_U_list[[2]])

table_summary_BSU<- human_readable_table(intensity,extended_FDA_approved_data,extended_non_FDA_approved_data)

#combine summare table with cv
add_to_summary <- BSU_cv[BSU_cv$uniprot_name %in% table_summary_BSU$uniprot_name,]
add_to_summary <- add_to_summary[order(add_to_summary$uniprot_name),]
table_summary_BSU <- table_summary_BSU[order(table_summary_BSU$uniprot_name),]
table_summary_BSU <- cbind(table_summary_BSU, add_to_summary[1:3])
rm(add_to_summary)

#cv_table(intensity)

#only_fda_intensity <- trim_and_add_FDA_protein_names(intensity, combinded_biomarker_names)

cv_under_percent <- cv_table_under_percent(cv_patients_biomarkers,25, 1, 9)

#table common saliva and blood
#BSU_1010_list <- binary_func(B_S_U_list, 1, onlyFDA = TRUE)
#BSU_trimmed <- trim_and_add_FDA_protein_names(BSU_median,combinded_biomarker_names )
BS_table <- table_summary_BSU[table_summary_BSU$`Saliva, protein accurance in different patients / in number of samples` != "0/0" & table_summary_BSU$`Blood, protein accurance in different patients / in number of samples`!= "0/0", ]
BS_table <- subset(BS_table, select= - c(`Urine, protein accurance in different patients / in number of samples`, `urine_cv`))
#table common urine and blood
BU_table <- table_summary_BSU[table_summary_BSU$`Urine, protein accurance in different patients / in number of samples` != "0/0" & table_summary_BSU$`Blood, protein accurance in different patients / in number of samples`!= "0/0", ]
BU_table <- subset(BU_table, select= - c(`Saliva, protein accurance in different patients / in number of samples`,`saliva_cv` ))

S_not_in_blood <- table_summary_BSU[table_summary_BSU$`Blood, protein accurance in different patients / in number of samples` == "0/0" &table_summary_BSU$`Saliva, protein accurance in different patients / in number of samples` != "0/0",]
S_not_in_blood <- subset(S_not_in_blood, select= -c(`Urine, protein accurance in different patients / in number of samples`, `urine_cv`, `blood_cv`))

U_not_in_blood <- table_summary_BSU[table_summary_BSU$`Blood, protein accurance in different patients / in number of samples` == "0/0" &table_summary_BSU$`Urine, protein accurance in different patients / in number of samples` != "0/0",]
U_not_in_blood <- subset(U_not_in_blood, select= -c(`Saliva, protein accurance in different patients / in number of samples`, `saliva_cv`, `blood_cv`))

SU_not_in_blood <- table_summary_BSU[table_summary_BSU$`Blood, protein accurance in different patients / in number of samples`== "0/0",]
SU_not_in_blood <- subset(SU_not_in_blood, select= -c(`blood_cv`, `blood_cv`))


#new stuff
intensity_FDA <-trim_and_add_FDA_protein_names(intensity, combinded_biomarker_names) 
View(intensity_FDA)
data_needed_for_graph <- intensity_FDA[,c(grep("Intensity.b", intensity_FDA), "gene_name", "uniprot_name", "protein_name_fda")]




