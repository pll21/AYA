library(readr)

parentdir <- "Data"

mut_df <- read_delim("Data/data_mutations_extended_3.0.0.txt", 
                     "\t", escape_double = FALSE, comment = "#", 
                     trim_ws = TRUE)
clin_df <- read_delim("Data/data_clinical_sample_3.0.0.txt", 
                      "\t", escape_double = FALSE, comment = "#", 
                      trim_ws = TRUE)

head(mut_df$Tumor_Sample_Barcode)
head(clin_df$SAMPLE_ID)

levels(factor(clin_df$CANCER_TYPE))
glioma_clin <- clin_df[clin_df$CANCER_TYPE== "Glioma",]
glioma_clin <- glioma_clin[glioma_clin$SAMPLE_TYPE== "Primary",]

mut_partdf <- mut_df[,1:16]
mut_limited <- mut_partdf[,c("Hugo_Symbol", "Tumor_Sample_Barcode")]
mut_list <- levels(factor(mut_partdf$Hugo_Symbol))

##
head(glioma_clin)

glioma_dup1 <- glioma_clin[duplicated(glioma_clin$PATIENT_ID),]
glioma_dup2 <- glioma_clin[duplicated(glioma_clin$PATIENT_ID, fromLast = T),]
glioma_dup <- rbind(glioma_dup1, glioma_dup2)
glioma_dup <- glioma_dup[order(glioma_dup$PATIENT_ID, glioma_dup$AGE_AT_SEQ_REPORT),]

glioma_nodup <- glioma_clin[!(duplicated(glioma_clin$PATIENT_ID) | duplicated(glioma_clin$PATIENT_ID, fromLast = TRUE)), ]
glioma_nodup <- glioma_clin[!duplicated(glioma_clin$PATIENT_ID),]
glioma_nodup <- as.data.frame(glioma_nodup)
glioma_nodup$AGE_AT_SEQ_REPORT[glioma_nodup$AGE_AT_SEQ_REPORT == "<18"] <- "17"
glioma_nodup$AGE_AT_SEQ_REPORT[glioma_nodup$AGE_AT_SEQ_REPORT == ">89"] <- "90"

glioma_nodup$AGE_AT_SEQ_REPORT <- as.numeric(glioma_nodup$AGE_AT_SEQ_REPORT)
glioma_nodup$AYA <- "Adult"
glioma_nodup$AYA[glioma_nodup$AGE_AT_SEQ_REPORT <40] <- "AYA"

get_muts <- function(clinicalsubgroup, mut_limited){
  ## Filter mutation table by clinical group
  mut_filtered <- merge(clinicalsubgroup, mut_limited, by.x = "SAMPLE_ID", by.y = "Tumor_Sample_Barcode", all.y = F)
  
  ## Get mutation counts for clinical group
  
  mut_table <- table(mut_filtered$PATIENT_ID, mut_filtered$Hugo_Symbol)
  
  mut_counts <- as.data.frame(margin.table(mut_table, 2))
  return(mut_counts)
}

## Generate PCT
gliomamuts <- get_muts(clinicalsubgroup = glioma_nodup, mut_limited = mut_limited)
gliomamuts <- levels(factor(gliomamuts$Var1))

mut_filtered <- merge(glioma_nodup, mut_limited, by.x = "SAMPLE_ID", by.y = "Tumor_Sample_Barcode", all.y = F)

#
pct <- glioma_nodup
pct[gliomamuts] <- NA

names(pct)
for (i in 10:length(pct)){
  gene <- names(pct)[i]
  for(j in 1:length(pct$SAMPLE_ID)){
    genemuts <- mut_limited[mut_limited$Hugo_Symbol == gene,]
    
    if(pct$SAMPLE_ID[j] %in% genemuts$Tumor_Sample_Barcode){
      pct[j,i] <- 1
    } else {
      pct[j,i] <- 0
    }
  }
}

# write.csv(pct, file = file.path("Data", "PCT_mutations_GENIE.csv"), quote = T, row.names = F)


##Analyses
# pct <- as.data.frame(read_csv(file.path("Data", "PCT_mutations_GENIE.csv")))


names(pct)[1:10]
levels(factor(pct$CANCER_TYPE_DETAILED))
gliomas <- c(1:4, 6,7,11,12,14:17)
astro <- c(1,6,7,15)
pct_gliomas <- pct[pct$CANCER_TYPE_DETAILED %in%  levels(factor(pct$CANCER_TYPE_DETAILED))[gliomas],]
pct_gliomas$AYA[pct_gliomas$AGE_AT_SEQ_REPORT <18] <- "Child"
pct_qual <- pct_gliomas[,9:length(pct)]

pct_IDHastro <- pct_qual[(pct_gliomas$IDH1 ==1 | pct_gliomas$IDH2 ==1) & pct_gliomas$CANCER_TYPE_DETAILED %in% levels(factor(pct$CANCER_TYPE_DETAILED))[astro], ]
pct_ODG <- pct_qual[pct_gliomas$CANCER_TYPE_DETAILED == "Oligodendroglioma" | pct_gliomas$CANCER_TYPE_DETAILED == "Anaplastic Oligodendroglioma", ]
pct_GBM <- pct_qual[!(pct_gliomas$IDH1 ==1 | pct_gliomas$IDH2 ==1) & pct_gliomas$CANCER_TYPE_DETAILED %in%  levels(factor(pct$CANCER_TYPE_DETAILED))[gliomas], ]


table(pct_IDHastro$AYA)
table(pct_ODG$AYA)
table(pct_GBM$AYA)

# df <- pct_ODG

get_qualtables <- function (df, parentdir, outputfile){
  tablelist <- list()
  tables_withNA <- list()
  # runlabel <- format(Sys.time(), "%F %H-%M")
  # ifelse(!dir.exists(file.path(parentdir, runlabel)), dir.create(file.path(parentdir, runlabel)), F)
  
  for(i in 1:(ncol(df))){
    tablelist[[i]] <- table(as.logical(df[,i]), df$AYA )
    tables_withNA[[i]] <- table(as.logical(df[,i]), df$AYA 
                                , useNA = "always"
    )
    names(tablelist)[i] <- colnames(df)[i]
  }
  
  
  chisq_results <- list()
  hits <- vector(length = length(tablelist))
  
  # for(i in 1:length(tablelist)){
  #   glm_results[[i]] <- chisq.test(tablelist[[i]])
  #   if (chisq_results[[i]]$p.value < 0.1){
  #     hits[i] <- TRUE
  #   }
  # }
  
  for(i in 1:length(tablelist)){
      if(nrow(tablelist[[i]]) ==2 ){
        chisq_results[[i]] <- fisher.test(tablelist[[i]])
        if (chisq_results[[i]]$p.value < 0.1){
          hits[i] <- TRUE
        }
        names(chisq_results)[i] <- colnames(df)[i]
      }
  }

  
  currfile <- file.path(parentdir,outputfile)
  
  for(i in 1: length(tablelist)){
    if (hits[i]){
      write.table(paste("Gene: ", names(tablelist)[i]), 
                  currfile, append = T, quote = F, sep = ",", col.names = F, row.names = F )
      # write.table(names(tablelist)[i], 
      #             currfile, append = T, quote = F, sep = ",", col.names = F, row.names = F )
      write.table(chisq_results[[i]]$method, 
                  currfile, append = T, quote = F, sep = ",", col.names = F, row.names = F)
      write.table(paste("p-val = ", formatC(chisq_results[[i]]$p.value, format ="e", digits=2), sep = ""), 
                  currfile, append = T, quote = F, sep = ",", col.names = F, row.names = F)
      write.table(tables_withNA[[i]], currfile,append = T, quote = F, sep = ",", col.names = NA)
      write.table("", currfile, append = T, quote = F, sep = ",", col.names = F, row.names = F )
    }
  }
  
  # sink(file = outputfile)
  #   for (i in 1: length(tablelist)){
  #     cat("Variable: ", names(tablelist)[i], sep = "")
  #     cat("\n", chisq_results[[i]]$method, "\n", "p = ", chisq_results[[i]]$p.value, sep = "")
  #     print(tablelist[[i]])
  #     cat("\n\n", sep = "")
  # 
  #   }
  # sink()
  # return(tablelist)
}

get_qualtables(df= pct_ODG, parentdir= "Data", outputfile = "ODG_sig_results.csv")

get_qualtables(df= pct_IDHastro, parentdir= "Data", outputfile = "IDH_sig_results.csv")

pct_IDHastro_nochild <- pct_IDHastro[pct_IDHastro$AYA != "Child",]
get_qualtables(df= pct_IDHastro_nochild, parentdir= "Data", outputfile = "IDH_sig_results_noChild.csv")

get_qualtables(df= pct_GBM, parentdir= "Data", outputfile = "GBM_sig_results.csv")
get_qualtables(df= pct_GBM_childAYA, parentdir= "Data", outputfile = "GBMChildAYA_sig_results.csv")

