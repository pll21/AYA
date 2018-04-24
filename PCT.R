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

write.csv(pct, file = file.path("Data", "PCT_mutations_GENIE.csv"), quote = T, row.names = F)

table(pct$AYA, pct$IDH1)
