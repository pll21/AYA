SuppTablePanGlioma <- read.csv("Data/SuppTablePanGlioma.csv",
                               skip = 1)

panMerge <- as.data.frame(SuppTablePanGlioma)
colnames(panMerge)
panMerge$Type <- "GBM" 
panMerge$Type[panMerge$Study=="Brain Lower Grade Glioma"] <- "LGG"

summary(panMerge)
panMerge$AYA <- "Adult"
panMerge$AYA[panMerge$Age..years.at.diagnosis. < 40] <- "AYA"
panMerge$AYA[is.na(panMerge$Age..years.at.diagnosis.)] <- NA

panMerge_select <- panMerge[,colnames(panMerge) %in% c("Type", "Histology", "Grade",
                                                       "Age..years.at.diagnosis.",
                                                       "Gender", "Survival..months.",
                                                       "Karnofsky.Performance.Score",
                                                       "IDH.codel.subtype",
                                                       "MGMT.promoter.status",
                                                       "TERT.expression.status",
                                                       "NEJM",
                                                       "Supervised.DNA.Methylation.Cluster",
                                                       "YA"
)]

table(panMerge_select$IDH.codel.subtype, panMerge_select$Histology)
