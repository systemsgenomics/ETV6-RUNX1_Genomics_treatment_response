# processing and filtering the IgCaller results

# sample IDs
samples <- c("GE6345", "GE7324", "GE6347", "GE5340", "GE7318", "GE5338", "GE7304", "GE7305", "GE9342", "GE7301", "GE9321", "GE6339", "GE5312", "GET19344", "GE9336","GE2328", "GE3309", "GET20301", "GE8321", "GE8329", "GE9347", "GE0335", "GE8322", "GE0310", "GE3317", "GEK20315", "GE0319", "GE1329", "GE9348", "GE1312", "GE3320", "GE8315", "GE9344")

# generating tables into which the results from all samples can be collected
collective_csr <- matrix(ncol=5, nrow=1)
collective_transl <- matrix(ncol=5, nrow=1)
collective <- matrix(ncol=12, nrow=1)

# conducting the following steps for each sample
for (s in samples) {
  
  # reading in the filtered data
  tab <- read.csv(paste0(s, ".tumor_output_filtered.tsv"), sep="\t")
  
  # separating the CSR and translocations
  csr <- subset(tab, tab$Analysis=="CSR")
  transl <- subset(tab, tab$Analysis=="Oncogenic IG rearrangement")
  tab <- subset(tab, tab$Analysis!="Oncogenic IG rearrangement" & tab$Analysis!="CSR")
  tab <- subset(tab, !is.na(tab$Score))
  
  # reading in the IGH, IGK, and IGL tables, and keeping variants with at least 5 split and/or insert size reads supporting them
  igh <- read.csv(paste0(s, ".tumor_output_IGH.tsv"), sep="\t")
  igh_pass <- subset(igh, igh[,3]>=5 | igh[,4]>=5)
  igk <- read.csv(paste0(s, ".tumor_output_IGK.tsv"), sep="\t")
  igk_pass <- subset(igk, igk[,3]>=5 | igk[,4]>=5)
  igl <- read.csv(paste0(s, ".tumor_output_IGL.tsv"), sep="\t")
  igl_pass <- subset(igl, igl[,3]>=5 | igl[,4]>=5)
  
  # keeping rearrangements that pass both the tool filters and the read support filter
  rows_to_keep <- vector()
  for (i in 1:nrow(tab)){
    for (j in 1:nrow(igh_pass)){if (tab[i,2]==igh_pass[j,1] & as.numeric(tab[i,4])==as.numeric(igh_pass[j,20])){rows_to_keep <- append(rows_to_keep, i)}} 
    for (j in 1:nrow(igk_pass)){if (tab[i,2]==igk_pass[j,1] & as.numeric(tab[i,4])==as.numeric(igk_pass[j,20])){rows_to_keep <- append(rows_to_keep, i)}} 
    for (j in 1:nrow(igl_pass)){if (tab[i,2]==igl_pass[j,1] & as.numeric(tab[i,4])==as.numeric(igl_pass[j,20])){rows_to_keep <- append(rows_to_keep, i)}}}
  tab2 <- tab[rows_to_keep,]
  
  # writing out the tables
  write.table(tab2, paste0(s, ".IGH_IGK_IGL_alterations.tab"), quote=F, sep="\t", col.names=T, row.names=F)
  write.table(csr, paste0(s, ".CSR_alterations.tab"), quote=F, sep="\t", col.names=T, row.names=F)
  write.table(transl, paste0(s, ".oncogenic_IG_rearrengements.tab"), quote=F, sep="\t", col.names=T, row.names=F)
  
  # adding the patient ID into the first column of each table, and combining them to the collective tables
  patientID <- matrix(ncol=1, nrow=nrow(tab2))
  patientID[,1] <- s
  tab2 <- cbind(patientID, tab2)
  colnames(collective) <- colnames(tab2)
  collective <- rbind(collective, tab2)
  patientID <- matrix(ncol=1, nrow=nrow(csr))
  patientID[,1] <- s
  csr <- cbind(patientID, csr[,1:4])
  colnames(collective_csr) <- colnames(csr)
  collective_csr <- rbind(collective_csr, csr)
  patientID <- matrix(ncol=1, nrow=nrow(transl))
  patientID[,1] <- s
  transl <- cbind(patientID, transl[,1:4])
  colnames(collective_transl) <- colnames(transl)
  collective_transl <- rbind(collective_transl, transl)}

# writing the collective tables out
write.table(collective, "collective.IGH_IGK_IGL_alterations.tab", quote=F, sep="\t", col.names=T, row.names=F)
write.table(collective_csr, "collective.CSR_alterations.tab", quote=F, sep="\t", col.names=T, row.names=F)
write.table(collective_transl, "collective.oncogenic_IG_rearrengements.tab", quote=F, sep="\t", col.names=T, row.names=F)
