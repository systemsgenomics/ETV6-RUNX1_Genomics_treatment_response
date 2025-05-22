# annotate CNV results and separate coding and non-coding events 

# collective SV results
tab <- read.table("collective_SV_annotated.tab", header=F)

# separating rows the variants of which do not affect any genes
intergenic <- subset(tab, tab[,8]=="full" & tab[,9]==".")
# removing the rows which have have multiple genes on the same row, i.e., full annotations
# since all the genes are handled separately in the split rows
tab <- subset(tab, tab[,8]!="full")

# adding a consequence column to the table
cons <- matrix(ncol=1, nrow=nrow(tab))
for (i in 1:nrow(tab)){
  # adding a loss consequence if a deletion overlaps 100% with the coding region
  if (tab[i,6]=="DEL" & as.numeric(tab[i,13])==100){
    cons[i,1] <- "loss"}
  # adding a coding sequence consequence if a deletion overlaps 1-99% with the coding region
  if (tab[i,6]=="DEL" & as.numeric(tab[i,13])<100 & as.numeric(tab[i,13])>0){
    cons[i,1] <- "coding_sequence_variant"}
  # adding a gain consequence if a duplication overlaps 100% with the coding region
  if (tab[i,6]=="DUP" & as.numeric(tab[i,13])==100){
    cons[i,1] <- "gain"}
  # adding a coding sequence consequence if a deletion overlaps 1-99% with the coding region
  if (tab[i,6]=="DUP" & as.numeric(tab[i,13])<100 & as.numeric(tab[i,13])>0){
    cons[i,1] <- "coding_sequence_variant"}
  # adding a coding sequence consequence if a breakend locates into the transcript region
  if (tab[i,6]=="BND" & as.numeric(tab[i,3])>as.numeric(tab[i,11]) & as.numeric(tab[i,3])<as.numeric(tab[i,12])){
    cons[i,1] <- "coding_sequence_variant"}
  # adding a coding sequence consequence if the whole gene affected does not locate into the inversion region
  # based on the txStart and txEnd coordinates
  if (tab[i,6]=="INV" & as.numeric(tab[i,11])<as.numeric(tab[i,3]) & as.numeric(tab[i,12])>as.numeric(tab[i,3])){
    cons[i,1] <- "coding_sequence_variant" }
  if (tab[i,6]=="INV" & as.numeric(tab[i,11])<as.numeric(tab[i,4]) & as.numeric(tab[i,12])>as.numeric(tab[i,4])){
    cons[i,1] <- "coding_sequence_variant" }} 
tab <- cbind(tab, cons)

# separating the coding and non-coding variants
coding <- subset(tab, is.na(tab[,15])==FALSE)
noncoding <- subset(tab, is.na(tab[,15])==TRUE)

# and including the relevant columns
coding <- coding[,c(1,2,3,4,6,7,9,14,15)] 
noncoding <- noncoding[,c(1,2,3,4,6,7,9,11,12,13,14)]
# determining column names
colnames(coding) <- c("patientID", "chrom", "start", "end", "SV_type", "callers", "gene", "location", "consequence")
colnames(noncoding) <- c("patientID", "chrom", "start", "end", "SV_type", "callers", "gene", "txStart", "txEnd", "CDS_overlap", "location")
# writing the tables out
write.table(coding, "collective_coding_SVs.tab", quote=F, sep="\t", col.names=T, row.names=F)
write.table(noncoding, "collective_noncoding_SVs.tab", quote=F, sep="\t", col.names=T, row.names=F)

# generating a separate file of the intergenic variants
for (i in 1:nrow(intergenic)){
  intergenic[i,2] <- paste0("chr", intergenic[i,2])} 
intergenic <- intergenic[,c(2,3,4,1,6)]  
write.table(intergenic, "collective_intergenic_SVs.bed", quote=F, sep="\t", col.names=F, row.names=F)
