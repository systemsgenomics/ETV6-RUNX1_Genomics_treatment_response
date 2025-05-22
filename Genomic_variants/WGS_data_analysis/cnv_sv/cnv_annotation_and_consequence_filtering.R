# annotate CNV results and separate coding and non-coding events 

# reading in the table with all the SVannot results
tab <- read.table("collective_CNV_annotated.tab", header=F)

# looking into 3 cases which have whole genome duplications, and therefore a different normal copy number
# reading in the ASCAT bed files of these cases
ge5340 <- read.table("GE5340.ascatCNV.filtered.with_dups.bed")
ge5340 <- subset(ge5340, ge5340[,6]==4)
ge9348 <- read.table("GE9348.A.ascatCNV.filtered.with_dups.bed")
ge9348 <- subset(ge9348, ge9348[,6]==4)
ge9336 <- read.table("GE9336.ascatCNV.filtered.with_dups.bed")
ge9336 <- subset(ge9336, ge9336[,6]==3)
# and excluding the baseline copy number events, i.e., cn 4 in GE9348.A and GE5340, and 3 in GE9336
rows_to_remove <- vector()
for (i in 1:nrow(tab)) {
  for (j in 1:nrow(ge5340)) {
    if (isTRUE(tab[i,1]=="GE5340") & isTRUE(tab[i,2]==ge5340[j,1]) & isTRUE(tab[i,3]==ge5340[j,2]) & isTRUE(tab[i,4]==ge5340[j,3])) {rows_to_remove <- append(rows_to_remove, i)}}
  for (j in 1:nrow(ge9348)) {
    if (isTRUE(tab[i,1]=="GE9348") & isTRUE(tab[i,2]==ge9348[j,1]) & isTRUE(tab[i,3]==ge9348[j,2]) & isTRUE(tab[i,4]==ge9348[j,3])) {rows_to_remove <- append(rows_to_remove, i)}}
  for (j in 1:nrow(ge9336)) {
    if (isTRUE(tab[i,1]=="GE9336") & isTRUE(tab[i,2]==ge9336[j,1]) & isTRUE(tab[i,3]==ge9336[j,2]) & isTRUE(tab[i,4]==ge9336[j,3])) {rows_to_remove <- append(rows_to_remove, i)}}}
tab <- tab[-rows_to_remove,]

# checking which events have the copy number 2:0, i.e., which have the possibility to be LOH events instead of deletions
# reading in the LOH events
loh <- read.table("collective_loh_events.tab")
# separating the end coordinate from the format column
for (i in 1:nrow(loh)) {
  vector <- unlist(strsplit(loh[i,3], "END="))
  loh[i,3] <- vector[2]}

# excluding the LOH events which were decided to be kept as deletions since they are in whole genome duplicated cases
# and events that are annotated to be duplications since they affect sex chromosomes
lohrows <- vector()
for (i in 1:nrow(tab)) {
  loh2 <- subset(loh, loh[,1]==tab[i,2] & loh[,2]==tab[i,3] & loh[,3]==tab[i,4] & loh[,6]==tab[i,1])
  if (nrow(loh2)>0) {lohrows <- append(lohrows, i)}}
lohcnv <- tab[lohrows,]
cnv <- tab[-lohrows,]
rows_to_keep <- vector()
for (i in 1:nrow(lohcnv)) {
  if (lohcnv[i,5]=="DUP") {rows_to_keep <- append(rows_to_keep, i)}
  if (lohcnv[i,1]=="GE9348" & lohcnv[i,6]!="GE9348A_ascat") {rows_to_keep <- append(rows_to_keep, i)}
  if (lohcnv[i,1]=="GE9336") {rows_to_keep <- append(rows_to_keep, i)}
  if (lohcnv[i,1]=="GE5340") {rows_to_keep <- append(rows_to_keep, i)}}
rows_to_keep <- unique(rows_to_keep)
notloh <- lohcnv[rows_to_keep,]
loh <- lohcnv[-rows_to_keep,]
# changing the deletion status to LOH for LOH events
loh[,5] <- c("LOH")
tab <- rbind(cnv, notloh, loh)

# separating the rows with full annotations and writing the table out
full <- subset(tab, tab[,7]=="full")
write.table(full, "collective_CNV_annotated_only_full.tab", col.names=F, row.names=F, quote=F, sep="\t")

# separating the rows, the variants of which do not directly affect any genes
intergenic <- subset(tab, tab[,7]=="full" & tab[,8]==".")

# removing the rows which have have multiple genes on the same row, i.e., full annotations
# since all the genes are handled separately in the split rows
tab <- subset(tab, tab[,7]!="full")
# adding a consequence column to the table
cons <- matrix(ncol=1, nrow=nrow(tab))
for (i in 1:nrow(tab)){
  # adding a loss consequence if a deletion overlaps 100% with the coding region
  if (isTRUE(tab[i,5]=="DEL") & isTRUE(as.numeric(tab[i,12])>=100)){
    cons[i,1] <- "loss" }
  # adding a coding sequence consequence if a deletion overlaps 1-99% with the coding region
  if (isTRUE(tab[i,5]=="DEL") & isTRUE(as.numeric(tab[i,12])<100) & isTRUE(as.numeric(tab[i,12])>0)){
    cons[i,1] <- "coding_sequence_variant"}
  # adding a loh consequence if the variant is a LOH event
  if (isTRUE(tab[i,5]=="LOH") & isTRUE(as.numeric(tab[i,12])>0)) {
    cons[i,1] <- "loh"}
  # adding a gain consequence if a duplication overlaps 100% with the coding region
  if (isTRUE(tab[i,5]=="DUP") & isTRUE(as.numeric(tab[i,12])>=100)){
    cons[i,1] <- "gain"}
  # adding a coding sequence consequence if a deletion overlaps 1-99% with the coding region
  if (isTRUE(tab[i,5]=="DUP") & isTRUE(as.numeric(tab[i,12])<100) & isTRUE(as.numeric(tab[i,12])>0)){
    cons[i,1] <- "coding_sequence_variant"}}
tab <- cbind(tab, cons)

# separating the coding and non-coding variants
coding <- subset(tab, is.na(tab[,14])==FALSE)
noncoding <- subset(tab, is.na(tab[,14])==TRUE)
# including the relevant columns
coding <- coding[,c(1,2,3,4,5,6,8,14)] 
noncoding <- noncoding[,c(1,2,3,4,5,6,8,10,11,12,13)]
# and adding column names
colnames(coding) <- c("patientID", "chrom", "start", "end", "CNV_type", "callers", "gene", "consequence")
colnames(noncoding) <- c("patientID", "chrom", "start", "end", "CNV_type", "callers", "gene", "txStart", "txEnd", "CDS_overlap", "location")
# writing the tables out
write.table(coding, "collective_coding_CNVs.tab", quote=F, sep="\t", col.names=T, row.names=F)
write.table(noncoding, "collective_noncoding_CNVs.tab", quote=F, sep="\t", col.names=T, row.names=F)

# generating a bed file of the intergenic variants
# adding a chr in front of the chromosome number
for (i in 1:nrow(intergenic)){
  intergenic[i,2] <- paste0("chr", intergenic[i,2]) }
intergenic <- intergenic[,c(2,3,4,1,5)]  
write.table(intergenic, "collective_intergenic_CNVs.bed", quote=F, sep="\t", col.names=F, row.names=F)
