# filter Vardict variants based on strand bias and VEP consequences
library(sjmisc)

# reading in the af 0.05 filtered variants
tab <- read.table("collective.SNV_indel.af0.05.vardict.filtered.tab")

# adding columns to the tables to indicate the strand ratio of the varint supporting reads
newcols <- matrix(ncol=2, nrow=nrow(tab))
for (i in 1:nrow(tab)) {
  strands <- unlist(strsplit(tab[i,16], ":"))
  newcols[i,1] <- (as.numeric(strands[1])/as.numeric(tab[i,12]))
  newcols[i,2] <- (as.numeric(strands[2])/as.numeric(tab[i,12]))}
tab <- cbind(tab, newcols)

# filtering out the variants with either strand ratio < 0.01
# or < 0.05 and less than 5 reads
remove <- vector()
for (i in 1:nrow(tab)) {
  strands <- unlist(strsplit(tab[i,16], ":"))
  if (tab[i,22]<0.01 | tab[i,23]<0.01 | tab[i,22]<0.05 & as.numeric(strands[1])<=5 |tab[i,23]<0.05 & as.numeric(strands[2])<=5) {remove <- append(remove, i)}}
tab <- tab[-remove,]

# rounding the strand ratio values
for (i in 1:nrow(tab)) {
  tab[i,22] <- round(as.numeric(tab[i,22]), digits=2)
  tab[i,23] <- round(as.numeric(tab[i,23]), digits=2)}

# subsetting to only SNVs and indels
tab <- subset(tab, tab[,6] != "DEL" & tab[,6] != "DUP" & tab[,6]!="INV")

# reading in the PON vardict results
pon <- read.table("collective_pon_SNV_indel.vardict.tab")

# determining the prevalence of variants in the tumor and PON samples
rec <- unique(tab[,1:4])
newcols <- matrix(ncol=2, nrow=nrow(rec))
for (i in 1:nrow(rec)) {
  tab2 <- subset(tab, tab[,1]==rec[i,1] & tab[,2]==rec[i,2] & tab[,3]==rec[i,3] & tab[,4]==rec[i,4])
  newcols[i,1] <- length(unique(tab2[,21]))
  pon2 <- subset(pon, pon[,1]==rec[i,1] & pon[,2]==rec[i,2] & pon[,3]==rec[i,3] & pon[,4]==rec[i,4])
  newcols[i,2] <- length(unique(pon2[,21]))}
rec <- cbind(rec, newcols)
newcols <- matrix(ncol=2, nrow=nrow(rec))
for (i in 1:nrow(rec)) {
  newcols[i,1] <- as.numeric(rec[i,5])/262
  newcols[i,2] <- as.numeric(rec[i,6])/38}
rec <- cbind(rec, newcols)
# determining the variants present in more than 1 case in normal panel
rec1 <- subset(rec, as.numeric(rec[,6]) > 1)
# removing those variants from the table
excl <- subset(tab, tab[,2] %in% rec1[,2])
excl2 <- unique(excl[,1:4])
keep <- vector()
for (i in 1:nrow(excl)) {
  test <- subset(rec1, rec1[,1]==excl[i,1] & rec1[,2]==excl[i,2] & rec1[,3]==excl[i,3] & rec1[,4]==excl[i,4])
  if (nrow(test)==0) {keep <- append(keep, i)}}
keep <- excl[keep,]
tab <- subset(tab, !tab[,2] %in% rec1[,2])
tab <- rbind(tab, keep)

# determining column names and writing out the table
colnames(tab) <- c("chr", "pos", "ref", "alt", "filter", "variant_type", "qual", "dp", "af", "hiaf", "adjaf", "vd", "hicnt", "pmean", "mq", "varbias", "sbf", "consequence", "impact", "gene", "caseid", "strand1_vd_ratio", "strand2_vd_ratio")
write.table(tab, "collective.SNV_indel.af0.05.vardict.annotated.tab", sep="\t", quote=F, col.names=T, row.names=F)

# determining the high/moderate impact variants
csq <- subset(tab, tab[,19]=="HIGH" | tab[,19]=="MODERATE")
write.table(csq, "collective.SNV_indel.af0.05.vardict.csq_filt.tab", sep="\t", quote=F, col.names=T, row.names=F)
