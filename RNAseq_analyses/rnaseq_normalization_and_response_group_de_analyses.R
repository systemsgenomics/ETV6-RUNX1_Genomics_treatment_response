# normalize the combined RNA-seq data
# conduct response group DE analyses

# required packages
library(limma)
library(Glimma)
library(edgeR)
library(Homo.sapiens)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(gridExtra)

# first, reading in the table with the response group info for both cohorts
resp_info <- read.table("combined_response_groups.txt", header=T)

# reading in the read counts for first sample cohort
counts_matrix <- read.table("ALL_GEX_n=315_read_counts.txt", header=T)
rownames(counts_matrix) <- counts_matrix[,1]
# and the response group information
response <- read.table("uu_rnaseq_mrd.txt", header=T)
samples <- as.vector(response[,1])
# subsetting the matrix to include only the ER cases with MRD data
counts_matrix <- counts_matrix[,samples]

# reading in the read counts for second sample cohort
counts_matrix2 <- read.table("full_salmon_gene_counts.tsv", header=T)
rownames(counts_matrix2) <- counts_matrix2[,1]
# reading in the sample info
sample_info <- read.table("rna_master-2.txt", header=T, sep="\t")
# choosing only ER samples from the table
er_samples <- subset(sample_info, sample_info$Subtype=="ER")
# saving the EGA IDs of ER cases into a vector
ega <- as.vector(er_samples[,1])
# subsetting the count matrix to only include ER cases
samples <- colnames(counts_matrix2)
for (s in 1:length(samples)) {
  if (samples[s] %in% ega) {
    er_cols <- append(er_cols, s)}}
counts_matrix2 <- counts_matrix2[,er_cols]
# then adding the gene IDs as the row names
rownames(counts_matrix2) <- counts_matrix2[,1]
counts_matrix2 <- counts_matrix2[,-1]

# making sure that the genes and row order are identical between the
identical(rownames(counts_matrix), rownames(counts_matrix2))
# which is TRUE

# saving the sample IDs as vectors
uusamples <- colnames(counts_matrix)
gesamples <- colnames(counts_matrix2)

# combining the two tables
counts_matrix <- cbind(counts_matrix, counts_matrix2)
# generating the DGE object from the combined count matrix
data <- DGEList(counts_matrix)

# organizing sample information, adding info on the sequencing batch, sex, and response to therapy 
samples <- rownames(data$samples)
batch <- vector()
uuinfo <- read.csv("ALL_GEX_n=315.csv")
for (s in samples) {
  if (s %in% uusamples) {
    uuinfo2 <- subset(uuinfo, uuinfo[,1]==s)
    batch <- append(batch, uuinfo2[1,22])}
  if (s %in% gesamples) {
    sample <- subset(sample_info, sample_info$ALT_ID_StudySpecific==s)
    batch <- append(batch, sample[1,6])}}
for (b in 1:length(batch)) {
  batch[b] <- gsub("-", "", batch[b])}
data$samples$batch <- batch

# reading in the table describing the sexes of the gepard and uu patients, and adding that info to the DGE object
phenotype <- read.csv("ALL_GEX_n=315.csv", header=T)
sexes <- read.table("sexes_gepard_samples.txt", sep="\t")
sex <- vector()
for (s in samples) {
  if (s %in% uusamples) {
    sample <- subset(phenotype, phenotype[,1]==s)
    sex <- append(sex, sample[1,24])}
  if (s %in% gesamples) {
    sample <- subset(sexes, sexes[,1]==s)
    if (sample[1,2]==1) {sex <- append(sex, "Male")}
    if (sample[1,2]==2) {sex <- append(sex, "Female")}}}
data$samples$sex <- sex

# adding d15, d29 and d79 response info
d15_response <- vector()
d29_response <- vector()
d79_response <- vector()
for (s in samples) {
  resp2 <- subset(resp_info, resp_info[,1]==s)
  if (is.na(resp2[1,2]) & s!="GE0325") {d15_response <- append(d15_response, "unknown")}
  if (is.na(resp2[1,3]) & s!="GE0325") {d29_response <- append(d29_response, "unknown")}
  if (is.na(resp2[1,4]) & s!="GE0325") {d79_response <- append(d79_response, "unknown")}
  if (!is.na(resp2[1,2])) {d15_response <- append(d15_response, resp2[1,2])}
  if (!is.na(resp2[1,3])) {d29_response <- append(d29_response, resp2[1,3])}
  if (!is.na(resp2[1,4])) {d79_response <- append(d79_response, resp2[1,4])}
  if (s=="GE0325") {
    d15_response <- append(d15_response, "slow")
    d29_response <- append(d29_response, "slow")
    d79_response <- append(d79_response, "unknown")}}
data$samples$d15_response <- d15_response
data$samples$d29_response <- d29_response
data$samples$d79_response <- d79_response

# also adding two additional columns for d29 comparisons: mrd pos vs neg, and fast+intermediate vs slow
d29_response2 <- vector()
d29_response3 <- vector()
for (s in samples) {
  resp2 <- subset(resp_info, resp_info[,1]==s)
  if (is.na(resp2[1,3]) & s!="GE0325") {d29_response2 <- append(d29_response2, "unknown")}
  if (is.na(resp2[1,3]) & s!="GE0325") {d29_response3 <- append(d29_response3, "unknown")}
  if (isTRUE(resp2[1,3]=="fast")) {d29_response2 <- append(d29_response2, "mrd_neg")}
  if (isTRUE(resp2[1,3]=="slow" | resp2[1,3]=="intermediate")) {d29_response2 <- append(d29_response2, "mrd_pos")}
  if (isTRUE(resp2[1,3]=="fast" | resp2[1,3]=="intermediate")) {d29_response3 <- append(d29_response3, "fastmed")}
  if (isTRUE(resp2[1,3]=="slow")) {d29_response3 <- append(d29_response3, "slow")}
  if (s=="GE0325") {
    d29_response2 <- append(d29_response2, "mrd_pos")
    d29_response3 <- append(d29_response3, "slow")}}
data$samples$d29_response2 <- d29_response2
data$samples$d29_response3 <- d29_response3

# adding symbols and geneIDs for the ENSEMBL genes in the counts matrix
geneid <- rownames(data)
kts <- keytypes(Homo.sapiens)
genes <- AnnotationDbi::select(Homo.sapiens, keys=geneid, columns=c("SYMBOL", "GENEID"), keytype="ENSEMBL")
# adding the genes into the DGE list
genes <- genes[!duplicated(genes$ENSEMBL),]
data$genes <- genes

# counts are transformed from the raw scale, using the counts per million function to normalize the read counts
cpm <- cpm(data)
# also computing the log-CPM values
lcpm <- cpm(data, log=T)
# checking the mean and median of library sizes
L <- mean(data$samples$lib.size) * 1e-6
M <- median(data$samples$lib.size) * 1e-6
c(L,M)
# and the lcpm summary
summary(lcpm)

# removing genes lowly expressed across the response group comparisons
table(rowSums(data$counts==0)==51)
keep.exprs1 <- filterByExpr(data, group=d15_response)
keep.exprs2 <- filterByExpr(data, group=d29_response)
keep.exprs3 <- filterByExpr(data, group=d79_response)
for (i in 1:length(keep.exprs1)) {
  if (isFALSE(keep.exprs1[i]) & isTRUE(keep.exprs2[i]) | isFALSE(keep.exprs1[i]) & isTRUE(keep.exprs3[i])) {keep.exprs1[i] <- TRUE}}
data <- data[keep.exprs1,,keep.lib.sizes=FALSE]

# normalising gene expression distributions, and visualizing the effects of normalization
pdf("TMM_normalization_effects_combo.pdf")
par(mfrow=c(2,1))
boxplot(lcpm, las=2, cex.axis=0.7)
title(main="A. Unnormalized data", ylab="log-cpm")
data <- calcNormFactors(data, method="TMM")
cpm <- cpm(data)
lcpm <- cpm(data, log=T)
boxplot(lcpm, las=2, cex.axis=0.7)
title(main="B. Normalized data", ylab="log-cpm")
dev.off()

# looking into unsupervised clustering of samples
pdf("MDS_response_groups_and_sex_and_batch_combo.pdf")
par(mfrow=c(3,2))
col.group1 <- as.factor(d15_response)
levels(col.group1) <- c("#a1d76a", "#b2182b", "grey68")
col.group2 <- as.factor(d29_response)
levels(col.group2) <- c("#a1d76a", "#f1a340", "#b2182b")
col.group3 <- as.factor(d79_response)
levels(col.group3) <- c("#a1d76a", "#b2182b", "grey68")
col.group1 <- as.character(col.group1)
col.group2 <- as.character(col.group2)
col.group3 <- as.character(col.group3)
col.sex <- as.factor(sex)
levels(col.sex) <- c("#01665e", "#984ea3")
col.sex <- as.character(col.sex)
col.batch <- as.factor(batch)
levels(col.batch) <- c("darkolivegreen3", "indianred1", "#5ab4ac", "powderblue", "chocolate", "goldenrod")
col.batch <- as.character(col.batch)
plotMDS(lcpm, labels=d15_response, col=col.group1, cex=0.5, cex.axis=0.7, cex.lab=0.7, cex.main=0.8, main="A. Day 15 response")
plotMDS(lcpm, labels=d29_response, col=col.group2, cex=0.5, cex.axis=0.7, cex.lab=0.7, cex.main=0.8, main="B. Day 29 response")
plotMDS(lcpm, labels=d79_response, col=col.group3, cex=0.5, cex.axis=0.7, cex.lab=0.7, cex.main=0.8, main="C. Day 79 response")
plotMDS(lcpm, labels=sex, col=col.sex, cex=0.5, cex.axis=0.7, cex.lab=0.7, cex.main=0.8, main="D. Sex")
plotMDS(lcpm, labels=batch, col=col.batch, cex=0.5, cex.axis=0.7, cex.lab=0.7, cex.main=0.8, main="E. Batch")
dev.off()

# writing out the lcpm and cpm matrices
write.table(lcpm, "lcpm_matrix_combo.tab", sep="\t", quote=F, col.names=T, row.names=T)
write.table(cpm, "cpm_matrix_combo.tab", sep="\t", quote=F, col.names=T, row.names=T)

# differential expression analyses
# adjusting for batch in each

# d29 response group comparison
# first creating a design matrix 
design <- model.matrix(~0+d29_response+batch)
# then generating a contrast matrix
contr.matrix <- makeContrasts(SlowVsFast=d29_responseslow-d29_responsefast, SlowVsIntermediate=d29_responseslow-d29_responseintermediate, IntermediateVsFast=d29_responseintermediate-d29_responsefast, levels=colnames(design))
# removing heteroscedascity from count data
pdf("mean-variance_trend_d29.pdf", width=10)
par(mfrow=c(1,2))
v <- voom(data, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit29 <- eBayes(vfit)
plotSA(efit29, main="Final model: Mean-variance trend")
dev.off()
# examining the top DE genes and the number of DE genes
topTable(efit29)
summary(decideTests(efit29))

# d15 response group comparison
# first creating a design matrix 
design <- model.matrix(~0+d15_response+batch)
# then generating a contrast matrix
contr.matrix <- makeContrasts(SlowVsFast=d15_responseslow-d15_responsefast, levels=colnames(design))
# removing heteroscedascity from count data
pdf("mean-variance_trend_d15.pdf", width=10)
par(mfrow=c(1,2))
v <- voom(data, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit15 <- eBayes(vfit)
plotSA(efit15, main="Final model: Mean-variance trend")
dev.off()
# examining the top DE genes and the number of DE genes
topTable(efit15)
summary(decideTests(efit15))

# d79 response group comparison
# first creating a design matrix 
design <- model.matrix(~0+d79_response+batch)
# then generating a contrast matrix
contr.matrix <- makeContrasts(PosVsNeg=d79_responsepositive-d79_responsenegative, levels=colnames(design))
# removing heteroscedascity from count data
pdf("mean-variance_trend_d79.pdf", width=10)
par(mfrow=c(1,2))
v <- voom(data, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit79 <- eBayes(vfit)
plotSA(efit79, main="Final model: Mean-variance trend")
dev.off()
# examining the top DE genes and the number of DE genes
topTable(efit79)
summary(decideTests(efit79))

# d29 MRD pos (slow and intermediate) vs neg (fast) response group comparison
# first creating a design matrix 
design <- model.matrix(~0+d29_response2+batch)
# then generating a contrast matrix
contr.matrix <- makeContrasts(PosVsNeg=d29_response2mrd_pos-d29_response2mrd_neg, levels=colnames(design))
# removing heteroscedascity from count data
pdf("mean-variance_trend_d29_posVSneg.pdf", width=10)
par(mfrow=c(1,2))
v <- voom(data, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit29_2 <- eBayes(vfit)
plotSA(efit29_2, main="Final model: Mean-variance trend")
dev.off()
# examining the top DE genes and the number of DE genes
topTable(efit29_2)
summary(decideTests(efit29_2))

# d29 slow vs fast + intermediate response group comparison
# first creating a design matrix 
design <- model.matrix(~0+d29_response3+batch)
# then generating a contrast matrix
contr.matrix <- makeContrasts(SlowVsFastMed=d29_response3slow-d29_response3fastmed, levels=colnames(design))
# removing heteroscedascity from count data
pdf("mean-variance_trend_d29_slowVSfast_med.pdf", width=10)
par(mfrow=c(1,2))
v <- voom(data, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit29_3 <- eBayes(vfit)
plotSA(efit29_3, main="Final model: Mean-variance trend")
dev.off()
# examining the top DE genes and the number of DE genes
topTable(efit29_3)
summary(decideTests(efit29_3))

# writing out the efit files
write.fit(efit29, "efit_d29.tab", results=NULL, sep="\t", adjust="BH")
write.fit(efit15, "efit_d15.tab", results=NULL, sep="\t", adjust="BH")
write.fit(efit79, "efit_d79.tab", results=NULL, sep="\t", adjust="BH")
write.fit(efit29_2, "efit_d29_posVSneg.tab", results=NULL, sep="\t", adjust="BH")
write.fit(efit29_3, "efit_d29_slowVSfast_med.tab", results=NULL, sep="\t", adjust="BH")
