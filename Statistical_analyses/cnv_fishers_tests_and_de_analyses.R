# Fisher's exact tests on the combined conumee (methylation) CNV and WGS CNV data
library(sjmisc)

# reading in the response table
resp <- read.table("sample_responses.txt", header=T)

# reading in the genes recurrently affected by CNVs the WGS cohort
genes1 <- read.table("recurrently_affected_genes_>3_by_coding_CNV_dups.tab", header=T)
genes2 <- read.table("recurrently_affected_genes_>3_by_coding_CNV_dels.tab", header=T)
genes3 <- read.table("recurrently_affected_genes_>5_by_coding_CNVs.tab", header=T)
genetab <- rbind(genes1, genes2, genes3)
genetab <- unique(genetab)
# combining all genes into one vector
genes <- append(genes1[,7], genes2[,7])
genes <- append(genes, genes3[,7])
genes <- unique(genes)

# WGS CNV results
nc <- read.table("CNVs_affecting_nc_genes.tab")
cnv <- read.table("collective_coding_CNVs.tab", header=T)
# combining the tables
gene_type <- matrix(ncol=1, nrow=nrow(cnv))
gene_type[,1] <- c("protein_coding")
cnv <- cbind(cnv[,c(1,2,5,7)], gene_type, cnv[,8])
nc <- nc[,1:6]
colnames(cnv) <- colnames(nc)
getab <- rbind(cnv, nc)
getab <- getab[,c(1:4,6)]

# panelseq and methylation array results
tab <- read.table("collective_CNV_filt_annot_per_gene.tab", header=T)
# simplifying the table
tab <- tab[,c(15,1,5,7,6)]

# combining the tables
colnames(getab) <- colnames(tab)
tab <- rbind(tab, getab)
# removing LOH events
tab <- subset(tab, tab[,3]!="LOH")
# correcting DUP into AMP
tab[,3] <- gsub("DUP", "AMP", tab[,3])
# including only coding consequences
tab <- subset(tab, tab[,5]=="loss" | tab[,5]=="nc_gene_loss" | tab[,5]=="coding_sequence_variant" | tab[,5]=="gain" | tab[,5]=="nc_gene_gain" | tab[,5]=="nc_gene_sequence_variant")
# subsetting it to only include the genes of interest
tab <- subset(tab, tab[,4] %in% genes)
tab <- unique(tab)

# adding response columns to the table
newcols <- matrix(ncol=4, nrow=nrow(tab))
for (i in 1:nrow(tab)) {
  resp2 <- subset(resp, resp[,1]==tab[i,1])
  newcols[i,1] <- resp2[1,2]
  newcols[i,2] <- resp2[1,3]
  newcols[i,3] <- resp2[1,4]
  if (isTRUE(str_contains(tab[i,1], "ALL_"))) {
    newcols[i,4] <- "array"}
  if (isTRUE(str_contains(tab[i,1], "GE"))) {
    newcols[i,4] <- "WGS"}
  if (isFALSE(str_contains(tab[i,1], "GE")) & isFALSE(str_contains(tab[i,1], "ALL_"))) {
    newcols[i,4] <- "panel"}}
colnames(newcols) <- c("d15_response", "d29_response", "d79_response", "cohort")
tab <- cbind(tab, newcols)

# generating vectors of the samples in each responder group
fast15 <- subset(resp, resp$d15_response=="fast")
slow15 <- subset(resp, resp$d15_response=="slow")
fast29 <- subset(resp, resp$d29_response=="fast")
med29 <- subset(resp, resp$d29_response=="intermediate")
slow29 <- subset(resp, resp$d29_response=="slow")
medslow29 <- subset(resp, resp$d29_response=="slow" | resp$d29_response=="intermediate")
fastmed29 <- subset(resp, resp$d29_response=="fast" | resp$d29_response=="intermediate")
fast79 <- subset(resp, resp$d79_response=="negative")
slow79 <- subset(resp, resp$d79_response=="positive")

# generating a table of the genes to include in the fisher's tests
genes <- tab[,c(4,2)]
genes <- unique(genes)
# excluding X chromosome genes
genes <- subset(genes, genes[,2]!="X")

## conducting the fisher's tests for deletions first

# generating a table into which the results can be collected
newcols <- matrix(ncol=24, nrow=nrow(genes))
colnames(newcols) <- c("d15_p", "d29all_p", "d29fm_p", "d29ms_p", "d79_p", "d15_top_prevalence", "d29all_top_prevalence", "d29fm_top_prevalence", "d29ms_top_prevalence", "d79_top_prevalence", "d15_frac_fast", "d15_frac_slow", "d29all_frac_fast", "d29all_frac_med", "d29all_frac_slow", "d29fm_frac_fast", "d29fm_frac_slow", "d29ms_frac_fast", "d29ms_frac_slow", "d79_frac_fast", "d79_frac_slow", "n_WGS", "n_array", "n_panel")
genes <- cbind(genes, newcols)
rownames(genes) <- genes[,1]
for (i in 1:nrow(genes)) {
  tab2 <- subset(tab, tab[,3]=="DEL" & tab[,4]==genes[i,1])
  # d15
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,6]=="slow")
  fast <- subset(tab2, tab2[,6]=="fast")
  test[1,1] <- length(unique(slow[,1]))
  test[2,1] <- length(unique(fast[,1]))
  test[1,2] <- (nrow(slow15)-length(unique(slow[,1])))
  test[2,2] <- (nrow(fast15)-length(unique(fast[,1])))
  fisher <- fisher.test(test)
  genes$d15_p[i] <- round(fisher$p, digits=4)
  genes$d15_frac_slow[i] <- round(length(unique(slow[,1]))/nrow(slow15), digits=3)
  genes$d15_frac_fast[i] <- round(length(unique(fast[,1]))/nrow(fast15), digits=3)
  # d29 all
  test <- matrix(ncol=2, nrow=3)
  slow <- subset(tab2, tab2[,7]=="slow")
  med <- subset(tab2, tab2[,7]=="intermediate")
  fast <- subset(tab2, tab2[,7]=="fast")
  test[1,1] <- length(unique(slow[,1]))
  test[2,1] <- length(unique(med[,1]))
  test[3,1] <- length(unique(fast[,1]))
  test[1,2] <- (nrow(slow29)-length(unique(slow[,1])))
  test[2,2] <- (nrow(med29)-length(unique(med[,1])))
  test[3,2] <- (nrow(fast29)-length(unique(fast[,1])))
  fisher <- fisher.test(test)
  genes$d29all_p[i] <- round(fisher$p, digits=4)
  genes$d29all_frac_slow[i] <- round(length(unique(slow[,1]))/nrow(slow29), digits=3)
  genes$d29all_frac_med[i] <- round(length(unique(med[,1]))/nrow(med29), digits=3)
  genes$d29all_frac_fast[i] <- round(length(unique(fast[,1]))/nrow(fast29), digits=3)
  # d29 fast + intermediate vs slow
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,7]=="slow")
  fast <- subset(tab2, tab2[,7]=="fast" | tab2[,7]=="intermediate")
  test[1,1] <- length(unique(slow[,1]))
  test[2,1] <- length(unique(fast[,1]))
  test[1,2] <- (nrow(slow29)-length(unique(slow[,1])))
  test[2,2] <- (nrow(fastmed29)-length(unique(fast[,1])))
  fisher <- fisher.test(test)
  genes$d29fm_p[i] <- round(fisher$p, digits=4)
  genes$d29fm_frac_slow[i] <- round(length(unique(slow[,1]))/nrow(slow29), digits=3)
  genes$d29fm_frac_fast[i] <- round(length(unique(fast[,1]))/nrow(fastmed29), digits=3)
  # d29 intermediate + slow vs fast
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,7]=="slow" | tab2[,7]=="intermediate")
  fast <- subset(tab2, tab2[,7]=="fast")
  test[1,1] <- length(unique(slow[,1]))
  test[2,1] <- length(unique(fast[,1]))
  test[1,2] <- (nrow(medslow29)-length(unique(slow[,1])))
  test[2,2] <- (nrow(fast29)-length(unique(fast[,1])))
  fisher <- fisher.test(test)
  genes$d29ms_p[i] <- round(fisher$p, digits=4)
  genes$d29ms_frac_slow[i] <- round(length(unique(slow[,1]))/nrow(medslow29), digits=3)
  genes$d29ms_frac_fast[i] <- round(length(unique(fast[,1]))/nrow(fast29), digits=3)
  # d79
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,8]=="positive")
  fast <- subset(tab2, tab2[,8]=="negative")
  test[1,1] <- length(unique(slow[,1]))
  test[2,1] <- length(unique(fast[,1]))
  test[1,2] <- (nrow(slow79)-length(unique(slow[,1])))
  test[2,2] <- (nrow(fast79)-length(unique(fast[,1])))
  fisher <- fisher.test(test)
  genes$d79_p[i] <- round(fisher$p, digits=4)
  genes$d79_frac_slow[i] <- round(length(unique(slow[,1]))/nrow(slow79), digits=3)
  genes$d79_frac_fast[i] <- round(length(unique(fast[,1]))/nrow(fast79), digits=3)}

# adding info to the columns describing the group with highest prevalence of alterations
for (i in 1:nrow(genes)) {
  if (isTRUE(genes$d15_frac_fast[i] < genes$d15_frac_slow[i])) {genes$d15_top_prevalence[i] <- c("slow")}
  if (isTRUE(genes$d15_frac_fast[i] > genes$d15_frac_slow[i])) {genes$d15_top_prevalence[i] <- c("fast")}
  if (isTRUE(genes$d29all_frac_fast[i] > genes$d29all_frac_slow[i] & genes$d29all_frac_fast[i] > genes$d29all_frac_med[i])) {genes$d29all_top_prevalence[i] <- c("fast")}
  if (isTRUE(genes$d29all_frac_med[i] > genes$d29all_frac_slow[i] & genes$d29all_frac_fast[i] < genes$d29all_frac_med[i])) {genes$d29all_top_prevalence[i] <- c("intermediate")}
  if (isTRUE(genes$d29all_frac_fast[i] < genes$d29all_frac_slow[i] & genes$d29all_frac_slow[i] > genes$d29all_frac_med[i])) {genes$d29all_top_prevalence[i] <- c("slow")}
  if (isTRUE(genes$d29fm_frac_fast[i] < genes$d29fm_frac_slow[i])) {genes$d29fm_top_prevalence[i] <- c("slow")}
  if (isTRUE(genes$d29fm_frac_fast[i] > genes$d29fm_frac_slow[i])) {genes$d29fm_top_prevalence[i] <- c("fast")}
  if (isTRUE(genes$d29ms_frac_fast[i] < genes$d29ms_frac_slow[i])) {genes$d29ms_top_prevalence[i] <- c("slow")}
  if (isTRUE(genes$d29ms_frac_fast[i] > genes$d29ms_frac_slow[i])) {genes$d29ms_top_prevalence[i] <- c("fast")}
  if (isTRUE(genes$d79_frac_fast[i] < genes$d79_frac_slow[i])) {genes$d79_top_prevalence[i] <- c("slow")}
  if (isTRUE(genes$d79_frac_fast[i] > genes$d79_frac_slow[i])) {genes$d79_top_prevalence[i] <- c("fast")}}

# adding info about the numbers in each dataset
for (i in 1:nrow(genes)) {
  tab2 <- subset(tab, tab[,3]=="DEL" & tab[,4]==genes[i,1])
  wgs <- subset(tab2, tab2$cohort=="WGS")
  array <- subset(tab2, tab2$cohort=="array")
  panel <- subset(tab2, tab2$cohort=="panel")
  genes$n_WGS[i] <- length(unique(wgs[,1]))
  genes$n_array[i] <- length(unique(array[,1]))
  genes$n_panel[i] <- length(unique(panel[,1]))}

# removing genes with 0 hits in any cohort, or less than 5 hits altogether
remove <- vector()
for (i in 1:nrow(genes)) {
  if (genes$n_WGS[i]==0 | genes$n_array[i]==0 | genes$n_panel[i]==0 | sum(genes[i,24:26]) <= 5) {
    remove <- append(remove, i)}}
del <- genes[-remove,]

# writing the table out
write.table(del, "del_fishers_test_results.tab", sep="\t", quote=F, col.names=T, row.names=F)


## continuing with amplifications

# genes of interest
genes <- tab[,c(4,2)]
genes <- unique(genes)
# excluding X chromosome genes
genes <- subset(genes, genes[,2]!="X")

# generating a table into which the results can be collected
newcols <- matrix(ncol=24, nrow=nrow(genes))
colnames(newcols) <- c("d15_p", "d29all_p", "d29fm_p", "d29ms_p", "d79_p", "d15_top_prevalence", "d29all_top_prevalence", "d29fm_top_prevalence", "d29ms_top_prevalence", "d79_top_prevalence", "d15_frac_fast", "d15_frac_slow", "d29all_frac_fast", "d29all_frac_med", "d29all_frac_slow", "d29fm_frac_fast", "d29fm_frac_slow", "d29ms_frac_fast", "d29ms_frac_slow", "d79_frac_fast", "d79_frac_slow", "n_WGS", "n_array", "n_panel")
genes <- cbind(genes, newcols)
rownames(genes) <- genes[,1]
for (i in 1:nrow(genes)) {
  tab2 <- subset(tab, tab[,3]=="AMP" & tab[,4]==genes[i,1])
  # d15
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,6]=="slow")
  fast <- subset(tab2, tab2[,6]=="fast")
  test[1,1] <- length(unique(slow[,1]))
  test[2,1] <- length(unique(fast[,1]))
  test[1,2] <- (nrow(slow15)-length(unique(slow[,1])))
  test[2,2] <- (nrow(fast15)-length(unique(fast[,1])))
  fisher <- fisher.test(test)
  genes$d15_p[i] <- round(fisher$p, digits=4)
  genes$d15_frac_slow[i] <- round(length(unique(slow[,1]))/nrow(slow15), digits=3)
  genes$d15_frac_fast[i] <- round(length(unique(fast[,1]))/nrow(fast15), digits=3)
  # d29 all
  test <- matrix(ncol=2, nrow=3)
  slow <- subset(tab2, tab2[,7]=="slow")
  med <- subset(tab2, tab2[,7]=="intermediate")
  fast <- subset(tab2, tab2[,7]=="fast")
  test[1,1] <- length(unique(slow[,1]))
  test[2,1] <- length(unique(med[,1]))
  test[3,1] <- length(unique(fast[,1]))
  test[1,2] <- (nrow(slow29)-length(unique(slow[,1])))
  test[2,2] <- (nrow(med29)-length(unique(med[,1])))
  test[3,2] <- (nrow(fast29)-length(unique(fast[,1])))
  fisher <- fisher.test(test)
  genes$d29all_p[i] <- round(fisher$p, digits=4)
  genes$d29all_frac_slow[i] <- round(length(unique(slow[,1]))/nrow(slow29), digits=3)
  genes$d29all_frac_med[i] <- round(length(unique(med[,1]))/nrow(med29), digits=3)
  genes$d29all_frac_fast[i] <- round(length(unique(fast[,1]))/nrow(fast29), digits=3)
  # d29 fast + intermediate vs slow
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,7]=="slow")
  fast <- subset(tab2, tab2[,7]=="fast" | tab2[,7]=="intermediate")
  test[1,1] <- length(unique(slow[,1]))
  test[2,1] <- length(unique(fast[,1]))
  test[1,2] <- (nrow(slow29)-length(unique(slow[,1])))
  test[2,2] <- (nrow(fastmed29)-length(unique(fast[,1])))
  fisher <- fisher.test(test)
  genes$d29fm_p[i] <- round(fisher$p, digits=4)
  genes$d29fm_frac_slow[i] <- round(length(unique(slow[,1]))/nrow(slow29), digits=3)
  genes$d29fm_frac_fast[i] <- round(length(unique(fast[,1]))/nrow(fastmed29), digits=3)
  # d29 fast vs intermediate + slow
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,7]=="slow" | tab2[,7]=="intermediate")
  fast <- subset(tab2, tab2[,7]=="fast")
  test[1,1] <- length(unique(slow[,1]))
  test[2,1] <- length(unique(fast[,1]))
  test[1,2] <- (nrow(medslow29)-length(unique(slow[,1])))
  test[2,2] <- (nrow(fast29)-length(unique(fast[,1])))
  fisher <- fisher.test(test)
  genes$d29ms_p[i] <- round(fisher$p, digits=4)
  genes$d29ms_frac_slow[i] <- round(length(unique(slow[,1]))/nrow(medslow29), digits=3)
  genes$d29ms_frac_fast[i] <- round(length(unique(fast[,1]))/nrow(fast29), digits=3)
  # d79
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,8]=="positive")
  fast <- subset(tab2, tab2[,8]=="negative")
  test[1,1] <- length(unique(slow[,1]))
  test[2,1] <- length(unique(fast[,1]))
  test[1,2] <- (nrow(slow79)-length(unique(slow[,1])))
  test[2,2] <- (nrow(fast79)-length(unique(fast[,1])))
  fisher <- fisher.test(test)
  genes$d79_p[i] <- round(fisher$p, digits=4)
  genes$d79_frac_slow[i] <- round(length(unique(slow[,1]))/nrow(slow79), digits=3)
  genes$d79_frac_fast[i] <- round(length(unique(fast[,1]))/nrow(fast79), digits=3)}

# adding info to the columns describing the group with highest prevalence of alterations
for (i in 1:nrow(genes)) {
  if (isTRUE(genes$d15_frac_fast[i] < genes$d15_frac_slow[i])) {genes$d15_top_prevalence[i] <- c("slow")}
  if (isTRUE(genes$d15_frac_fast[i] > genes$d15_frac_slow[i])) {genes$d15_top_prevalence[i] <- c("fast")}
  if (isTRUE(genes$d29all_frac_fast[i] > genes$d29all_frac_slow[i] & genes$d29all_frac_fast[i] > genes$d29all_frac_med[i])) {genes$d29all_top_prevalence[i] <- c("fast")}
  if (isTRUE(genes$d29all_frac_med[i] > genes$d29all_frac_slow[i] & genes$d29all_frac_fast[i] < genes$d29all_frac_med[i])) {genes$d29all_top_prevalence[i] <- c("intermediate")}
  if (isTRUE(genes$d29all_frac_fast[i] < genes$d29all_frac_slow[i] & genes$d29all_frac_slow[i] > genes$d29all_frac_med[i])) {genes$d29all_top_prevalence[i] <- c("slow")}
  if (isTRUE(genes$d29fm_frac_fast[i] < genes$d29fm_frac_slow[i])) {genes$d29fm_top_prevalence[i] <- c("slow")}
  if (isTRUE(genes$d29fm_frac_fast[i] > genes$d29fm_frac_slow[i])) {genes$d29fm_top_prevalence[i] <- c("fast")}
  if (isTRUE(genes$d29ms_frac_fast[i] < genes$d29ms_frac_slow[i])) {genes$d29ms_top_prevalence[i] <- c("slow")}
  if (isTRUE(genes$d29ms_frac_fast[i] > genes$d29ms_frac_slow[i])) {genes$d29ms_top_prevalence[i] <- c("fast")}
  if (isTRUE(genes$d79_frac_fast[i] < genes$d79_frac_slow[i])) {genes$d79_top_prevalence[i] <- c("slow")}
  if (isTRUE(genes$d79_frac_fast[i] > genes$d79_frac_slow[i])) {genes$d79_top_prevalence[i] <- c("fast")}}

# adding info about the numbers in each dataset
for (i in 1:nrow(genes)) {
  tab2 <- subset(tab, tab[,3]=="AMP" & tab[,4]==genes[i,1])
  wgs <- subset(tab2, tab2$cohort=="WGS")
  array <- subset(tab2, tab2$cohort=="array")
  panel <- subset(tab2, tab2$cohort=="panel")
  genes$n_WGS[i] <- length(unique(wgs[,1]))
  genes$n_array[i] <- length(unique(array[,1]))
  genes$n_panel[i] <- length(unique(panel[,1]))}

# removing genes with 0 hits in any cohort, or less than 5 hits altogether
remove <- vector()
for (i in 1:nrow(genes)) {
  if (genes$n_WGS[i]==0 | genes$n_array[i]==0 | genes$n_panel[i]==0 | sum(genes[i,24:26]) <= 5) {
    remove <- append(remove, i)}}
amp <- genes[-remove,]

# writing the table out
write.table(amp, "amp_fishers_test_results.tab", sep="\t", quote=F, col.names=T, row.names=F)

# checking the numbers of significant genes in each timepoint comparison
# d15
del15 <- subset(del, as.numeric(del$d15_p) < 0.1)
amp15 <- subset(amp, as.numeric(amp$d15_p) < 0.1)
genes15 <- append(del15[,1], amp15[,1])
genes15 <- unique(genes15)
# d29
del29 <- subset(del, as.numeric(del$d29fm_p) < 0.1 | as.numeric(del$d29ms_p) < 0.1)
amp29 <- subset(amp, as.numeric(amp$d29fm_p) < 0.1 | as.numeric(amp$d29ms_p) < 0.1)
genes29 <- append(del29[,1], amp29[,1])
genes29 <- unique(genes29)
# d79 
del79 <- subset(del, as.numeric(del$d79_p) < 0.1)
amp79 <- subset(amp, as.numeric(amp$d79_p) < 0.1)
genes79 <- append(del79[,1], amp79[,1])
genes79 <- unique(genes79)

## combining RNAseq info with the CN results
library(Homo.sapiens)
library(stats)
library(sjmisc)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(limma)
library(Glimma)
library(edgeR)

# preparing the RNAseq data first
# the table with the response group info for both cohorts
resp_info <- read.table("combined_response_groups.txt", header=T)

# read counts
counts_matrix <- read.table("ALL_GEX_read_counts.txt", header=T)
rownames(counts_matrix) <- counts_matrix[,1]
counts_matrix2 <- read.table("full_salmon_gene_counts.tsv", header=T)
rownames(counts_matrix2) <- counts_matrix2[,1]

# making sure that the genes and row order are identical between the
identical(rownames(counts_matrix), rownames(counts_matrix2))
# which is TRUE

# saving the sample IDs as vectors
uusamples <- colnames(counts_matrix)
gesamples <- colnames(counts_matrix2)

# combining the two tables
counts_matrix <- cbind(counts_matrix, counts_matrix2)

# generating the DGE object from the count matrix
data <- DGEList(counts_matrix)

# organizing sample information, adding info on the sequencing batch, sex, and response to therapy 
samples <- rownames(data$samples)
batch <- vector()
uuinfo <- read.csv("ALL_GEX.csv")
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
# table describing the sexes
phenotype <- uuinfo
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
# adding response info
d15_response <- vector()
d29_response <- vector()
d79_response <- vector()
for (s in samples) {
  resp2 <- subset(resp_info, resp_info[,1]==s)
  if (is.na(resp2[1,2])) {d15_response <- append(d15_response, "unknown")}
  if (is.na(resp2[1,3])) {d29_response <- append(d29_response, "unknown")}
  if (is.na(resp2[1,4])) {d79_response <- append(d79_response, "unknown")}
  if (!is.na(resp2[1,2])) {d15_response <- append(d15_response, resp2[1,2])}
  if (!is.na(resp2[1,3])) {d29_response <- append(d29_response, resp2[1,3])}
  if (!is.na(resp2[1,4])) {d79_response <- append(d79_response, resp2[1,4])}}
data$samples$d15_response <- d15_response
data$samples$d29_response <- d29_response
data$samples$d79_response <- d79_response
# adding two additional columns for d29 comparisons: mrd pos vs neg, and fast+intermediate vs slow
d29_response2 <- vector()
d29_response3 <- vector()
for (s in samples) {
  resp2 <- subset(resp_info, resp_info[,1]==s)
  if (is.na(resp2[1,3])) {d29_response2 <- append(d29_response2, "unknown")}
  if (is.na(resp2[1,3])) {d29_response3 <- append(d29_response3, "unknown")}
  if (isTRUE(resp2[1,3]=="fast")) {d29_response2 <- append(d29_response2, "mrd_neg")}
  if (isTRUE(resp2[1,3]=="slow" | resp2[1,3]=="intermediate")) {d29_response2 <- append(d29_response2, "mrd_pos")}
  if (isTRUE(resp2[1,3]=="fast" | resp2[1,3]=="intermediate")) {d29_response3 <- append(d29_response3, "fastmed")}
  if (isTRUE(resp2[1,3]=="slow")) {d29_response3 <- append(d29_response3, "slow")}}
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
# computing the log-CPM values
lcpm <- cpm(data, log=T)
# checking the mean and median of library sizes
L <- mean(data$samples$lib.size) * 1e-6
M <- median(data$samples$lib.size) * 1e-6
c(L,M)

# removing genes lowly expressed across the response group comparisons
table(rowSums(data$counts==0)==50)
keep.exprs1 <- filterByExpr(data, group=d15_response)
keep.exprs2 <- filterByExpr(data, group=d29_response)
keep.exprs3 <- filterByExpr(data, group=d79_response)
for (i in 1:length(keep.exprs1)) {
  if (isFALSE(keep.exprs1[i]) & isTRUE(keep.exprs2[i]) | isFALSE(keep.exprs1[i]) & isTRUE(keep.exprs3[i])) {keep.exprs1[i] <- TRUE}
}
data <- data[keep.exprs1,,keep.lib.sizes=FALSE]

# calculating the normalization factors
data <- calcNormFactors(data, method="TMM")
# and generating the lcpm data
lcpm <- cpm(data, log=T)
# saving the samples as a vector
samples <- colnames(lcpm)

# preparing the CNV data for the comparisons

# WGS CNV results
nc <- read.table("/Volumes/groups/allseq/students/sanni_wrk/ETV6-RUNX1_cases/balsamic_results/data_files/lncRNA_genes/CNVs_affecting_nc_genes.tab")
cnv <- read.table("/Volumes/groups/allseq/students/sanni_wrk/ETV6-RUNX1_cases/balsamic_results/data_files/GEPARD_CNV/annotations/collective_coding_CNVs.tab", header=T)
# combining the tables
cnv <- cbind(cnv[,c(1,2,5,7)])
nc <- nc[,1:4]
colnames(cnv) <- colnames(nc)
getab <- rbind(cnv, nc)
getab[,3] <- gsub("DUP", "AMP", getab[,3])
getab <- subset(getab, getab[,3]!="LOH")

# methylation & panelseq CNV results
cnv2 <- read.table("collective_CNV_filt_annot_per_gene.tab", header=T)
# subsetting to gain/loss/coding seq consequences
cnv2 <- subset(cnv2, cnv2[,6]=="loss" | cnv2[,6]=="nc_gene_loss" | cnv2[,6]=="coding_sequence_variant" | cnv2[,6]=="gain" | cnv2[,6]=="nc_gene_gain" | cnv2[,6]=="nc_gene_sequence_variant")
# ordering the table into the same order as the WGS data
cnv2 <- cnv2[,c(15,1,5,7)]
# combining the tables
colnames(getab) <- colnames(cnv2)
cnv <- rbind(getab, cnv2)

# subsetting the CNV data to only cases with available RNAseq data
cnv <- subset(cnv, cnv[,1] %in% samples)

# deletions
dels <- subset(cnv, cnv[,3]=="DEL" & cnv[,4] %in% del[,1])
# and adding columns that describe the geneID of the genes
geneid <- matrix(ncol=2, nrow=nrow(dels))
genes <- as.vector(dels[,4])
kts <- keytypes(Homo.sapiens)
genetab <- AnnotationDbi::select(Homo.sapiens, keys=genes, columns=c("GENEID", "ENSEMBL"), keytype="SYMBOL")
genetab <- genetab[!duplicated(genetab$SYMBOL),]
for (i in 1:nrow(dels)) {
  genetab2 <- subset(genetab, genetab[,1]==dels[i,4])
  geneid[i,1] <- genetab2[1,2]
  geneid[i,2] <- genetab2[1,3]}
dels <- cbind(dels, geneid)

# amplifications
amps <- subset(cnv, cnv[,3]=="AMP" & cnv[,4] %in% amp[,1])
# and adding columns that describe the geneID of the genes
geneid <- matrix(ncol=2, nrow=nrow(amps))
genes <- as.vector(amps[,4])
kts <- keytypes(Homo.sapiens)
genetab <- AnnotationDbi::select(Homo.sapiens, keys=genes, columns=c("GENEID", "ENSEMBL"), keytype="SYMBOL")
genetab <- genetab[!duplicated(genetab$SYMBOL),]
for (i in 1:nrow(amps)) {
  genetab2 <- subset(genetab, genetab[,1]==amps[i,4])
  geneid[i,1] <- genetab2[1,2]
  geneid[i,2] <- genetab2[1,3]}
amps <- cbind(amps, geneid)

# starting with day 29 Fisher's test results for deletions
d29 <- del29
# adding columns that describe the geneID of the genes
geneid <- matrix(ncol=2, nrow=nrow(d29))
genes <- as.vector(d29[,1])
kts <- keytypes(Homo.sapiens)
genetab <- AnnotationDbi::select(Homo.sapiens, keys=genes, columns=c("GENEID", "ENSEMBL"), keytype="SYMBOL")
genetab <- genetab[!duplicated(genetab$SYMBOL),]
for (i in 1:nrow(d29)) {
  genetab2 <- subset(genetab, genetab[,1]==d29[i,1])
  geneid[i,1] <- genetab2[1,2]
  geneid[i,2] <- genetab2[1,3]} 
d29 <- cbind(d29, geneid)

# adding a columns to the table which summarizes the expression in cases with a deletion vs not
# in cases with deletions vs not
newcols <- matrix(ncol=2, nrow=nrow(d29))
for (i in 1:nrow(d29)) {
  del2 <- subset(dels, dels[,5]==d29[i,27])
  affected <- unique(del2[,1])
  unaffected <- setdiff(samples, affected)
  if (isTRUE(d29[i,27] %in% rownames(lcpm))) {
    exp1 <- lcpm[d29[i,27],affected]
    exp2 <- lcpm[d29[i,27],unaffected]
    if (isFALSE(is.na(exp2[1])) & length(samples)!=length(unaffected)) {
      newcols[i,1] <- mean(as.numeric(exp1))
      newcols[i,2] <- mean(as.numeric(exp2))}}}
colnames(newcols) <- c("average_expression_affected", "average_expression_unaffected")
d29 <- cbind(d29, newcols)

# filtering out the genes not expressed, and ones for which deleted cases express higher levels
d29 <- subset(d29, !is.na(d29[,29]))
remove <- vector()
for (i in 1:nrow(d29)) {
  if (isTRUE(as.numeric(d29[i,29])>as.numeric(d29[i,30]))) {remove <- append(remove, i)}}
d29 <- d29[-remove,]
del29 <- d29
# also generating a table into which all the genes of interest will be combined into
goi <- d29[,c(1,27)]

# continuing with day 15 Fisher's test results for deletions
d15 <- del15
# adding columns that describe the geneID of the genes
geneid <- matrix(ncol=2, nrow=nrow(d15))
genes <- as.vector(d15[,1])
kts <- keytypes(Homo.sapiens)
genetab <- AnnotationDbi::select(Homo.sapiens, keys=genes, columns=c("GENEID", "ENSEMBL"), keytype="SYMBOL")
genetab <- genetab[!duplicated(genetab$SYMBOL),]
for (i in 1:nrow(d15)) {
  genetab2 <- subset(genetab, genetab[,1]==d15[i,1])
  geneid[i,1] <- genetab2[1,2]
  geneid[i,2] <- genetab2[1,3]} 
d15 <- cbind(d15, geneid)

# adding columns to the table which summarizes the expression in cases with a deletion vs not
newcols <- matrix(ncol=2, nrow=nrow(d15))
for (i in 1:nrow(d15)) {
  del2 <- subset(dels, dels[,5]==d15[i,27])
  affected <- unique(del2[,1])
  unaffected <- setdiff(samples, affected)
  if (isTRUE(d15[i,27] %in% rownames(lcpm))) {
    exp1 <- lcpm[d15[i,27],affected]
    exp2 <- lcpm[d15[i,27],unaffected]
    if (isFALSE(is.na(exp2[1])) & length(samples)!=length(unaffected)) {
      newcols[i,1] <- mean(as.numeric(exp1))
      newcols[i,2] <- mean(as.numeric(exp2))}}}
colnames(newcols) <- c("average_expression_affected", "average_expression_unaffected")
d15 <- cbind(d15, newcols)

# filtering out the genes not expressed, and ones for which deleted cases express higher levels
d15 <- subset(d15, !is.na(d15[,29]))
remove <- vector()
for (i in 1:nrow(d15)) {
  if (isTRUE(as.numeric(d15[i,29])>as.numeric(d15[i,30]))) {remove <- append(remove, i)}}
d15 <- d15[-remove,]
del15 <- d15
goi <- rbind(goi, d15[,c(1,27)])

# continuing with day 79 Fisher's test results for deletions
d79 <- del79
# and adding columns that describe the geneID of the genes
geneid <- matrix(ncol=2, nrow=nrow(d79))
genes <- as.vector(d79[,1])
kts <- keytypes(Homo.sapiens)
genetab <- AnnotationDbi::select(Homo.sapiens, keys=genes, columns=c("GENEID", "ENSEMBL"), keytype="SYMBOL")
genetab <- genetab[!duplicated(genetab$SYMBOL),]
for (i in 1:nrow(d79)) {
  genetab2 <- subset(genetab, genetab[,1]==d79[i,1])
  geneid[i,1] <- genetab2[1,2]
  geneid[i,2] <- genetab2[1,3]} 
d79 <- cbind(d79, geneid)

# adding columns to the table which summarizes the expression in cases with a deletion vs not
newcols <- matrix(ncol=2, nrow=nrow(d79))
for (i in 1:nrow(d79)) {
  del2 <- subset(dels, dels[,5]==d79[i,27])
  affected <- unique(del2[,1])
  unaffected <- setdiff(samples, affected)
  if (isTRUE(d79[i,27] %in% rownames(lcpm))) {
    exp1 <- lcpm[d79[i,27],affected]
    exp2 <- lcpm[d79[i,27],unaffected]
    if (isFALSE(is.na(exp2[1])) & length(samples)!=length(unaffected)) {
      newcols[i,1] <- mean(as.numeric(exp1))
      newcols[i,2] <- mean(as.numeric(exp2))}}}
colnames(newcols) <- c("average_expression_affected", "average_expression_unaffected")
d79 <- cbind(d79, newcols)

# filtering out the genes not expressed, and ones for which deleted cases express higher levels
d79 <- subset(d79, !is.na(d79[,29]))
remove <- vector()
for (i in 1:nrow(d79)) {
  if (isTRUE(as.numeric(d79[i,29])>as.numeric(d79[i,30]))) {remove <- append(remove, i)}}
d79 <- d79[-remove,]
del79 <- d79
goi <- rbind(goi, d79[,c(1,27)])

# continuing with day 29 Fisher's test results for amplifications
d29 <- amp29
# adding columns that describe the geneID of the genes
geneid <- matrix(ncol=2, nrow=nrow(d29))
genes <- as.vector(d29[,1])
kts <- keytypes(Homo.sapiens)
genetab <- AnnotationDbi::select(Homo.sapiens, keys=genes, columns=c("GENEID", "ENSEMBL"), keytype="SYMBOL")
genetab <- genetab[!duplicated(genetab$SYMBOL),]
for (i in 1:nrow(d29)) {
  genetab2 <- subset(genetab, genetab[,1]==d29[i,1])
  geneid[i,1] <- genetab2[1,2]
  geneid[i,2] <- genetab2[1,3]} 
d29 <- cbind(d29, geneid)

# adding columns to the table which summarizes the expression in cases with a amplification vs not
newcols <- matrix(ncol=2, nrow=nrow(d29))
for (i in 1:nrow(d29)) {
  amp2 <- subset(amps, amps[,5]==d29[i,27])
  affected <- unique(amp2[,1])
  unaffected <- setdiff(samples, affected)
  if (isTRUE(d29[i,27] %in% rownames(lcpm))) {
    exp1 <- lcpm[d29[i,27],affected]
    exp2 <- lcpm[d29[i,27],unaffected]
    if (isFALSE(is.na(exp2[1])) & length(samples)!=length(unaffected)) {
      newcols[i,1] <- mean(as.numeric(exp1))
      newcols[i,2] <- mean(as.numeric(exp2))}}}
colnames(newcols) <- c("average_expression_affected", "average_expression_unaffected")
d29 <- cbind(d29, newcols)

# filtering out the genes not expressed, and ones for which ampeted cases express higher levels
d29 <- subset(d29, !is.na(d29[,29]))
remove <- vector()
for (i in 1:nrow(d29)) {
  if (isTRUE(as.numeric(d29[i,29])<as.numeric(d29[i,30]))) {remove <- append(remove, i)}}
d29 <- d29[-remove,]
amp29 <- d29
goi <- rbind(goi, d29[,c(1,27)])

# continuing with day 15 Fisher's test results for amplifications
d15 <- amp15
# adding columns that describe the geneID of the genes
geneid <- matrix(ncol=2, nrow=nrow(d15))
genes <- as.vector(d15[,1])
kts <- keytypes(Homo.sapiens)
genetab <- AnnotationDbi::select(Homo.sapiens, keys=genes, columns=c("GENEID", "ENSEMBL"), keytype="SYMBOL")
genetab <- genetab[!duplicated(genetab$SYMBOL),]
for (i in 1:nrow(d15)) {
  genetab2 <- subset(genetab, genetab[,1]==d15[i,1])
  geneid[i,1] <- genetab2[1,2]
  geneid[i,2] <- genetab2[1,3]} 
d15 <- cbind(d15, geneid)

# adding a columns to the table which summarizes the expression in cases with a amplification vs not
newcols <- matrix(ncol=2, nrow=nrow(d15))
for (i in 1:nrow(d15)) {
  amp2 <- subset(amps, amps[,5]==d15[i,27])
  affected <- unique(amp2[,1])
  unaffected <- setdiff(samples, affected)
  if (isTRUE(d15[i,27] %in% rownames(lcpm))) {
    exp1 <- lcpm[d15[i,27],affected]
    exp2 <- lcpm[d15[i,27],unaffected]
    if (isFALSE(is.na(exp2[1])) & length(samples)!=length(unaffected)) {
      newcols[i,1] <- mean(as.numeric(exp1))
      newcols[i,2] <- mean(as.numeric(exp2))}}}
colnames(newcols) <- c("average_expression_affected", "average_expression_unaffected")
d15 <- cbind(d15, newcols)

# filtering out the genes not expressed, and ones for which ampeted cases express higher levels
d15 <- subset(d15, !is.na(d15[,29]))
remove <- vector()
for (i in 1:nrow(d15)) {
  if (isTRUE(as.numeric(d15[i,29])<as.numeric(d15[i,30]))) {remove <- append(remove, i)}}
d15 <- d15[-remove,]
amp15 <- d15
goi <- rbind(goi, d15[,c(1,27)])

# continuing with day 79 Fisher's test results for amplifications
d79 <- amp79
# adding columns that describe the geneID of the genes
geneid <- matrix(ncol=2, nrow=nrow(d79))
genes <- as.vector(d79[,1])
kts <- keytypes(Homo.sapiens)
genetab <- AnnotationDbi::select(Homo.sapiens, keys=genes, columns=c("GENEID", "ENSEMBL"), keytype="SYMBOL")
genetab <- genetab[!duplicated(genetab$SYMBOL),]
for (i in 1:nrow(d79)) {
  genetab2 <- subset(genetab, genetab[,1]==d79[i,1])
  geneid[i,1] <- genetab2[1,2]
  geneid[i,2] <- genetab2[1,3]} 
d79 <- cbind(d79, geneid)

# adding a columns to the table which summarizes the expression in cases with a amplification vs not
newcols <- matrix(ncol=2, nrow=nrow(d79))
for (i in 1:nrow(d79)) {
  amp2 <- subset(amps, amps[,5]==d79[i,27])
  affected <- unique(amp2[,1])
  unaffected <- setdiff(samples, affected)
  if (isTRUE(d79[i,27] %in% rownames(lcpm))) {
    exp1 <- lcpm[d79[i,27],affected]
    exp2 <- lcpm[d79[i,27],unaffected]
    if (isFALSE(is.na(exp2[1])) & length(samples)!=length(unaffected)) {
      newcols[i,1] <- mean(as.numeric(exp1))
      newcols[i,2] <- mean(as.numeric(exp2))}}}
colnames(newcols) <- c("average_expression_affected", "average_expression_unaffected")
d79 <- cbind(d79, newcols)

# filtering out the genes not expressed, and ones for which ampeted cases express higher levels
d79 <- subset(d79, !is.na(d79[,29]))
remove <- vector()
for (i in 1:nrow(d79)) {
  if (isTRUE(as.numeric(d79[i,29])<as.numeric(d79[i,30]))) {remove <- append(remove, i)}}
d79 <- d79[-remove,]
amp79 <- d79
goi <- rbind(goi, d79[,c(1,27)])

## simplifying the genes of interest table to only include unique entries
goi <- unique(goi)

# copy number DE analyses for the genes of interest
newcols <- matrix(ncol=5, nrow=nrow(goi))
colnames(newcols) <- c("ave_exp", "lfc_amp", "lfc_del", "pval", "adj.pval")
for (g in 1:nrow(goi)) {
  
  # determining the copy number status and conducting DE analyses accordingly
  dels2 <- subset(dels, dels[,4]==goi[g,1])
  amps2 <- subset(amps, amps[,4]==goi[g,1])
  cn_status <- vector()
  del_cases <- unique(dels2[,1])
  amp_cases <- unique(amps2[,1])
  neutral_cases <- setdiff(samples, del_cases)
  neutral_cases <- setdiff(neutral_cases, amp_cases)
  for (s in samples) {
    if (s %in% del_cases & !s %in% amp_cases) {cn_status <- append(cn_status, "del")}
    if (s %in% amp_cases & !s %in% del_cases) {cn_status <- append(cn_status, "amp")}
    if (s %in% del_cases & s %in% amp_cases) {
      cn_status <- append(cn_status, "aa_neutral")
      amp_cases <- amp_cases[-which(amp_cases==s)]
      del_cases <- del_cases[-which(del_cases==s)]}
    if (s %in% neutral_cases) {cn_status <- append(cn_status, "aa_neutral")}}
  data$samples$cn_status <- cn_status
  
  # DE analysis adjusting for batch
  design <- model.matrix(~cn_status+batch)
  v <- voom(data, design, plot=TRUE)
  vfit <- lmFit(v, design)
  efit <- eBayes(vfit)
  if (length(del_cases)==0) {
    toptab <- topTable(efit, coef="cn_statusamp", adjust.method = "BH", n=Inf)
    write.table(toptab, paste0("cn_efits/efit_", goi[g,1], ".tab"), sep="\t", quote=F, col.names=T, row.names = F)}
  if (length(amp_cases)==0) {
    toptab <- topTable(efit, coef="cn_statusdel", adjust.method = "BH", n=Inf)
    write.table(toptab, paste0("cn_efits/efit_", goi[g,1], ".tab"), sep="\t", quote=F, col.names=T, row.names = F)}
  if (length(amp_cases)!=0 & length(del_cases)!=0) {
    toptab <- topTable(efit, coef=c("cn_statusamp", "cn_statusdel"), adjust.method = "BH", n=Inf)
    write.table(toptab, paste0("cn_efits/efit_", goi[g,1], ".tab"), sep="\t", quote=F, col.names=T, row.names = F)}
  # adding the info for the specific gene in the goi table
  fit <- read.table(paste0("cn_efits/efit_", goi[g,1], ".tab"), header=T)
  fit <- subset(fit, fit$ENSEMBL==goi[g,2])
  if (colnames(fit)[9]=="adj.P.Val") {
    newcols[g,1] <- fit[1,6]
    newcols[g,2] <- fit[1,4]
    newcols[g,3] <- fit[1,5]
    newcols[g,4] <- fit[1,8]
    newcols[g,5] <- fit[1,9]}
  if (colnames(fit)[9]!="adj.P.Val" & length(amp_cases)==0) {
    newcols[g,1] <- fit[1,5]
    newcols[g,2] <- "."
    newcols[g,3] <- fit[1,4]
    newcols[g,4] <- fit[1,7]
    newcols[g,5] <- fit[1,8]}
  if (colnames(fit)[9]!="adj.P.Val" & length(del_cases)==0) {
    newcols[g,1] <- fit[1,5]
    newcols[g,2] <- fit[1,4]
    newcols[g,3] <- "."
    newcols[g,4] <- fit[1,7]
    newcols[g,5] <- fit[1,8]}}
goi <- cbind(goi, newcols)

# adding the annotations to the deletion/amplification results
rownames(goi) <- goi[,1]
goi2 <- goi[del29[,1],]
del29 <- cbind(del29, goi2[,3:7])
goi2 <- goi[amp29[,1],]
amp29 <- cbind(amp29, goi2[,3:7])
goi2 <- goi[del15[,1],]
del15 <- cbind(del15, goi2[,3:7])
goi2 <- goi[amp15[,1],]
amp15 <- cbind(amp15, goi2[,3:7])
goi2 <- goi[del79[,1],]
del79 <- cbind(del79, goi2[,3:7])
goi2 <- goi[amp79[,1],]
amp79 <- cbind(amp79, goi2[,3:7])
# and writing the tables out
write.table(del29, "d29_del_fishers_test_results_combined_de_info.tab", sep="\t", quote=F, col.names=T, row.names=F)
write.table(amp29, "d29_amp_fishers_test_results_combined_de_info.tab", sep="\t", quote=F, col.names=T, row.names=F)
write.table(del15, "d15_del_fishers_test_results_combined_de_info.tab", sep="\t", quote=F, col.names=T, row.names=F)
write.table(amp15, "d15_amp_fishers_test_results_combined_de_info.tab", sep="\t", quote=F, col.names=T, row.names=F)
write.table(del79, "d79_del_fishers_test_results_combined_de_info.tab", sep="\t", quote=F, col.names=T, row.names=F)
write.table(amp79, "d79_amp_fishers_test_results_combined_de_info.tab", sep="\t", quote=F, col.names=T, row.names=F)
