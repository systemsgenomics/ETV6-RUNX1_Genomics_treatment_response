# process the cBioPortal output and conduct Fisher's exact tests on functionally filtered SNV/indel results
library(sjmisc)

# reading in the mutationmapper output
mut <- read.table("collective_cbioportal_output.tab", sep="\t", header=T)

# reformatting the results so that the annotations and functional impact predictions are in their own columns
newcols <- matrix(ncol=8, nrow=nrow(mut))
for (i in 1:nrow(mut)) {
  vec1 <- unlist(strsplit(mut[i,4], ";"))
  vec2 <- unlist(strsplit(mut[i,5], ";"))
  newcols[i,1] <- vec1[1]
  newcols[i,2] <- vec1[2]
  newcols[i,3] <- vec1[3]
  newcols[i,4] <- vec1[4]
  newcols[i,5] <- vec1[5]
  newcols[i,6] <- vec2[1]
  newcols[i,7] <- vec2[2]
  newcols[i,8] <- vec2[3]}
mut <- cbind(mut[,c(1,3,6,7,11,12,13,14,15,26)], newcols[,c(1,4,6:8)])

# reading in the original input table
input <- read.table("cbioportal_input.txt", header=T)
# adding a column describing the affected gene
gene <- matrix(ncol=1, nrow=nrow(mut))
for (i in 1:nrow(mut)) {
  input2 <- subset(input, input[,1]==mut[i,5] & input[,2]==mut[i,6]  & input[,4]==mut[i,8]  & input[,5]==mut[i,9] & input[,8]==mut[i,1])
  gene[i,1] <- input2[1,7]}
mut <- cbind(mut[,1:4], gene, mut[,5:15])

# writing the table out
colnames(mut) <- c("caseid", "aa_change", "mut_type", "var_type", "gene", "chr", "start", "end", "ref", "alt", "clinvar", "oncokb", "hotspot", "mutationassessor", "sift", "polyphen-2")
write.table(mut, "func_annotated_snv_indel_results.tab", col.names = T, row.names=F, quote=F, sep="\t")

# then filtering the list to only functionally relevant
keep <- vector()
for (i in 1:nrow(mut)) {
  if (isTRUE(str_contains(mut[i,14], "high") | str_contains(mut[i,15], "deleterious") | str_contains(mut[i,16], "damaging") | mut[i,3]=="Frame_Shift_Del" | mut[i,3]=="Nonsense_Mutation" | mut[i,3]=="Frame_Shift_Ins" | mut[i,3]=="Splice_Site")) {keep <- append(keep, i)}}
keep <- mut[keep,]
# writing the table out
write.table(keep, "func_filtered_snv_indel_results.tab", col.names = T, row.names=F, quote=F, sep="\t")
tab <- keep

# reading in the response info
resp <- read.table("sample_responses.txt", header=T)
# generating tables of the samples in each responder group
fast15 <- subset(resp, resp$d15_response=="fast")
slow15 <- subset(resp, resp$d15_response=="slow")
fast29 <- subset(resp, resp$d29_response=="fast")
med29 <- subset(resp, resp$d29_response=="intermediate")
slow29 <- subset(resp, resp$d29_response=="slow")
medslow29 <- subset(resp, resp$d29_response=="slow" | resp$d29_response=="intermediate")
fastmed29 <- subset(resp, resp$d29_response=="fast" | resp$d29_response=="intermediate")
fast79 <- subset(resp, resp$d79_mrd_status=="negative")
slow79 <- subset(resp, resp$d79_mrd_status=="positive")

# adding columns to the table describing the response groups and cohort
newcols <- matrix(ncol=4, nrow=nrow(tab))
for (i in 1:nrow(tab)) {
  resp2 <- subset(resp, resp[,1]==tab[i,1])
  newcols[i,1] <- resp2[1,2]
  newcols[i,2] <- resp2[1,3]
  newcols[i,3] <- resp2[1,4]
  if (isTRUE(str_contains(tab[i,1], "GE"))) {
    newcols[i,4] <- "WGS"}
  if (isFALSE(str_contains(tab[i,1], "GE"))) {
    newcols[i,4] <- "panel"}}
colnames(newcols) <- c("d15_response", "d29_response", "d79_response", "cohort")
tab <- cbind(tab, newcols)

# determining the unique genes
genes <- unique(tab$gene)

# generating a table into which the results can be collected
newcols <- matrix(ncol=23, nrow=length(genes))
colnames(newcols) <- c("d15_p", "d29all_p", "d29fm_p", "d29ms_p", "d79_p", "d15_top_prevalence", "d29all_top_prevalence", "d29fm_top_prevalence", "d29ms_top_prevalence", "d79_top_prevalence", "d15_frac_fast", "d15_frac_slow", "d29all_frac_fast", "d29all_frac_med", "d29all_frac_slow", "d29fm_frac_fast", "d29fm_frac_slow", "d29ms_frac_fast", "d29ms_frac_slow", "d79_frac_fast", "d79_frac_slow", "n_WGS", "n_panel")
genes <- cbind(genes, newcols)
rownames(genes) <- genes[,1]
genes <- as.data.frame(genes)
for (i in 1:nrow(genes)) {
  tab2 <- subset(tab, tab[,5]==genes[i,1])
  # d15
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,17]=="slow")
  fast <- subset(tab2, tab2[,17]=="fast")
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
  slow <- subset(tab2, tab2[,18]=="slow")
  med <- subset(tab2, tab2[,18]=="intermediate")
  fast <- subset(tab2, tab2[,18]=="fast")
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
  slow <- subset(tab2, tab2[,18]=="slow")
  fast <- subset(tab2, tab2[,18]=="fast" | tab2[,18]=="intermediate")
  test[1,1] <- length(unique(slow[,1]))
  test[2,1] <- length(unique(fast[,1]))
  test[1,2] <- (nrow(slow29)-length(unique(slow[,1])))
  test[2,2] <- (nrow(fastmed29)-length(unique(fast[,1])))
  fisher <- fisher.test(test)
  genes$d29fm_p[i] <- round(fisher$p, digits=4)
  genes$d29fm_frac_slow[i] <- round(length(unique(slow[,1]))/nrow(slow29), digits=3)
  genes$d29fm_frac_fast[i] <- round(length(unique(fast[,1]))/nrow(fastmed29), digits=3)
  # d29 fast vs intermediate and slow
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,18]=="slow" | tab2[,18]=="intermediate")
  fast <- subset(tab2, tab2[,18]=="fast")
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
  slow <- subset(tab2, tab2[,19]=="positive")
  fast <- subset(tab2, tab2[,19]=="negative")
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
  tab2 <- subset(tab, tab[,5]==genes[i,1])
  wgs <- subset(tab2, tab2$cohort=="WGS")
  panel <- subset(tab2, tab2$cohort=="panel")
  genes$n_WGS[i] <- length(unique(wgs[,1]))
  genes$n_panel[i] <- length(unique(panel[,1]))}
# excluding variants with less than 3 mutations altogether
genes <- as.matrix(genes)
genes <- subset(genes, (as.numeric(genes[,23])+as.numeric(genes[,24])) > 2)

# writing the table out
write.table(genes, "func_filt_panel_wgs_snv_fishers_test_results.tab", col.names = T, row.names = F, sep="\t", quote=F)

# distinguishing the more significant genes using the p-value cutoff of 0.05
genes <- as.matrix(genes)
genes <- subset(genes, as.numeric(genes[,2]) < 0.05 | as.numeric(genes[,3]) < 0.05 | as.numeric(genes[,4]) < 0.05 | as.numeric(genes[,5]) < 0.05 | as.numeric(genes[,6]) < 0.05)

# writing the table out
write.table(genes, "sig_func_filt_panel_wgs_snv_fishers_test_results.tab", col.names = T, row.names = F, sep="\t", quote=F)
