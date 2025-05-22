# compare the CNV results to the drug screen genes
library(ggplot2)
library(sjmisc)

# response table
resp <- read.table("sample_responses_cnv_cohort.txt", header=T)

# WGS CNV results
cnv <- read.table("/Volumes/groups/allseq/students/sanni_wrk/ETV6-RUNX1_cases/balsamic_results/data_files/GEPARD_CNV/annotations/collective_CNV_annotated_only_full.tab", header=F)
# simplifying the table
getab <- cnv[,c(2:5,1)]
# panelseq and methylation array results
tab <- read.table("collective_CNV_filt_annot_per_event.tab", header=T, sep="\t")
# simplifying the table
tab <- tab[,c(1:3,5,8)]
# combining the tables
colnames(getab) <- colnames(tab)
tab <- rbind(tab, getab)
# subsetting it to only include the deletions
tab <- subset(tab, tab$CNV_type=="DEL")

# adding response columns to the table
newcols <- matrix(ncol=4, nrow=nrow(tab))
for (i in 1:nrow(tab)) {
  resp2 <- subset(resp, resp[,1]==tab[i,5])
  newcols[i,1] <- resp2[1,2]
  newcols[i,2] <- resp2[1,3]
  newcols[i,3] <- resp2[1,4]
  if (isTRUE(str_contains(tab[i,5], "ALL_"))) {
    newcols[i,4] <- "array"}
  if (isTRUE(str_contains(tab[i,5], "GE"))) {
    newcols[i,4] <- "WGS"}
  if (isFALSE(str_contains(tab[i,5], "GE")) & isFALSE(str_contains(tab[i,5], "ALL_"))) {
    newcols[i,4] <- "panel"}}
colnames(newcols) <- c("d15_response", "d29_response", "d79_response", "cohort")
tab <- cbind(tab, newcols)

# reading in the top sensitizing and resistance inducing CRISPR screen results of the different drugs
vcr <- read.table("vincristine.csv", header=T)
vcr_res <- subset(vcr, vcr$pos.fdr <= 0.25)
vcr_sens <- subset(vcr, vcr$neg.fdr <= 0.25)
mp <- read.table("6mercaptopurine.csv", header=T)
mp_res <- subset(mp, mp$pos.fdr <= 0.25)
mp_sens <- subset(mp, mp$neg.fdr <= 0.25)
arac <- read.table("cytarabine.csv", header=T)
arac_res <- subset(arac, arac$pos.fdr <= 0.25)
arac_sens <- subset(arac, arac$neg.fdr <= 0.25)
dnr <- read.table("daunorubicin.csv", header=T)
dnr_res <- subset(dnr, dnr$pos.fdr <= 0.25)
dnr_sens <- subset(dnr, dnr$neg.fdr <= 0.25)
lasp <- read.table("Lasparaginase.csv", header=T)
lasp_res <- subset(lasp, lasp$pos.fdr <= 0.25)
lasp_sens <- subset(lasp, lasp$neg.fdr <= 0.25)
maf <- read.table("maphosamide.csv", header=T)
maf_res <- subset(maf, maf$pos.fdr <= 0.25)
maf_sens <- subset(maf, maf$neg.fdr <= 0.25)
mtx <- read.table("methotrexate.csv", header=T)
mtx_res <- subset(mtx, mtx$pos.fdr <= 0.25)
mtx_sens <- subset(mtx, mtx$neg.fdr <= 0.25)
dex <- read.csv("dexamethasone_sensitivity_significant.csv", header=T)
dex_res <- subset(dex, dex$Rho.Phenotype > 0)
dex_sens <- subset(dex, dex$Rho.Phenotype < 0)

# all deletions vs all CRISPR FDR < 0.25 genes

# first, determining the coordinates for all the genes included in the CRISPR analyses
genes <- vcr[,1]
original_genes <- genes
genetab <- matrix(ncol=3, nrow=length(genes))
# fixing the miRNA gene names
for (g in 1:length(genes)) {
  genes[g] <- gsub("hsa-mir-", "MIR", genes[g])}
for (g in 1:length(genes)) {
  genes[g] <- toupper(genes[g])}
for (g in 1:length(genes)) {
  genes[g] <- gsub("ORF", "orf", genes[g])}
rownames(genetab) <- genes

# determining the gene coordinates based on the same RefSeq file as used in AnnotSV
annotsv <- read.table("genes.RefSeq.sorted.bed")
for (g in genes) {
  annot2 <- subset(annotsv, annotsv[,5]==g)
  if (nrow(annot2)!=0) {
    genetab[g,1] <- annot2[1,1]
    genetab[g,2] <- annot2[1,2]
    genetab[g,3] <- annot2[1,3]}}

# substituting the rownames with the original gene names that match the CRISPR data
rownames(genetab) <- original_genes
genetab <- genetab[-remove,]
colnames(genetab) <- c("chr", "start", "end")

# generating a panel of the genes of interest
# i.e., genes passing the FDR cutoff
genes <- append(vcr_res[,1], vcr_sens[,1])
genes <- append(genes, mp_res[,1])
genes <- append(genes, mp_sens[,1])
genes <- append(genes, arac_res[,1])
genes <- append(genes, arac_sens[,1])
genes <- append(genes, dnr_res[,1])
genes <- append(genes, dnr_sens[,1])
genes <- append(genes, lasp_res[,1])
genes <- append(genes, lasp_sens[,1])
genes <- append(genes, maf_res[,1])
genes <- append(genes, maf_sens[,1])
genes <- append(genes, mtx_res[,1])
genes <- append(genes, mtx_sens[,1])
genes <- append(genes, dex_res[,2])
genes <- append(genes, dex_sens[,2])
genes <- unique(genes)
genes <- genes[which(genes %in% rownames(genetab))]
panel <- genetab[genes,]

# generating a table for the Fisher's test results from the analyses comparing numbers of cases with any gene hits vs not
fisher15 <- matrix(ncol=6, nrow=8)
colnames(fisher15) <- c("pval_top100_res", "pval_top100_sens", "pval_top200_res", "pval_top200_sens", "pval_fdr_res", "pval_fdr_sens")
rownames(fisher15) <- c("VCR", "6MP", "ARAC", "DNR", "LASP", "MAF", "MTX", "DEX")
fisher29fs <- matrix(ncol=6, nrow=8)
colnames(fisher29fs) <- c("pval_top100_res", "pval_top100_sens", "pval_top200_res", "pval_top200_sens", "pval_fdr_res", "pval_fdr_sens")
rownames(fisher29fs) <- c("VCR", "6MP", "ARAC", "DNR", "LASP", "MAF", "MTX", "DEX")
fisher29fm <- matrix(ncol=6, nrow=8)
colnames(fisher29fm) <- c("pval_top100_res", "pval_top100_sens", "pval_top200_res", "pval_top200_sens", "pval_fdr_res", "pval_fdr_sens")
rownames(fisher29fm) <- c("VCR", "6MP", "ARAC", "DNR", "LASP", "MAF", "MTX", "DEX")
fisher29ms <- matrix(ncol=6, nrow=8)
colnames(fisher29ms) <- c("pval_top100_res", "pval_top100_sens", "pval_top200_res", "pval_top200_sens", "pval_fdr_res", "pval_fdr_sens")
rownames(fisher29ms) <- c("VCR", "6MP", "ARAC", "DNR", "LASP", "MAF", "MTX", "DEX")
fisher79 <- matrix(ncol=6, nrow=8)
colnames(fisher79) <- c("pval_top100_res", "pval_top100_sens", "pval_top200_res", "pval_top200_sens", "pval_fdr_res", "pval_fdr_sens")
rownames(fisher79) <- c("VCR", "6MP", "ARAC", "DNR", "LASP", "MAF", "MTX", "DEX")

# generating a table into which the info about the screen genes can be collected
delmat <- matrix(ncol=10, nrow=1)
colnames(delmat) <- c("chrom", "start", "end", "cnv_type", "sample_id", "d15_response", "d29_response", "d79_response", "cohort", "gene")
for (g in genes) {
  # checking from the del table which cases have deletions overlapping the gene region
  del2 <- subset(del, del[,1]==panel[g,1])
  loss <- vector()
  for (i in 1:nrow(del2)) {
    if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
    if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
    if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,2])<=as.numeric(panel[g,3]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
    if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}}
  del2 <- del2[loss,]
  gene <- matrix(ncol=1, nrow=nrow(del2))
  gene[,1] <- g
  loss <- cbind(del2, gene)
  colnames(loss) <- colnames(delmat)
  delmat <- rbind(delmat, loss)}
delmat <- delmat[-1,]

# generating a table for plotting the numbers of hits
samples <- resp[,1]
caseid <- vector()
d15_response <- vector()
d29_response <- vector()
d79_response <- vector()
impact <- rep(c("res", "sens"), 2864)
drug <- rep(c("VCR", "VCR", "6MP", "6MP", "ARAC", "ARAC", "DNR", "DNR", "LASP", "LASP", "MAF", "MAF", "MTX", "MTX", "DEX", "DEX"), 358)
n_del <- vector()
for (s in samples) {
  caseid <- append(caseid, rep(s, 16))
  resp2 <- subset(resp, resp[,1]==s)
  d15_response <- append(d15_response, rep(resp2[1,2], 16))
  d29_response <- append(d29_response, rep(resp2[1,3], 16))
  d79_response <- append(d79_response, rep(resp2[1,4], 16))
  del <- subset(delmat, delmat[,5]==s)
  hits <- subset(del, del[,10] %in% vcr_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% vcr_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mp_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mp_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% arac_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% arac_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% dnr_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% dnr_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% lasp_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% lasp_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% maf_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% maf_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mtx_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mtx_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% dex_top[,2])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% dex_bottom[,2])
  n_del <- append(n_del, length(unique(hits[,10])))}
plot <- data.frame(caseid, d15_response, d29_response, d79_response, impact, drug, n_del)
res <- subset(plot, plot$impact=="res") 
sens <- subset(plot, plot$impact=="sens")

# generating summary tables for the medians and means of hits, and wilcox values
drugs <- unique(drug)
summary15 <- matrix(ncol=12, nrow=length(drugs))
rownames(summary15) <- drugs
colnames(summary15) <- c("total_res_genes", "total_sens_genes", "wilcox_pval_res", "wilcox_pval_sens", "res_del_median_fast", "res_del_median_slow", "res_del_mean_fast", "res_del_mean_slow", "sens_del_median_fast", "sens_del_median_slow", "sens_del_mean_fast", "sens_del_mean_slow")
summary15[1,1] <- nrow(vcr_top)
summary15[1,2] <- (nrow(vcr_bottom)-1)
summary15[2,1] <- nrow(mp_top)
summary15[2,2] <- nrow(mp_bottom)
summary15[3,1] <- nrow(arac_top)
summary15[3,2] <- nrow(arac_bottom)
summary15[4,1] <- nrow(dnr_top)
summary15[4,2] <- nrow(dnr_bottom)
summary15[5,1] <- nrow(lasp_top)
summary15[5,2] <- nrow(lasp_bottom)
summary15[6,1] <- nrow(maf_top)
summary15[6,2] <- nrow(maf_bottom)
summary15[7,1] <- (nrow(mtx_top)-3)
summary15[7,2] <- nrow(mtx_bottom)
summary15[8,1] <- (nrow(dex_top)-1)
summary15[8,2] <- (nrow(dex_bottom)-1)
summary29 <- matrix(ncol=24, nrow=length(drugs))
rownames(summary29) <- drugs
colnames(summary29) <- c("total_res_genes", "total_sens_genes", "wilcox_pval_res_f_s",  "wilcox_pval_res_fm_s",  "wilcox_pval_res_f_ms", "wilcox_pval_sens_f_s",  "wilcox_pval_sens_fm_s", "wilcox_pval_sens_f_ms",  "res_del_median_fast", "res_del_median_fastmed", "res_del_median_medslow", "res_del_median_slow", "res_del_mean_fast", "res_del_mean_fastmed", "res_del_mean_medslow", "res_del_mean_slow", "sens_del_median_fast", "sens_del_median_fastmed", "sens_del_median_medslow", "sens_del_median_slow", "sens_del_mean_fast", "sens_del_mean_fastmed", "sens_del_mean_medslow", "sens_del_mean_slow")
summary29[1,1] <- nrow(vcr_top)
summary29[1,2] <- (nrow(vcr_bottom)-1)
summary29[2,1] <- nrow(mp_top)
summary29[2,2] <- nrow(mp_bottom)
summary29[3,1] <- nrow(arac_top)
summary29[3,2] <- nrow(arac_bottom)
summary29[4,1] <- nrow(dnr_top)
summary29[4,2] <- nrow(dnr_bottom)
summary29[5,1] <- nrow(lasp_top)
summary29[5,2] <- nrow(lasp_bottom)
summary29[6,1] <- nrow(maf_top)
summary29[6,2] <- nrow(maf_bottom)
summary29[7,1] <- (nrow(mtx_top)-3)
summary29[7,2] <- nrow(mtx_bottom)
summary29[8,1] <- (nrow(dex_top)-1)
summary29[8,2] <- (nrow(dex_bottom)-1)
summary79 <- matrix(ncol=12, nrow=length(drugs))
rownames(summary79) <- drugs
colnames(summary79) <- c("total_res_genes", "total_sens_genes", "wilcox_pval_res", "wilcox_pval_sens", "res_del_median_fast", "res_del_median_slow", "res_del_mean_fast", "res_del_mean_slow", "sens_del_median_fast", "sens_del_median_slow", "sens_del_mean_fast", "sens_del_mean_slow")
summary79[1,1] <- nrow(vcr_top)
summary79[1,2] <- (nrow(vcr_bottom)-1)
summary79[2,1] <- nrow(mp_top)
summary79[2,2] <- nrow(mp_bottom)
summary79[3,1] <- nrow(arac_top)
summary79[3,2] <- nrow(arac_bottom)
summary79[4,1] <- nrow(dnr_top)
summary79[4,2] <- nrow(dnr_bottom)
summary79[5,1] <- nrow(lasp_top)
summary79[5,2] <- nrow(lasp_bottom)
summary79[6,1] <- nrow(maf_top)
summary79[6,2] <- nrow(maf_bottom)
summary79[7,1] <- (nrow(mtx_top)-3)
summary79[7,2] <- nrow(mtx_bottom)
summary79[8,1] <- (nrow(dex_top)-1)
summary79[8,2] <- (nrow(dex_bottom)-1)

# drawing the plots
# resistance hits
plot15 <- subset(res, !is.na(res$d15_response) & res$d15_response!="unknown")
pdf("crispr_res_hit_boxplot_d15.pdf")
ggplot(plot15, aes(x = drug, y = n_del, col = d15_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
plot29 <- subset(res, !is.na(res$d29_response) & res$d29_response!="unknown")
pdf("crispr_res_hit_boxplot_d29.pdf")
ggplot(plot29, aes(x = drug, y = n_del, col = d29_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#f1a340", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
plot79 <- subset(res, !is.na(res$d79_response) & res$d79_response!="unknown")
pdf("crispr_res_hit_boxplot_d79.pdf")
ggplot(plot79, aes(x = drug, y = n_del, col = d79_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

# collecting the info to the summary tables as well
for (d in drugs) {
  fast <- subset(plot15, plot15$d15_response=="fast" & plot15$drug==d)
  slow <- subset(plot15, plot15$d15_response=="slow" & plot15$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary15[d,3] <- round(test$p.value, digits=3)
  summary15[d,5] <- median(fast$n_del)
  summary15[d,6] <- median(slow$n_del)
  summary15[d,7] <- round(mean(fast$n_del), digits=3)
  summary15[d,8] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d)
  med <- subset(plot29, plot29$d29_response=="intermediate" & plot29$drug==d)
  slow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d)
  fastmed <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  medslow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  test2 <- wilcox.test(fastmed$n_del, slow$n_del)
  test3 <- wilcox.test(fast$n_del, medslow$n_del)
  summary29[d,3] <- round(test$p.value, digits=3)
  summary29[d,4] <- round(test2$p.value, digits=3)
  summary29[d,5] <- round(test3$p.value, digits=3)
  summary29[d,9] <- median(fast$n_del)
  summary29[d,10] <- median(fastmed$n_del)
  summary29[d,11] <- median(medslow$n_del)
  summary29[d,12] <- median(slow$n_del)
  summary29[d,13] <- round(mean(fast$n_del), digits=3)
  summary29[d,14] <- round(mean(fastmed$n_del), digits=3)
  summary29[d,15] <- round(mean(medslow$n_del), digits=3)
  summary29[d,16] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot79, plot79$d79_response=="negative" & plot79$drug==d)
  slow <- subset(plot79, plot79$d79_response=="positive" & plot79$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary79[d,3] <- round(test$p.value, digits=3)
  summary79[d,5] <- median(fast$n_del)
  summary79[d,6] <- median(slow$n_del)
  summary79[d,7] <- round(mean(fast$n_del), digits=3)
  summary79[d,8] <- round(mean(slow$n_del), digits=3)}
# conducting the Fisher's tests comparing the response groups
for (d in drugs) {
  n_fast_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del==0)
  n_fast_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del > 0)
  n_slow_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del==0)
  n_slow_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher15[d,5] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fs[d,5] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fm[d,5] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29ms[d,5] <- test$p.value
  n_fast_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del==0)
  n_fast_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del > 0)
  n_slow_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del==0)
  n_slow_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher79[d,5] <- test$p.value}

# sensitivity hits
plot15 <- subset(sens, !is.na(sens$d15_response) & sens$d15_response!="unknown")
pdf("crispr_sens_hit_boxplot_d15.pdf")
ggplot(plot15, aes(x = drug, y = n_del, col = d15_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
plot29 <- subset(sens, !is.na(sens$d29_response) & sens$d29_response!="unknown")
pdf("crispr_sens_hit_boxplot_d29.pdf")
ggplot(plot29, aes(x = drug, y = n_del, col = d29_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#f1a340", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
plot79 <- subset(sens, !is.na(sens$d79_response) & sens$d79_response!="unknown")
pdf("crispr_sens_hit_boxplot_d79.pdf")
ggplot(plot79, aes(x = drug, y = n_del, col = d79_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

# collecting the info to the summary tables as well
for (d in drugs) {
  fast <- subset(plot15, plot15$d15_response=="fast" & plot15$drug==d)
  slow <- subset(plot15, plot15$d15_response=="slow" & plot15$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary15[d,4] <- round(test$p.value, digits=3)
  summary15[d,9] <- median(fast$n_del)
  summary15[d,10] <- median(slow$n_del)
  summary15[d,11] <- round(mean(fast$n_del), digits=3)
  summary15[d,12] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d)
  med <- subset(plot29, plot29$d29_response=="intermediate" & plot29$drug==d)
  slow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d)
  fastmed <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  medslow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  test2 <- wilcox.test(fastmed$n_del, slow$n_del)
  test3 <- wilcox.test(fast$n_del, medslow$n_del)
  summary29[d,6] <- round(test$p.value, digits=3)
  summary29[d,7] <- round(test2$p.value, digits=3)
  summary29[d,8] <- round(test3$p.value, digits=3)
  summary29[d,17] <- median(fast$n_del)
  summary29[d,18] <- median(fastmed$n_del)
  summary29[d,19] <- median(medslow$n_del)
  summary29[d,20] <- median(slow$n_del)
  summary29[d,21] <- round(mean(fast$n_del), digits=3)
  summary29[d,22] <- round(mean(fastmed$n_del), digits=3)
  summary29[d,23] <- round(mean(medslow$n_del), digits=3)
  summary29[d,24] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot79, plot79$d79_response=="negative" & plot79$drug==d)
  slow <- subset(plot79, plot79$d79_response=="positive" & plot79$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary79[d,4] <- round(test$p.value, digits=3)
  summary79[d,9] <- median(fast$n_del)
  summary79[d,10] <- median(slow$n_del)
  summary79[d,11] <- round(mean(fast$n_del), digits=3)
  summary79[d,12] <- round(mean(slow$n_del), digits=3)}

# conducting the Fisher's tests comparing the response groups
for (d in drugs) {
  n_fast_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del==0)
  n_fast_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del > 0)
  n_slow_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del==0)
  n_slow_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher15[d,6] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fs[d,6] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fm[d,6] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29ms[d,6] <- test$p.value
  n_fast_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del==0)
  n_fast_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del > 0)
  n_slow_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del==0)
  n_slow_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher79[d,6] <- test$p.value}

# writing the summary tables out
summary15 <- cbind(drugs, summary15)
write.table(summary15, "crispr_hit_summary_d15.tab", sep="\t", quote=F, col.names=T, row.names=F)
summary29 <- cbind(drugs, summary29)
write.table(summary29, "crispr_hit_summary_d29.tab", sep="\t", quote=F, col.names=T, row.names=F)
summary79 <- cbind(drugs, summary79)
write.table(summary79, "crispr_hit_summary_d79.tab", sep="\t", quote=F, col.names=T, row.names=F)

## using the top 200 genes for each drug
vcr_res <- subset(vcr, vcr$pos.rank <= 200)
vcr_sens <- subset(vcr, vcr$neg.rank > 20266)
dnr_res <- subset(dnr, dnr$pos.rank <= 200)
dnr_sens <- subset(dnr, dnr$neg.rank > 20266)
mtx_res <- subset(mtx, mtx$pos.rank <= 200)
mtx_sens <- subset(mtx, mtx$neg.rank > 20266)
maf_res <- subset(maf, maf$pos.rank <= 200)
maf_sens <- subset(maf, maf$neg.rank > 20265)
mp_res <- subset(mp, mp$pos.rank <= 200)
mp_sens <- subset(mp, mp$neg.rank > 20266)
lasp_res <- subset(lasp, lasp$pos.rank <= 200)
lasp_sens <- subset(lasp, lasp$neg.rank > 20266)
arac_res <- subset(arac, arac$pos.rank <= 200)
arac_sens <- subset(arac, arac$neg.rank > 20266)
genes <- append(vcr_res[,1], vcr_sens[,1])
genes <- append(genes, mp_res[,1])
genes <- append(genes, mp_sens[,1])
genes <- append(genes, arac_res[,1])
genes <- append(genes, arac_sens[,1])
genes <- append(genes, dnr_res[,1])
genes <- append(genes, dnr_sens[,1])
genes <- append(genes, lasp_res[,1])
genes <- append(genes, lasp_sens[,1])
genes <- append(genes, maf_res[,1])
genes <- append(genes, maf_sens[,1])
genes <- append(genes, mtx_res[,1])
genes <- append(genes, mtx_sens[,1])
genes <- unique(genes)
genes <- genes[which(genes %in% rownames(genetab))]
panel <- genetab[genes,]

# generating a table into which the info about the CRISPR genes can be collected
delmat <- matrix(ncol=10, nrow=1)
colnames(delmat) <- c("chrom", "start", "end", "cnv_type", "sample_id", "d15_response", "d29_response", "d79_response", "cohort", "gene")
for (g in genes) {
  # checking from the del table, which cases have deletions overlapping the gene region
  del2 <- subset(del, del[,1]==panel[g,1])
  loss <- vector()
  for (i in 1:nrow(del2)) {
    if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
    if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
    if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,2])<=as.numeric(panel[g,3]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
    if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}}
  del2 <- del2[loss,]
  gene <- matrix(ncol=1, nrow=nrow(del2))
  gene[,1] <- g
  loss <- cbind(del2, gene)
  colnames(loss) <- colnames(delmat)
  delmat <- rbind(delmat, loss)}
delmat <- delmat[-1,]

# generating a table for plotting the numbers of hits
samples <- resp[,1]
caseid <- vector()
d15_response <- vector()
d29_response <- vector()
d79_response <- vector()
impact <- rep(c("res", "sens"), 2506)
drug <- rep(c("VCR", "VCR", "6MP", "6MP", "ARAC", "ARAC", "DNR", "DNR", "LASP", "LASP", "MAF", "MAF", "MTX", "MTX"), 358)
n_del <- vector()
for (s in samples) {
  caseid <- append(caseid, rep(s, 14))
  resp2 <- subset(resp, resp[,1]==s)
  d15_response <- append(d15_response, rep(resp2[1,2], 14))
  d29_response <- append(d29_response, rep(resp2[1,3], 14))
  d79_response <- append(d79_response, rep(resp2[1,4], 14))
  del <- subset(delmat, delmat[,5]==s)
  hits <- subset(del, del[,10] %in% vcr_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% vcr_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mp_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mp_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% arac_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% arac_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% dnr_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% dnr_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% lasp_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% lasp_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% maf_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% maf_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mtx_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mtx_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))}
plot <- data.frame(caseid, d15_response, d29_response, d79_response, impact, drug, n_del)
res <- subset(plot, plot$impact=="res") 
sens <- subset(plot, plot$impact=="sens")

# generating summary tables for the medians and means of hits, and wilcox values
drugs <- unique(drug)
summary15 <- matrix(ncol=12, nrow=length(drugs))
rownames(summary15) <- drugs
colnames(summary15) <- c("total_res_genes", "total_sens_genes", "wilcox_pval_res", "wilcox_pval_sens", "res_del_median_fast", "res_del_median_slow", "res_del_mean_fast", "res_del_mean_slow", "sens_del_median_fast", "sens_del_median_slow", "sens_del_mean_fast", "sens_del_mean_slow")
summary15[1,1] <- nrow(vcr_top)
summary15[1,2] <- (nrow(vcr_bottom)-1)
summary15[2,1] <- nrow(mp_top)
summary15[2,2] <- nrow(mp_bottom)
summary15[3,1] <- nrow(arac_top)
summary15[3,2] <- nrow(arac_bottom)
summary15[4,1] <- nrow(dnr_top)
summary15[4,2] <- nrow(dnr_bottom)
summary15[5,1] <- nrow(lasp_top)
summary15[5,2] <- nrow(lasp_bottom)
summary15[6,1] <- nrow(maf_top)
summary15[6,2] <- nrow(maf_bottom)
summary15[7,1] <- (nrow(mtx_top)-3)
summary15[7,2] <- nrow(mtx_bottom)
summary29 <- matrix(ncol=24, nrow=length(drugs))
rownames(summary29) <- drugs
colnames(summary29) <- c("total_res_genes", "total_sens_genes", "wilcox_pval_res_f_s",  "wilcox_pval_res_fm_s",  "wilcox_pval_res_f_ms", "wilcox_pval_sens_f_s",  "wilcox_pval_sens_fm_s", "wilcox_pval_sens_f_ms",  "res_del_median_fast", "res_del_median_fastmed", "res_del_median_medslow", "res_del_median_slow", "res_del_mean_fast", "res_del_mean_fastmed", "res_del_mean_medslow", "res_del_mean_slow", "sens_del_median_fast", "sens_del_median_fastmed", "sens_del_median_medslow", "sens_del_median_slow", "sens_del_mean_fast", "sens_del_mean_fastmed", "sens_del_mean_medslow", "sens_del_mean_slow")
summary29[1,1] <- nrow(vcr_top)
summary29[1,2] <- (nrow(vcr_bottom)-1)
summary29[2,1] <- nrow(mp_top)
summary29[2,2] <- nrow(mp_bottom)
summary29[3,1] <- nrow(arac_top)
summary29[3,2] <- nrow(arac_bottom)
summary29[4,1] <- nrow(dnr_top)
summary29[4,2] <- nrow(dnr_bottom)
summary29[5,1] <- nrow(lasp_top)
summary29[5,2] <- nrow(lasp_bottom)
summary29[6,1] <- nrow(maf_top)
summary29[6,2] <- nrow(maf_bottom)
summary29[7,1] <- (nrow(mtx_top)-3)
summary29[7,2] <- nrow(mtx_bottom)
summary79 <- matrix(ncol=12, nrow=length(drugs))
rownames(summary79) <- drugs
colnames(summary79) <- c("total_res_genes", "total_sens_genes", "wilcox_pval_res", "wilcox_pval_sens", "res_del_median_fast", "res_del_median_slow", "res_del_mean_fast", "res_del_mean_slow", "sens_del_median_fast", "sens_del_median_slow", "sens_del_mean_fast", "sens_del_mean_slow")
summary79[1,1] <- nrow(vcr_top)
summary79[1,2] <- (nrow(vcr_bottom)-1)
summary79[2,1] <- nrow(mp_top)
summary79[2,2] <- nrow(mp_bottom)
summary79[3,1] <- nrow(arac_top)
summary79[3,2] <- nrow(arac_bottom)
summary79[4,1] <- nrow(dnr_top)
summary79[4,2] <- nrow(dnr_bottom)
summary79[5,1] <- nrow(lasp_top)
summary79[5,2] <- nrow(lasp_bottom)
summary79[6,1] <- nrow(maf_top)
summary79[6,2] <- nrow(maf_bottom)
summary79[7,1] <- (nrow(mtx_top)-3)
summary79[7,2] <- nrow(mtx_bottom)

# resistance hits
plot15 <- subset(res, !is.na(res$d15_response) & res$d15_response!="unknown")
pdf("crispr_top200_res_hit_boxplot_d15.pdf")
ggplot(plot15, aes(x = drug, y = n_del, col = d15_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()
plot29 <- subset(res, !is.na(res$d29_response) & res$d29_response!="unknown")
pdf("crispr_top200_res_hit_boxplot_d29.pdf")
ggplot(plot29, aes(x = drug, y = n_del, col = d29_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#f1a340", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()
plot79 <- subset(res, !is.na(res$d79_response) & res$d79_response!="unknown")
pdf("crispr_top200_res_hit_boxplot_d79.pdf")
ggplot(plot79, aes(x = drug, y = n_del, col = d79_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()

# collecting the info to the summary tables as well
for (d in drugs) {
  fast <- subset(plot15, plot15$d15_response=="fast" & plot15$drug==d)
  slow <- subset(plot15, plot15$d15_response=="slow" & plot15$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary15[d,3] <- round(test$p.value, digits=3)
  summary15[d,5] <- median(fast$n_del)
  summary15[d,6] <- median(slow$n_del)
  summary15[d,7] <- round(mean(fast$n_del), digits=3)
  summary15[d,8] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d)
  med <- subset(plot29, plot29$d29_response=="intermediate" & plot29$drug==d)
  slow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d)
  fastmed <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  medslow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  test2 <- wilcox.test(fastmed$n_del, slow$n_del)
  test3 <- wilcox.test(fast$n_del, medslow$n_del)
  summary29[d,3] <- round(test$p.value, digits=3)
  summary29[d,4] <- round(test2$p.value, digits=3)
  summary29[d,5] <- round(test3$p.value, digits=3)
  summary29[d,9] <- median(fast$n_del)
  summary29[d,10] <- median(fastmed$n_del)
  summary29[d,11] <- median(medslow$n_del)
  summary29[d,12] <- median(slow$n_del)
  summary29[d,13] <- round(mean(fast$n_del), digits=3)
  summary29[d,14] <- round(mean(fastmed$n_del), digits=3)
  summary29[d,15] <- round(mean(medslow$n_del), digits=3)
  summary29[d,16] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot79, plot79$d79_response=="negative" & plot79$drug==d)
  slow <- subset(plot79, plot79$d79_response=="positive" & plot79$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary79[d,3] <- round(test$p.value, digits=3)
  summary79[d,5] <- median(fast$n_del)
  summary79[d,6] <- median(slow$n_del)
  summary79[d,7] <- round(mean(fast$n_del), digits=3)
  summary79[d,8] <- round(mean(slow$n_del), digits=3)}

# conducting the Fisher's tests comparing the response groups
for (d in drugs) {
  n_fast_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del==0)
  n_fast_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del > 0)
  n_slow_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del==0)
  n_slow_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher15[d,3] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fs[d,3] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fm[d,3] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29ms[d,3] <- test$p.value
  n_fast_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del==0)
  n_fast_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del > 0)
  n_slow_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del==0)
  n_slow_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher79[d,3] <- test$p.value}

# sensitivity hits
plot15 <- subset(sens, !is.na(sens$d15_response) & sens$d15_response!="unknown")
pdf("crispr_top200_sens_hit_boxplot_d15.pdf")
ggplot(plot15, aes(x = drug, y = n_del, col = d15_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()
plot29 <- subset(sens, !is.na(sens$d29_response) & sens$d29_response!="unknown")
pdf("crispr_top200_sens_hit_boxplot_d29.pdf")
ggplot(plot29, aes(x = drug, y = n_del, col = d29_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#f1a340", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()
plot79 <- subset(sens, !is.na(sens$d79_response) & sens$d79_response!="unknown")
pdf("crispr_top200_sens_hit_boxplot_d79.pdf")
ggplot(plot79, aes(x = drug, y = n_del, col = d79_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()

# collecting the info to the summary tables as well
for (d in drugs) {
  fast <- subset(plot15, plot15$d15_response=="fast" & plot15$drug==d)
  slow <- subset(plot15, plot15$d15_response=="slow" & plot15$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary15[d,4] <- round(test$p.value, digits=3)
  summary15[d,9] <- median(fast$n_del)
  summary15[d,10] <- median(slow$n_del)
  summary15[d,11] <- round(mean(fast$n_del), digits=3)
  summary15[d,12] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d)
  med <- subset(plot29, plot29$d29_response=="intermediate" & plot29$drug==d)
  slow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d)
  fastmed <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  medslow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  test2 <- wilcox.test(fastmed$n_del, slow$n_del)
  test3 <- wilcox.test(fast$n_del, medslow$n_del)
  summary29[d,6] <- round(test$p.value, digits=3)
  summary29[d,7] <- round(test2$p.value, digits=3)
  summary29[d,8] <- round(test3$p.value, digits=3)
  summary29[d,17] <- median(fast$n_del)
  summary29[d,18] <- median(fastmed$n_del)
  summary29[d,19] <- median(medslow$n_del)
  summary29[d,20] <- median(slow$n_del)
  summary29[d,21] <- round(mean(fast$n_del), digits=3)
  summary29[d,22] <- round(mean(fastmed$n_del), digits=3)
  summary29[d,23] <- round(mean(medslow$n_del), digits=3)
  summary29[d,24] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot79, plot79$d79_response=="negative" & plot79$drug==d)
  slow <- subset(plot79, plot79$d79_response=="positive" & plot79$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary79[d,4] <- round(test$p.value, digits=3)
  summary79[d,9] <- median(fast$n_del)
  summary79[d,10] <- median(slow$n_del)
  summary79[d,11] <- round(mean(fast$n_del), digits=3)
  summary79[d,12] <- round(mean(slow$n_del), digits=3)}

# conducting the Fisher's tests comparing the response groups
for (d in drugs) {
  n_fast_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del==0)
  n_fast_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del > 0)
  n_slow_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del==0)
  n_slow_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher15[d,4] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fs[d,4] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fm[d,4] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29ms[d,4] <- test$p.value
  n_fast_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del==0)
  n_fast_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del > 0)
  n_slow_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del==0)
  n_slow_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher79[d,4] <- test$p.value}

# writing the summary tables out
summary15 <- cbind(drugs, summary15[,3:12])
write.table(summary15, "crispr_top200_hit_summary_d15.tab", sep="\t", quote=F, col.names=T, row.names=F)
summary29 <- cbind(drugs, summary29[,3:20])
write.table(summary29, "crispr_top200_hit_summary_d29.tab", sep="\t", quote=F, col.names=T, row.names=F)
summary79 <- cbind(drugs, summary79[,3:12])
write.table(summary79, "crispr_top200_hit_summary_d79.tab", sep="\t", quote=F, col.names=T, row.names=F)

## top 100 genes for each drug
vcr_res <- subset(vcr, vcr$pos.rank <= 100)
vcr_sens <- subset(vcr, vcr$neg.rank > 20366)
dnr_res <- subset(dnr, dnr$pos.rank <= 100)
dnr_sens <- subset(dnr, dnr$neg.rank > 20366)
mtx_res <- subset(mtx, mtx$pos.rank <= 100)
mtx_sens <- subset(mtx, mtx$neg.rank > 20366)
maf_res <- subset(maf, maf$pos.rank <= 100)
maf_sens <- subset(maf, maf$neg.rank > 20365)
mp_res <- subset(mp, mp$pos.rank <= 100)
mp_sens <- subset(mp, mp$neg.rank > 20366)
lasp_res <- subset(lasp, lasp$pos.rank <= 100)
lasp_sens <- subset(lasp, lasp$neg.rank > 20366)
arac_res <- subset(arac, arac$pos.rank <= 100)
arac_sens <- subset(arac, arac$neg.rank > 20366)
genes <- append(vcr_res[,1], vcr_sens[,1])
genes <- append(genes, mp_res[,1])
genes <- append(genes, mp_sens[,1])
genes <- append(genes, arac_res[,1])
genes <- append(genes, arac_sens[,1])
genes <- append(genes, dnr_res[,1])
genes <- append(genes, dnr_sens[,1])
genes <- append(genes, lasp_res[,1])
genes <- append(genes, lasp_sens[,1])
genes <- append(genes, maf_res[,1])
genes <- append(genes, maf_sens[,1])
genes <- append(genes, mtx_res[,1])
genes <- append(genes, mtx_sens[,1])
genes <- unique(genes)
genes <- genes[which(genes %in% rownames(genetab))]
panel <- genetab[genes,]

# generating a table into which the info about the CRISPR genes can be collected
delmat <- matrix(ncol=10, nrow=1)
colnames(delmat) <- c("chrom", "start", "end", "cnv_type", "sample_id", "d15_response", "d29_response", "d79_response", "cohort", "gene")
for (g in genes) {
  # checking from the del table, which cases have deletions overlapping the gene region
  del2 <- subset(del, del[,1]==panel[g,1])
  loss <- vector()
  for (i in 1:nrow(del2)) {
    if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
    if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
    if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,2])<=as.numeric(panel[g,3]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
    if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}}
  del2 <- del2[loss,]
  gene <- matrix(ncol=1, nrow=nrow(del2))
  gene[,1] <- g
  loss <- cbind(del2, gene)
  colnames(loss) <- colnames(delmat)
  delmat <- rbind(delmat, loss)}
delmat <- delmat[-1,]

# generating a table for plotting the numbers of hits
samples <- resp[,1]
caseid <- vector()
d15_response <- vector()
d29_response <- vector()
d79_response <- vector()
impact <- rep(c("res", "sens"), 2506)
drug <- rep(c("VCR", "VCR", "6MP", "6MP", "ARAC", "ARAC", "DNR", "DNR", "LASP", "LASP", "MAF", "MAF", "MTX", "MTX"), 358)
n_del <- vector()
for (s in samples) {
  caseid <- append(caseid, rep(s, 14))
  resp2 <- subset(resp, resp[,1]==s)
  d15_response <- append(d15_response, rep(resp2[1,2], 14))
  d29_response <- append(d29_response, rep(resp2[1,3], 14))
  d79_response <- append(d79_response, rep(resp2[1,4], 14))
  del <- subset(delmat, delmat[,5]==s)
  hits <- subset(del, del[,10] %in% vcr_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% vcr_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mp_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mp_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% arac_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% arac_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% dnr_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% dnr_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% lasp_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% lasp_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% maf_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% maf_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mtx_top[,1])
  n_del <- append(n_del, length(unique(hits[,10])))
  hits <- subset(del, del[,10] %in% mtx_bottom[,1])
  n_del <- append(n_del, length(unique(hits[,10])))}
plot <- data.frame(caseid, d15_response, d29_response, d79_response, impact, drug, n_del)
res <- subset(plot, plot$impact=="res") 
sens <- subset(plot, plot$impact=="sens")

# generating summary tables for the medians and means of hits, and wilcox values
drugs <- unique(drug)
summary15 <- matrix(ncol=12, nrow=length(drugs))
rownames(summary15) <- drugs
colnames(summary15) <- c("total_res_genes", "total_sens_genes", "wilcox_pval_res", "wilcox_pval_sens", "res_del_median_fast", "res_del_median_slow", "res_del_mean_fast", "res_del_mean_slow", "sens_del_median_fast", "sens_del_median_slow", "sens_del_mean_fast", "sens_del_mean_slow")
summary15[1,1] <- nrow(vcr_res)
summary15[1,2] <- nrow(vcr_sens)
summary15[2,1] <- nrow(mp_res)
summary15[2,2] <- nrow(mp_sens)
summary15[3,1] <- nrow(arac_res)
summary15[3,2] <- nrow(arac_sens)
summary15[4,1] <- nrow(dnr_res)
summary15[4,2] <- nrow(dnr_sens)
summary15[5,1] <- nrow(lasp_res)
summary15[5,2] <- nrow(lasp_sens)
summary15[6,1] <- nrow(maf_res)
summary15[6,2] <- nrow(maf_sens)
summary15[7,1] <- nrow(mtx_res)
summary15[7,2] <- nrow(mtx_sens)
summary29 <- matrix(ncol=24, nrow=length(drugs))
rownames(summary29) <- drugs
colnames(summary29) <- c("total_res_genes", "total_sens_genes", "wilcox_pval_res_f_s",  "wilcox_pval_res_fm_s",  "wilcox_pval_res_f_ms", "wilcox_pval_sens_f_s",  "wilcox_pval_sens_fm_s", "wilcox_pval_sens_f_ms",  "res_del_median_fast", "res_del_median_fastmed", "res_del_median_medslow", "res_del_median_slow", "res_del_mean_fast", "res_del_mean_fastmed", "res_del_mean_medslow", "res_del_mean_slow", "sens_del_median_fast", "sens_del_median_fastmed", "sens_del_median_medslow", "sens_del_median_slow", "sens_del_mean_fast", "sens_del_mean_fastmed", "sens_del_mean_medslow", "sens_del_mean_slow")
summary29[1,1] <- nrow(vcr_res)
summary29[1,2] <- nrow(vcr_sens)
summary29[2,1] <- nrow(mp_res)
summary29[2,2] <- nrow(mp_sens)
summary29[3,1] <- nrow(arac_res)
summary29[3,2] <- nrow(arac_sens)
summary29[4,1] <- nrow(dnr_res)
summary29[4,2] <- nrow(dnr_sens)
summary29[5,1] <- nrow(lasp_res)
summary29[5,2] <- nrow(lasp_sens)
summary29[6,1] <- nrow(maf_res)
summary29[6,2] <- nrow(maf_sens)
summary29[7,1] <- nrow(mtx_res)
summary29[7,2] <- nrow(mtx_sens)
summary79 <- matrix(ncol=12, nrow=length(drugs))
rownames(summary79) <- drugs
colnames(summary79) <- c("total_res_genes", "total_sens_genes", "wilcox_pval_res", "wilcox_pval_sens", "res_del_median_fast", "res_del_median_slow", "res_del_mean_fast", "res_del_mean_slow", "sens_del_median_fast", "sens_del_median_slow", "sens_del_mean_fast", "sens_del_mean_slow")
summary79[1,1] <- nrow(vcr_res)
summary79[1,2] <- nrow(vcr_sens)
summary79[2,1] <- nrow(mp_res)
summary79[2,2] <- nrow(mp_sens)
summary79[3,1] <- nrow(arac_res)
summary79[3,2] <- nrow(arac_sens)
summary79[4,1] <- nrow(dnr_res)
summary79[4,2] <- nrow(dnr_sens)
summary79[5,1] <- nrow(lasp_res)
summary79[5,2] <- nrow(lasp_sens)
summary79[6,1] <- nrow(maf_res)
summary79[6,2] <- nrow(maf_sens)
summary79[7,1] <- nrow(mtx_res)
summary79[7,2] <- nrow(mtx_sens)

# resistance hits
plot15 <- subset(res, !is.na(res$d15_response) & res$d15_response!="unknown")
pdf("crispr_top100_res_hit_boxplot_d15.pdf")
ggplot(plot15, aes(x = drug, y = n_del, col = d15_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()
plot29 <- subset(res, !is.na(res$d29_response) & res$d29_response!="unknown")
pdf("crispr_top100_res_hit_boxplot_d29.pdf")
ggplot(plot29, aes(x = drug, y = n_del, col = d29_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#f1a340", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()
plot79 <- subset(res, !is.na(res$d79_response) & res$d79_response!="unknown")
pdf("crispr_top100_res_hit_boxplot_d79.pdf")
ggplot(plot79, aes(x = drug, y = n_del, col = d79_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()

# collecting the info to the summary tables as well
for (d in drugs) {
  fast <- subset(plot15, plot15$d15_response=="fast" & plot15$drug==d)
  slow <- subset(plot15, plot15$d15_response=="slow" & plot15$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary15[d,3] <- round(test$p.value, digits=3)
  summary15[d,5] <- median(fast$n_del)
  summary15[d,6] <- median(slow$n_del)
  summary15[d,7] <- round(mean(fast$n_del), digits=3)
  summary15[d,8] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d)
  med <- subset(plot29, plot29$d29_response=="intermediate" & plot29$drug==d)
  slow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d)
  fastmed <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  medslow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  test2 <- wilcox.test(fastmed$n_del, slow$n_del)
  test3 <- wilcox.test(fast$n_del, medslow$n_del)
  summary29[d,3] <- round(test$p.value, digits=3)
  summary29[d,4] <- round(test2$p.value, digits=3)
  summary29[d,5] <- round(test3$p.value, digits=3)
  summary29[d,9] <- median(fast$n_del)
  summary29[d,10] <- median(fastmed$n_del)
  summary29[d,11] <- median(medslow$n_del)
  summary29[d,12] <- median(slow$n_del)
  summary29[d,13] <- round(mean(fast$n_del), digits=3)
  summary29[d,14] <- round(mean(fastmed$n_del), digits=3)
  summary29[d,15] <- round(mean(medslow$n_del), digits=3)
  summary29[d,16] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot79, plot79$d79_response=="negative" & plot79$drug==d)
  slow <- subset(plot79, plot79$d79_response=="positive" & plot79$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary79[d,3] <- round(test$p.value, digits=3)
  summary79[d,5] <- median(fast$n_del)
  summary79[d,6] <- median(slow$n_del)
  summary79[d,7] <- round(mean(fast$n_del), digits=3)
  summary79[d,8] <- round(mean(slow$n_del), digits=3)}

# conducting the Fisher's tests comparing the response groups
for (d in drugs) {
  n_fast_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del==0)
  n_fast_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del > 0)
  n_slow_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del==0)
  n_slow_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher15[d,1] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fs[d,1] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fm[d,1] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29ms[d,1] <- test$p.value
  n_fast_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del==0)
  n_fast_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del > 0)
  n_slow_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del==0)
  n_slow_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher79[d,1] <- test$p.value}

# sensitivity hits
plot15 <- subset(sens, !is.na(sens$d15_response) & sens$d15_response!="unknown")
pdf("crispr_top100_sens_hit_boxplot_d15.pdf")
ggplot(plot15, aes(x = drug, y = n_del, col = d15_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()
plot29 <- subset(sens, !is.na(sens$d29_response) & sens$d29_response!="unknown")
pdf("crispr_top100_sens_hit_boxplot_d29.pdf")
ggplot(plot29, aes(x = drug, y = n_del, col = d29_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#f1a340", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()
plot79 <- subset(sens, !is.na(sens$d79_response) & sens$d79_response!="unknown")
pdf("crispr_top100_sens_hit_boxplot_d79.pdf")
ggplot(plot79, aes(x = drug, y = n_del, col = d79_response)) + geom_boxplot() + scale_color_manual(values=c("#a1d76a", "#b2182b")) + ylab("Number of genes deleted") + xlab("Drug")
dev.off()

# collecting the info to the summary tables as well
for (d in drugs) {
  fast <- subset(plot15, plot15$d15_response=="fast" & plot15$drug==d)
  slow <- subset(plot15, plot15$d15_response=="slow" & plot15$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary15[d,4] <- round(test$p.value, digits=3)
  summary15[d,9] <- median(fast$n_del)
  summary15[d,10] <- median(slow$n_del)
  summary15[d,11] <- round(mean(fast$n_del), digits=3)
  summary15[d,12] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d)
  med <- subset(plot29, plot29$d29_response=="intermediate" & plot29$drug==d)
  slow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d)
  fastmed <- subset(plot29, plot29$d29_response=="fast" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  medslow <- subset(plot29, plot29$d29_response=="slow" & plot29$drug==d | plot29$d29_response=="intermediate" & plot29$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  test2 <- wilcox.test(fastmed$n_del, slow$n_del)
  test3 <- wilcox.test(fast$n_del, medslow$n_del)
  summary29[d,6] <- round(test$p.value, digits=3)
  summary29[d,7] <- round(test2$p.value, digits=3)
  summary29[d,8] <- round(test3$p.value, digits=3)
  summary29[d,17] <- median(fast$n_del)
  summary29[d,18] <- median(fastmed$n_del)
  summary29[d,19] <- median(medslow$n_del)
  summary29[d,20] <- median(slow$n_del)
  summary29[d,21] <- round(mean(fast$n_del), digits=3)
  summary29[d,22] <- round(mean(fastmed$n_del), digits=3)
  summary29[d,23] <- round(mean(medslow$n_del), digits=3)
  summary29[d,24] <- round(mean(slow$n_del), digits=3)
  fast <- subset(plot79, plot79$d79_response=="negative" & plot79$drug==d)
  slow <- subset(plot79, plot79$d79_response=="positive" & plot79$drug==d)
  test <- wilcox.test(fast$n_del, slow$n_del)
  summary79[d,4] <- round(test$p.value, digits=3)
  summary79[d,9] <- median(fast$n_del)
  summary79[d,10] <- median(slow$n_del)
  summary79[d,11] <- round(mean(fast$n_del), digits=3)
  summary79[d,12] <- round(mean(slow$n_del), digits=3)}

# conducting the Fisher's tests comparing the response groups
for (d in drugs) {
  n_fast_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del==0)
  n_fast_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="fast" & plot15$n_del > 0)
  n_slow_neg <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del==0)
  n_slow_pos <- subset(plot15, plot15$drug==d & plot15$d15_response=="slow" & plot15$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher15[d,2] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fs[d,2] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29fm[d,2] <- test$p.value
  n_fast_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del==0)
  n_fast_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="fast" & plot29$n_del > 0)
  n_slow_neg <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del==0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del==0)
  n_slow_pos <- subset(plot29, plot29$drug==d & plot29$d29_response=="slow" & plot29$n_del > 0 | plot29$drug==d & plot29$d29_response=="intermediate" & plot29$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher29ms[d,2] <- test$p.value
  n_fast_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del==0)
  n_fast_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="negative" & plot79$n_del > 0)
  n_slow_neg <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del==0)
  n_slow_pos <- subset(plot79, plot79$drug==d & plot79$d79_response=="positive" & plot79$n_del > 0)
  fortest <- matrix(ncol=2, nrow=2)
  fortest[1,1] <- nrow(n_fast_neg)
  fortest[1,2] <- nrow(n_fast_pos)
  fortest[2,1] <- nrow(n_slow_neg)
  fortest[2,2] <- nrow(n_slow_pos)
  test <- fisher.test(fortest)
  fisher79[d,2] <- test$p.value}

# writing the summary tables out
summary15 <- cbind(drugs, summary15[,3:12])
write.table(summary15, "crispr_top100_hit_summary_d15.tab", sep="\t", quote=F, col.names=T, row.names=F)
summary29 <- cbind(drugs, summary29[,3:20])
write.table(summary29, "crispr_top100_hit_summary_d29.tab", sep="\t", quote=F, col.names=T, row.names=F)
summary79 <- cbind(drugs, summary79[,3:12])
write.table(summary79, "crispr_top100_hit_summary_d79.tab", sep="\t", quote=F, col.names=T, row.names=F)

# writing out the Fisher's test results
for (i in 1:nrow(fisher15)) {
  for (j in 1:ncol(fisher15)) {
    fisher15[i,j] <- round(as.numeric(fisher15[i,j]), digits=3)}}
for (i in 1:nrow(fisher29fs)) {
  for (j in 1:ncol(fisher29fs)) {
    fisher29fs[i,j] <- round(as.numeric(fisher29fs[i,j]), digits=3)}}
for (i in 1:nrow(fisher29fm)) {
  for (j in 1:ncol(fisher29fm)) {
    fisher29fm[i,j] <- round(as.numeric(fisher29fm[i,j]), digits=3)}}
for (i in 1:nrow(fisher29ms)) {
  for (j in 1:ncol(fisher29ms)) {
    fisher29ms[i,j] <- round(as.numeric(fisher29ms[i,j]), digits=3)}}
for (i in 1:nrow(fisher79)) {
  for (j in 1:ncol(fisher79)) {
    fisher79[i,j] <- round(as.numeric(fisher79[i,j]), digits=3)}}
drugs <- rownames(fisher15)
fisher15 <- cbind(drugs, fisher15)
write.table(fisher15, "top_crispr_hits_neg_vs_pos_fishers_test_results_d15.tab", col.names=T, row.names=F, quote=F, sep="\t")
fisher29fs <- cbind(drugs, fisher29fs)
write.table(fisher29fs, "top_crispr_hits_neg_vs_pos_fishers_test_results_d29_fast_vs_slow.tab", col.names=T, row.names=F, quote=F, sep="\t")
fisher29fm <- cbind(drugs, fisher29fm)
write.table(fisher29fm, "top_crispr_hits_neg_vs_pos_fishers_test_results_d29_fastmed_vs_slow.tab", col.names=T, row.names=F, quote=F, sep="\t")
fisher29ms <- cbind(drugs, fisher29ms)
write.table(fisher29ms, "top_crispr_hits_neg_vs_pos_fishers_test_results_d29_fast_vs_medslow.tab", col.names=T, row.names=F, quote=F, sep="\t")
fisher79 <- cbind(drugs, fisher79)
write.table(fisher79, "top_crispr_hits_neg_vs_pos_fishers_test_results_d79.tab", col.names=T, row.names=F, quote=F, sep="\t")


## generating empirical p-values for the drug screen gene deletion hits

# response table
resp <- read.table("sample_responses_cnv_cohort.txt", header=T)
# WGS CNV results
cnv <- read.table("/Volumes/groups/allseq/students/sanni_wrk/ETV6-RUNX1_cases/balsamic_results/data_files/GEPARD_CNV/annotations/collective_CNV_annotated_only_full.tab", header=F)
# simplifying the table
getab <- cnv[,c(2:5,1)]
# panelseq and methylation array results
tab <- read.table("collective_CNV_filt_annot_per_event.tab", header=T, sep="\t")
# simplifying the table
tab <- tab[,c(1:3,5,8)]
# combining the tables
colnames(getab) <- colnames(tab)
tab <- rbind(tab, getab)
# subsetting it to only include the deletions
tab <- subset(tab, tab$CNV_type=="DEL")

# adding response columns to the table
newcols <- matrix(ncol=4, nrow=nrow(tab))
for (i in 1:nrow(tab)) {
  resp2 <- subset(resp, resp[,1]==tab[i,5])
  newcols[i,1] <- resp2[1,2]
  newcols[i,2] <- resp2[1,3]
  newcols[i,3] <- resp2[1,4]
  if (isTRUE(str_contains(tab[i,5], "ALL_"))) {
    newcols[i,4] <- "array"}
  if (isTRUE(str_contains(tab[i,5], "GE"))) {
    newcols[i,4] <- "WGS"}
  if (isFALSE(str_contains(tab[i,5], "GE")) & isFALSE(str_contains(tab[i,5], "ALL_"))) {
    newcols[i,4] <- "panel"}}
colnames(newcols) <- c("d15_response", "d29_response", "d79_response", "cohort")
del <- cbind(tab, newcols)

## conducting the 1000 randomization iterations, utilizing 200 random genes

# generating a table into which the results can be collected
tab <- matrix(ncol=27, nrow=1000)
colnames(tab) <- c("iteration", "res_del_median_fast15", "res_del_median_slow15", "res_del_median_diff15", "res_del_median_fast29", "res_del_median_fastmed29", "res_del_median_medslow29", "res_del_median_slow29", "res_del_median_diff29_f_s", "res_del_median_diff29_fm_s", "res_del_median_diff29_f_ms","res_del_median_fast79", "res_del_median_slow79", "res_del_median_diff79", "sens_del_median_fast15", "sens_del_median_slow15", "sens_del_median_diff15", "sens_del_median_fast29", "sens_del_median_fastmed29", "sens_del_median_medslow29", "sens_del_median_slow29", "sens_del_median_diff29_f_s", "sens_del_median_diff29_fm_s", "sens_del_median_diff29_f_ms","sens_del_median_fast79", "sens_del_median_slow79", "sens_del_median_diff79")
tab[,1] <- 1:1000
# conducting 1000 random iterations
for (j in 1:1000) {
  
  # sampling 200 random genes to be the resistance and sensitivity genes
  resgenes <- sample(vcr[,1], 200)
  sensgenes <- sample(vcr[,1], 200)
  genes <- append(resgenes, sensgenes)
  # generating a table into which the info about the CRISPR genes can be collected
  delmat <- matrix(ncol=10, nrow=1)
  colnames(delmat) <- c("chrom", "start", "end", "size", "cnv_type", "sample_id", "d15_response", "d29_response", "d79_response", "gene")
  for (g in genes) {
    # checking from the del table, which cases have deletions overlapping the gene region
    del2 <- subset(del, del[,1]==panel[g,1])
    loss <- vector()
    for (i in 1:nrow(del2)) {
      if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
      if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
      if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,2])<=as.numeric(panel[g,3]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
      if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}}
    del2 <- del2[loss,]
    gene <- matrix(ncol=1, nrow=nrow(del2))
    gene[,1] <- g
    loss <- cbind(del2, gene)
    colnames(loss) <- colnames(delmat)
    delmat <- rbind(delmat, loss)}
  delmat <- delmat[-1,]
  
  # generating a table summarizing the numbers of hits to resistance and sensitivity genes
  samples <- d29[,2]
  samples <- append(d29[,2], gepard[,1])
  d15_response <- vector()
  d29_response <- vector()
  d79_response <- vector()
  n_del_res <- vector()
  n_del_sens <- vector()
  for (s in samples) {
    if (s %in% d29[,2]) {
      resp <- subset(d15, d15[,2]==s)
      d15_response <- append(d15_response, resp[1,4])
      resp <- subset(d29, d29[,2]==s)
      d29_response <- append(d29_response, resp[1,4])
      resp <- subset(d79, d79[,2]==s)
      d79_response <- append(d79_response, resp[1,4])}
    if (s %in% gepard[,1]) {
      resp <- subset(gepard, gepard[,1]==s)
      d15_response <- append(d15_response, resp[1,2])
      d29_response <- append(d29_response, resp[1,3])
      d79_response <- append(d79_response, resp[1,4])}
    dels <- subset(delmat, delmat[,6]==s)
    hits <- subset(dels, dels[,10] %in% resgenes)
    n_del_res <- append(n_del_res, nrow(hits))
    hits <- subset(dels, dels[,10] %in% sensgenes)
    n_del_sens <- append(n_del_sens, nrow(hits))}
  data <- data.frame(samples, d15_response, d29_response, d79_response, n_del_res, n_del_sens)
  
  # checking the medians for each response group
  fast <- subset(data, data$d15_response=="fast")
  slow <- subset(data, data$d15_response=="slow")
  tab[j,2] <- median(fast$n_del_res)
  tab[j,3] <- median(slow$n_del_res)
  tab[j,4] <- (median(fast$n_del_res)-median(slow$n_del_res))
  tab[j,15] <- median(fast$n_del_sens)
  tab[j,16] <- median(slow$n_del_sens)
  tab[j,17] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
  fast <- subset(data, data$d29_response=="fast")
  slow <- subset(data, data$d29_response=="slow")
  tab[j,5] <- median(fast$n_del_res)
  tab[j,8] <- median(slow$n_del_res)
  tab[j,9] <- (median(fast$n_del_res)-median(slow$n_del_res))
  tab[j,18] <- median(fast$n_del_sens)
  tab[j,21] <- median(slow$n_del_sens)
  tab[j,22] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
  fast <- subset(data, data$d29_response=="fast" | data$d29_response=="intermediate")
  slow <- subset(data, data$d29_response=="slow")
  tab[j,6] <- median(fast$n_del_res)
  tab[j,10] <- (median(fast$n_del_res)-median(slow$n_del_res))
  tab[j,19] <- median(fast$n_del_sens)
  tab[j,23] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
  fast <- subset(data, data$d29_response=="fast")
  slow <- subset(data, data$d29_response=="slow" | data$d29_response=="intermediate")
  tab[j,7] <- median(slow$n_del_res)
  tab[j,11] <- (median(fast$n_del_res)-median(slow$n_del_res))
  tab[j,20] <- median(slow$n_del_sens)
  tab[j,24] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
  fast <- subset(data, data$d79_response=="negative")
  slow <- subset(data, data$d79_response=="positive")
  tab[j,12] <- median(fast$n_del_res)
  tab[j,13] <- median(slow$n_del_res)
  tab[j,14] <- (median(fast$n_del_res)-median(slow$n_del_res))
  tab[j,25] <- median(fast$n_del_sens)
  tab[j,26] <- median(slow$n_del_sens)
  tab[j,27] <- (median(fast$n_del_sens)-median(slow$n_del_sens))}

# writing out the summary table
write.table(tab, "randomized_200_crispr_gene_medians.tab", sep="\t", quote=F, col.names=T, row.names=F)

## utilizing 100 random genes

# generating a table into which the results can be collected
tab <- matrix(ncol=27, nrow=1000)
colnames(tab) <- c("iteration", "res_del_median_fast15", "res_del_median_slow15", "res_del_median_diff15", "res_del_median_fast29", "res_del_median_fastmed29", "res_del_median_medslow29", "res_del_median_slow29", "res_del_median_diff29_f_s", "res_del_median_diff29_fm_s", "res_del_median_diff29_f_ms","res_del_median_fast79", "res_del_median_slow79", "res_del_median_diff79", "sens_del_median_fast15", "sens_del_median_slow15", "sens_del_median_diff15", "sens_del_median_fast29", "sens_del_median_fastmed29", "sens_del_median_medslow29", "sens_del_median_slow29", "sens_del_median_diff29_f_s", "sens_del_median_diff29_fm_s", "sens_del_median_diff29_f_ms","sens_del_median_fast79", "sens_del_median_slow79", "sens_del_median_diff79")
tab[,1] <- 1:1000
# conducting 1000 random iterations
for (j in 1:1000) {
  # sampling 100 random genes to be the resistance and sensitivity genes
  resgenes <- sample(vcr[,1], 100)
  sensgenes <- sample(vcr[,1], 100)
  genes <- append(resgenes, sensgenes)
  # generating a table into which the info about the CRISPR genes can be collected
  delmat <- matrix(ncol=10, nrow=1)
  colnames(delmat) <- c("chrom", "start", "end", "size", "cnv_type", "sample_id", "d15_response", "d29_response", "d79_response", "gene")
  for (g in genes) {
    # checking from the del table, which cases have deletions overlapping the gene region
    del2 <- subset(del, del[,1]==panel[g,1])
    loss <- vector()
    for (i in 1:nrow(del2)) {
      if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
      if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
      if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,2])<=as.numeric(panel[g,3]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
      if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}}
    del2 <- del2[loss,]
    gene <- matrix(ncol=1, nrow=nrow(del2))
    gene[,1] <- g
    loss <- cbind(del2, gene)
    colnames(loss) <- colnames(delmat)
    delmat <- rbind(delmat, loss)}
  delmat <- delmat[-1,]
  
  # generating a table summarizing the numbers of hits to resistance and sensitivity genes
  samples <- d29[,2]
  samples <- append(d29[,2], gepard[,1])
  d15_response <- vector()
  d29_response <- vector()
  d79_response <- vector()
  n_del_res <- vector()
  n_del_sens <- vector()
  for (s in samples) {
    if (s %in% d29[,2]) {
      resp <- subset(d15, d15[,2]==s)
      d15_response <- append(d15_response, resp[1,4])
      resp <- subset(d29, d29[,2]==s)
      d29_response <- append(d29_response, resp[1,4])
      resp <- subset(d79, d79[,2]==s)
      d79_response <- append(d79_response, resp[1,4])}
    if (s %in% gepard[,1]) {
      resp <- subset(gepard, gepard[,1]==s)
      d15_response <- append(d15_response, resp[1,2])
      d29_response <- append(d29_response, resp[1,3])
      d79_response <- append(d79_response, resp[1,4])}
    dels <- subset(delmat, delmat[,6]==s)
    hits <- subset(dels, dels[,10] %in% resgenes)
    n_del_res <- append(n_del_res, nrow(hits))
    hits <- subset(dels, dels[,10] %in% sensgenes)
    n_del_sens <- append(n_del_sens, nrow(hits))}
  data <- data.frame(samples, d15_response, d29_response, d79_response, n_del_res, n_del_sens)
  
  # checking the medians for each response group
  fast <- subset(data, data$d15_response=="fast")
  slow <- subset(data, data$d15_response=="slow")
  tab[j,2] <- median(fast$n_del_res)
  tab[j,3] <- median(slow$n_del_res)
  tab[j,4] <- (median(fast$n_del_res)-median(slow$n_del_res))
  tab[j,15] <- median(fast$n_del_sens)
  tab[j,16] <- median(slow$n_del_sens)
  tab[j,17] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
  fast <- subset(data, data$d29_response=="fast")
  slow <- subset(data, data$d29_response=="slow")
  tab[j,5] <- median(fast$n_del_res)
  tab[j,8] <- median(slow$n_del_res)
  tab[j,9] <- (median(fast$n_del_res)-median(slow$n_del_res))
  tab[j,18] <- median(fast$n_del_sens)
  tab[j,21] <- median(slow$n_del_sens)
  tab[j,22] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
  fast <- subset(data, data$d29_response=="fast" | data$d29_response=="intermediate")
  slow <- subset(data, data$d29_response=="slow")
  tab[j,6] <- median(fast$n_del_res)
  tab[j,10] <- (median(fast$n_del_res)-median(slow$n_del_res))
  tab[j,19] <- median(fast$n_del_sens)
  tab[j,23] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
  fast <- subset(data, data$d29_response=="fast")
  slow <- subset(data, data$d29_response=="slow" | data$d29_response=="intermediate")
  tab[j,7] <- median(slow$n_del_res)
  tab[j,11] <- (median(fast$n_del_res)-median(slow$n_del_res))
  tab[j,20] <- median(slow$n_del_sens)
  tab[j,24] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
  fast <- subset(data, data$d79_response=="negative")
  slow <- subset(data, data$d79_response=="positive")
  tab[j,12] <- median(fast$n_del_res)
  tab[j,13] <- median(slow$n_del_res)
  tab[j,14] <- (median(fast$n_del_res)-median(slow$n_del_res))
  tab[j,25] <- median(fast$n_del_sens)
  tab[j,26] <- median(slow$n_del_sens)
  tab[j,27] <- (median(fast$n_del_sens)-median(slow$n_del_sens))}

# writing out the summary table
write.table(tab, "randomized_100_crispr_gene_medians.tab", sep="\t", quote=F, col.names=T, row.names=F)

## utilizing a matching number of random genes to the FDR cutoff filtered drug screen genes

# generating tables into which the results can be collected
tab <- matrix(ncol=81, nrow=1000)
drugs <- c("VCR", "6MP", "ARAC", "DNR", "LASP", "MAF", "MTX", "DEX")
cols <- c("iteration")
for (d in drugs) {
  cols <- append(cols, paste0("res_del_median_diff15_", d))
  cols <- append(cols, paste0("sens_del_median_diff15_", d))
  cols <- append(cols, paste0("res_del_median_diff29_f_s_", d))
  cols <- append(cols, paste0("sens_del_median_diff29_f_s_", d))
  cols <- append(cols, paste0("res_del_median_diff29_fm_s_", d))
  cols <- append(cols, paste0("sens_del_median_diff29_fm_s_", d))
  cols <- append(cols, paste0("res_del_median_diff29_f_ms_", d))
  cols <- append(cols, paste0("sens_del_median_diff29_f_ms_", d))
  cols <- append(cols, paste0("res_del_median_diff79_", d))
  cols <- append(cols, paste0("sens_del_median_diff79_", d))}
colnames(tab) <- cols
tab[,1] <- 1:1000

# reading in the table with the numbers of genes per drug passing the FDR cutoff
ngenes <- read.table("crispr_hit_summary_d15.tab", header=T)
ngenes <- ngenes[,1:3]
rownames(ngenes) <- ngenes[,1]
ngenes <- ngenes[,-1]

# conducting 1000 random iterations
for (j in 1:1000) {
  # sampling 1435 random genes to be the resistance and sensitivity genes
  resgenes <- sample(vcr[,1], 1435)
  sensgenes <- sample(vcr[,1], 1435)
  genes <- append(resgenes, sensgenes)
  # generating a table into which the info about the CRISPR genes can be collected
  delmat <- matrix(ncol=10, nrow=1)
  colnames(delmat) <- c("chrom", "start", "end", "size", "cnv_type", "sample_id", "d15_response", "d29_response", "d79_response", "gene")
  for (g in genes) {
    # checking from the del table, which cases have deletions overlapping the gene region
    del2 <- subset(del, del[,1]==panel[g,1])
    loss <- vector()
    for (i in 1:nrow(del2)) {
      if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
      if (isTRUE(as.numeric(del2[i,2])<=as.numeric(panel[g,2]) & as.numeric(del2[i,3])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
      if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,2])<=as.numeric(panel[g,3]) & as.numeric(del2[i,3])>=as.numeric(panel[g,3]))) {loss <- append(loss, i)}
      if (isTRUE(as.numeric(del2[i,2])>=as.numeric(panel[g,2]) & as.numeric(del2[i,3])<=as.numeric(panel[g,3]))) {loss <- append(loss, i)}}
    del2 <- del2[loss,]
    gene <- matrix(ncol=1, nrow=nrow(del2))
    gene[,1] <- g
    loss <- cbind(del2, gene)
    colnames(loss) <- colnames(delmat)
    delmat <- rbind(delmat, loss)}
  delmat <- delmat[-1,]
  
  # conducting the following steps for each drug
  # i.e., choosing the correct amount of genes that matches the FDR filtered number of genes
  # and calculating the means for those
  for (d in drugs) {
    resgenes2 <- sample(resgenes, ngenes[d,1])
    sensgenes2 <- sample(sensgenes, ngenes[d,2])
    genes2 <- append(resgenes2, sensgenes2)
    delmat2 <- subset(delmat, delmat[,10] %in% genes2)
    
    # generating a table summarizing the numbers of hits to resistance and sensitivity genes
    samples <- d29[,2]
    samples <- append(d29[,2], gepard[,1])
    d15_response <- vector()
    d29_response <- vector()
    d79_response <- vector()
    n_del_res <- vector()
    n_del_sens <- vector()
    for (s in samples) {
      if (s %in% d29[,2]) {
        resp <- subset(d15, d15[,2]==s)
        d15_response <- append(d15_response, resp[1,4])
        resp <- subset(d29, d29[,2]==s)
        d29_response <- append(d29_response, resp[1,4])
        resp <- subset(d79, d79[,2]==s)
        d79_response <- append(d79_response, resp[1,4])}
      if (s %in% gepard[,1]) {
        resp <- subset(gepard, gepard[,1]==s)
        d15_response <- append(d15_response, resp[1,2])
        d29_response <- append(d29_response, resp[1,3])
        d79_response <- append(d79_response, resp[1,4])}
      dels <- subset(delmat2, delmat2[,6]==s)
      hits <- subset(dels, dels[,10] %in% resgenes2)
      n_del_res <- append(n_del_res, nrow(hits))
      hits <- subset(dels, dels[,10] %in% sensgenes2)
      n_del_sens <- append(n_del_sens, nrow(hits))}
    data <- data.frame(samples, d15_response, d29_response, d79_response, n_del_res, n_del_sens)
    
    # checking the medians for each response group
    # d15
    fast <- subset(data, data$d15_response=="fast")
    slow <- subset(data, data$d15_response=="slow")
    tab[j,paste0("res_del_median_diff15_", d)] <- (median(fast$n_del_res)-median(slow$n_del_res))
    tab[j,paste0("sens_del_median_diff15_", d)] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
    # d29 fast vs slow
    fast <- subset(data, data$d29_response=="fast")
    slow <- subset(data, data$d29_response=="slow")
    tab[j,paste0("res_del_median_diff29_f_s_", d)] <- (median(fast$n_del_res)-median(slow$n_del_res))
    tab[j,paste0("sens_del_median_diff29_f_s_", d)] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
    # d29 fastmed vs slow
    fast <- subset(data, data$d29_response=="fast" | data$d29_response=="intermediate")
    slow <- subset(data, data$d29_response=="slow")
    tab[j,paste0("res_del_median_diff29_fm_s_", d)] <- (median(fast$n_del_res)-median(slow$n_del_res))
    tab[j,paste0("sens_del_median_diff29_fm_s_", d)] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
    # d29 fast vs medslow
    fast <- subset(data, data$d29_response=="fast")
    slow <- subset(data, data$d29_response=="slow" | data$d29_response=="intermediate")
    tab[j,paste0("res_del_median_diff29_f_ms_", d)] <- (median(fast$n_del_res)-median(slow$n_del_res))
    tab[j,paste0("sens_del_median_diff29_f_ms_", d)] <- (median(fast$n_del_sens)-median(slow$n_del_sens))
    # d79
    fast <- subset(data, data$d79_response=="negative")
    slow <- subset(data, data$d79_response=="positive")
    tab[j,paste0("res_del_median_diff79_", d)] <- (median(fast$n_del_res)-median(slow$n_del_res))
    tab[j,paste0("sens_del_median_diff79_", d)] <- (median(fast$n_del_sens)-median(slow$n_del_sens))}}

# writing out the summary table
write.table(tab, "randomized_fdrcutoff_crispr_gene_medians.tab", sep="\t", quote=F, col.names=T, row.names=F)

# computing the FDR rates
# generating an empty file into which the results can be collected into
fdr <- matrix(ncol=30, nrow=8)
colnames(fdr) <- c("fdr_100_res_d15", "fdr_100_sens_d15", "fdr_200_res_d15", "fdr_200_sens_d15", "fdr_fdr_res_d15", "fdr_fdr_sens_d15", "fdr_100_res_d29fs", "fdr_100_sens_d29fs", "fdr_200_res_d29fs", "fdr_200_sens_d29fs", "fdr_fdr_res_d29fs", "fdr_fdr_sens_d29fs", "fdr_100_res_d29fm", "fdr_100_sens_d29fm", "fdr_200_res_d29fm", "fdr_200_sens_d29fm", "fdr_fdr_res_d29fm", "fdr_fdr_sens_d29fm", "fdr_100_res_d29ms", "fdr_100_sens_d29ms", "fdr_200_res_d29ms", "fdr_200_sens_d29ms", "fdr_fdr_res_d29ms", "fdr_fdr_sens_d29ms", "fdr_100_res_d79", "fdr_100_sens_d79", "fdr_200_res_d79", "fdr_200_sens_d79", "fdr_fdr_res_d79", "fdr_fdr_sens_d79")
fdr <- as.data.frame(fdr)

# 200 genes
random <- read.table("randomized_200_crispr_gene_medians.tab", header=T)
d15 <- read.table("crispr_top200_hit_summary_d15.tab", header=T)
d29 <- read.table("crispr_top200_hit_summary_d29.tab", header=T)
d79 <- read.table("crispr_top200_hit_summary_d79.tab", header=T)
drugs <- d15[,1]
for (d in 1:length(drugs)) {
  # d15 res
  diff <- (as.numeric(d15$res_del_median_fast[d])-as.numeric(d15$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff15) >= diff)
    fdr$fdr_200_res_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff15) <= diff)
    fdr$fdr_200_res_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff15) != 0)
    fdr$fdr_200_res_d15[d] <- (nrow(random2)/1000)}
  # d15 sens
  diff <- (as.numeric(d15$sens_del_median_fast[d])-as.numeric(d15$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff15) >= diff)
    fdr$fdr_200_sens_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff15) <= diff)
    fdr$fdr_200_sens_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff15) != 0)
    fdr$fdr_200_sens_d15[d] <- (nrow(random2)/1000)}
  # d29 fs res
  diff <- (as.numeric(d29$res_del_median_fast[d])-as.numeric(d29$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_s) >= diff)
    fdr$fdr_200_res_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_s) <= diff)
    fdr$fdr_200_res_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_s) != 0)
    fdr$fdr_200_res_d29fs[d] <- (nrow(random2)/1000)}
  # d29 fs sens
  diff <- (as.numeric(d29$sens_del_median_fast[d])-as.numeric(d29$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_s) >= diff)
    fdr$fdr_200_sens_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_s) <= diff)
    fdr$fdr_200_sens_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_s) != 0)
    fdr$fdr_200_sens_d29fs[d] <- (nrow(random2)/1000)}
  # d29 fm res
  diff <- (as.numeric(d29$res_del_median_fastmed[d])-as.numeric(d29$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_fm_s) >= diff)
    fdr$fdr_200_res_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_fm_s) <= diff)
    fdr$fdr_200_res_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_fm_s) != 0)
    fdr$fdr_200_res_d29fm[d] <- (nrow(random2)/1000)}
  # d29 fm sens
  diff <- (as.numeric(d29$sens_del_median_fastmed[d])-as.numeric(d29$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_fm_s) >= diff)
    fdr$fdr_200_sens_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_fm_s) <= diff)
    fdr$fdr_200_sens_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_fm_s) != 0)
    fdr$fdr_200_sens_d29fm[d] <- (nrow(random2)/1000)}
  # d29 ms res
  diff <- (as.numeric(d29$res_del_median_fast[d])-as.numeric(d29$res_del_median_medslow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_ms) >= diff)
    fdr$fdr_200_res_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_ms) <= diff)
    fdr$fdr_200_res_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_ms) != 0)
    fdr$fdr_200_res_d29ms[d] <- (nrow(random2)/1000)}
  # d29 ms sens
  diff <- (as.numeric(d29$sens_del_median_fast[d])-as.numeric(d29$sens_del_median_medslow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_ms) >= diff)
    fdr$fdr_200_sens_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_ms) <= diff)
    fdr$fdr_200_sens_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_ms) != 0)
    fdr$fdr_200_sens_d29ms[d] <- (nrow(random2)/1000)}
  # d79 res
  diff <- (as.numeric(d79$res_del_median_fast[d])-as.numeric(d79$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff79) >= diff)
    fdr$fdr_200_res_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff79) <= diff)
    fdr$fdr_200_res_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff79) != 0)
    fdr$fdr_200_res_d79[d] <- (nrow(random2)/1000)}
  # d79 sens
  diff <- (as.numeric(d79$sens_del_median_fast[d])-as.numeric(d79$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff79) >= diff)
    fdr$fdr_200_sens_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff79) <= diff)
    fdr$fdr_200_sens_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff79) != 0)
    fdr$fdr_200_sens_d79[d] <- (nrow(random2)/1000)}}

# 100 genes
random <- read.table("randomized_100_crispr_gene_medians.tab", header=T)
d15 <- read.table("crispr_top100_hit_summary_d15.tab", header=T)
d29 <- read.table("crispr_top100_hit_summary_d29.tab", header=T)
d79 <- read.table("crispr_top100_hit_summary_d79.tab", header=T)
drugs <- d15[,1]
for (d in 1:length(drugs)) {
  # d15 res
  diff <- (as.numeric(d15$res_del_median_fast[d])-as.numeric(d15$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff15) >= diff)
    fdr$fdr_100_res_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff15) <= diff)
    fdr$fdr_100_res_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff15) != 0)
    fdr$fdr_100_res_d15[d] <- (nrow(random2)/1000)}
  # d15 sens
  diff <- (as.numeric(d15$sens_del_median_fast[d])-as.numeric(d15$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff15) >= diff)
    fdr$fdr_100_sens_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff15) <= diff)
    fdr$fdr_100_sens_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff15) != 0)
    fdr$fdr_100_sens_d15[d] <- (nrow(random2)/1000)}
  # d29 fs res
  diff <- (as.numeric(d29$res_del_median_fast[d])-as.numeric(d29$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_s) >= diff)
    fdr$fdr_100_res_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_s) <= diff)
    fdr$fdr_100_res_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_s) != 0)
    fdr$fdr_100_res_d29fs[d] <- (nrow(random2)/1000)}
  # d29 fs sens
  diff <- (as.numeric(d29$sens_del_median_fast[d])-as.numeric(d29$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_s) >= diff)
    fdr$fdr_100_sens_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_s) <= diff)
    fdr$fdr_100_sens_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_s) != 0)
    fdr$fdr_100_sens_d29fs[d] <- (nrow(random2)/1000)}
  # d29 fm res
  diff <- (as.numeric(d29$res_del_median_fastmed[d])-as.numeric(d29$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_fm_s) >= diff)
    fdr$fdr_100_res_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_fm_s) <= diff)
    fdr$fdr_100_res_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_fm_s) != 0)
    fdr$fdr_100_res_d29fm[d] <- (nrow(random2)/1000)}
  # d29 fm sens
  diff <- (as.numeric(d29$sens_del_median_fastmed[d])-as.numeric(d29$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_fm_s) >= diff)
    fdr$fdr_100_sens_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_fm_s) <= diff)
    fdr$fdr_100_sens_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_fm_s) != 0)
    fdr$fdr_100_sens_d29fm[d] <- (nrow(random2)/1000)}
  # d29 ms res
  diff <- (as.numeric(d29$res_del_median_fast[d])-as.numeric(d29$res_del_median_medslow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_ms) >= diff)
    fdr$fdr_100_res_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_ms) <= diff)
    fdr$fdr_100_res_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff29_f_ms) != 0)
    fdr$fdr_100_res_d29ms[d] <- (nrow(random2)/1000)}
  # d29 ms sens
  diff <- (as.numeric(d29$sens_del_median_fast[d])-as.numeric(d29$sens_del_median_medslow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_ms) >= diff)
    fdr$fdr_100_sens_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_ms) <= diff)
    fdr$fdr_100_sens_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff29_f_ms) != 0)
    fdr$fdr_100_sens_d29ms[d] <- (nrow(random2)/1000)}
  # d79 res
  diff <- (as.numeric(d79$res_del_median_fast[d])-as.numeric(d79$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff79) >= diff)
    fdr$fdr_100_res_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff79) <= diff)
    fdr$fdr_100_res_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$res_del_median_diff79) != 0)
    fdr$fdr_100_res_d79[d] <- (nrow(random2)/1000)}
  # d79 sens
  diff <- (as.numeric(d79$sens_del_median_fast[d])-as.numeric(d79$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff79) >= diff)
    fdr$fdr_100_sens_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff79) <= diff)
    fdr$fdr_100_sens_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random$sens_del_median_diff79) != 0)
    fdr$fdr_100_sens_d79[d] <- (nrow(random2)/1000)}}

# FDR cutoff gene version of the analysis
random <- read.table("randomized_fdrcutoff_crispr_gene_medians.tab", header=T)
d15 <- read.table("crispr_hit_summary_d15.tab", header=T)
d29 <- read.table("crispr_hit_summary_d29.tab", header=T)
d79 <- read.table("crispr_hit_summary_d79.tab", header=T)
drugs <- d15[,1]
rownames(fdr) <- drugs
for (d in 1:length(drugs)) {
  # d15 res
  diff <- (as.numeric(d15$res_del_median_fast[d])-as.numeric(d15$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff15_", rownames(fdr)[d])]) >= diff)
    fdr$fdr_fdr_res_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff15_", rownames(fdr)[d])]) <= diff)
    fdr$fdr_fdr_res_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff15_", rownames(fdr)[d])]) != 0)
    fdr$fdr_fdr_res_d15[d] <- (nrow(random2)/1000)}
  # d15 sens
  diff <- (as.numeric(d15$sens_del_median_fast[d])-as.numeric(d15$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff15_", rownames(fdr)[d])]) >= diff)
    fdr$fdr_fdr_sens_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff15_", rownames(fdr)[d])]) <= diff)
    fdr$fdr_fdr_sens_d15[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff15_", rownames(fdr)[d])]) != 0)
    fdr$fdr_fdr_sens_d15[d] <- (nrow(random2)/1000)}
  # d29 fs res
  diff <- (as.numeric(d29$res_del_median_fast[d])-as.numeric(d29$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff29_f_s_", rownames(fdr)[d])]) >= diff)
    fdr$fdr_fdr_res_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff29_f_s_", rownames(fdr)[d])]) <= diff)
    fdr$fdr_fdr_res_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff29_f_s_", rownames(fdr)[d])]) != 0)
    fdr$fdr_fdr_res_d29fs[d] <- (nrow(random2)/1000)}
  # d29 fs sens
  diff <- (as.numeric(d29$sens_del_median_fast[d])-as.numeric(d29$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff29_f_s_", rownames(fdr)[d])]) >= diff)
    fdr$fdr_fdr_sens_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff29_f_s_", rownames(fdr)[d])]) <= diff)
    fdr$fdr_fdr_sens_d29fs[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff29_f_s_", rownames(fdr)[d])]) != 0)
    fdr$fdr_fdr_sens_d29fs[d] <- (nrow(random2)/1000)}
  # d29 fm res
  diff <- (as.numeric(d29$res_del_median_fastmed[d])-as.numeric(d29$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff29_fm_s_", rownames(fdr)[d])]) >= diff)
    fdr$fdr_fdr_res_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff29_fm_s_", rownames(fdr)[d])]) <= diff)
    fdr$fdr_fdr_res_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff29_fm_s_", rownames(fdr)[d])]) != 0)
    fdr$fdr_fdr_res_d29fm[d] <- (nrow(random2)/1000)}
  # d29 fm sens
  diff <- (as.numeric(d29$sens_del_median_fastmed[d])-as.numeric(d29$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff29_fm_s_", rownames(fdr)[d])]) >= diff)
    fdr$fdr_fdr_sens_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff29_fm_s_", rownames(fdr)[d])]) <= diff)
    fdr$fdr_fdr_sens_d29fm[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff29_fm_s_", rownames(fdr)[d])]) != 0)
    fdr$fdr_fdr_sens_d29fm[d] <- (nrow(random2)/1000)}
  # d29 ms res
  diff <- (as.numeric(d29$res_del_median_fast[d])-as.numeric(d29$res_del_median_medslow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff29_f_ms_", rownames(fdr)[d])]) >= diff)
    fdr$fdr_fdr_res_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff29_f_ms_", rownames(fdr)[d])]) <= diff)
    fdr$fdr_fdr_res_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff29_f_ms_", rownames(fdr)[d])]) != 0)
    fdr$fdr_fdr_res_d29ms[d] <- (nrow(random2)/1000)}
  # d29 ms sens
  diff <- (as.numeric(d29$sens_del_median_fast[d])-as.numeric(d29$sens_del_median_medslow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff29_f_ms_", rownames(fdr)[d])]) >= diff)
    fdr$fdr_fdr_sens_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff29_f_ms_", rownames(fdr)[d])]) <= diff)
    fdr$fdr_fdr_sens_d29ms[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff29_f_ms_", rownames(fdr)[d])]) != 0)
    fdr$fdr_fdr_sens_d29ms[d] <- (nrow(random2)/1000)}
  # d79 res
  diff <- (as.numeric(d79$res_del_median_fast[d])-as.numeric(d79$res_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff79_", rownames(fdr)[d])]) >= diff)
    fdr$fdr_fdr_res_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff79_", rownames(fdr)[d])]) <= diff)
    fdr$fdr_fdr_res_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("res_del_median_diff79_", rownames(fdr)[d])]) != 0)
    fdr$fdr_fdr_res_d79[d] <- (nrow(random2)/1000)}
  # d79 sens
  diff <- (as.numeric(d79$sens_del_median_fast[d])-as.numeric(d79$sens_del_median_slow[d]))
  if (isTRUE(diff > 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff79_", rownames(fdr)[d])]) >= diff)
    fdr$fdr_fdr_sens_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff < 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff79_", rownames(fdr)[d])]) <= diff)
    fdr$fdr_fdr_sens_d79[d] <- (nrow(random2)/1000)}
  if (isTRUE(diff == 0)) {
    random2 <- subset(random, as.numeric(random[,paste0("sens_del_median_diff79_", rownames(fdr)[d])]) != 0)
    fdr$fdr_fdr_sens_d79[d] <- (nrow(random2)/1000)}}

# writing the table out
fdr <- cbind(drugs, fdr)
write.table(fdr, "crispr_gene_median_diff_fdr_summary.tab", sep="\t", quote=F, col.names=T, row.names=F)
