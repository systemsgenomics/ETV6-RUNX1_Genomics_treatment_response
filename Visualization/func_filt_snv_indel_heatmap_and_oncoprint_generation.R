# describing the functionality filtered SNV/indel results
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(colorRamp2)
library(sjmisc)

# the significant results from panel + WGS analysis
tab <- read.table("sig_func_filt_panel_wgs_snv_fishers_test_results.tab", header=T)

# the table with the gene annotation information
geneannot <- read.table("gene_functions.txt", header=T, sep="\t")

# determining genes not arising from some of the d29 comparisons
notfm <- c("ETV6", "SETD1B", "IKZF2", "PTPRD", "DOT1L")
notms <- c("MSH2", "ARPP21", "ERG", "PTPRD")
notall <- c("ARPP21")

## day 29 fast + intermediate vs slow comparison

# a table for heatmap generation
genes <- setdiff(tab[,1], notfm)
plot <- matrix(ncol=6, nrow=length(genes))
rownames(plot) <- genes
colnames(plot) <- c("d15_fast", "d15_slow", "d29_fast_intermediate", "d29_slow", "d79_fast", "d79_slow")
# adding the fraction to the table
for (g in genes) {
  panel2 <- subset(tab, tab$genes==g)
  plot[g,1] <- panel2$d15_frac_fast[1]
  plot[g,2] <- panel2$d15_frac_slow[1]
  plot[g,3] <- panel2$d29fm_frac_fast[1]
  plot[g,4] <- panel2$d29fm_frac_slow[1]
  plot[g,5] <- panel2$d79_frac_fast[1]
  plot[g,6] <- panel2$d79_frac_slow[1]}

# table describing the p-value levels
plot2 <- matrix(ncol=6, nrow=length(genes))
rownames(plot2) <- genes
colnames(plot2) <- c("pval_d15", "pval_d152", "pval_d79", "pval_d792", "pval_d79", "pval_d792")
for (g in genes) {
  panel2 <- subset(tab, tab$genes==g)
  plot2[g,1] <- panel2$d15_p[1]
  plot2[g,2] <- panel2$d15_p[1]
  plot2[g,3] <- panel2$d29fm_p[1]
  plot2[g,4] <- panel2$d29fm_p[1]
  plot2[g,5] <- panel2$d79_p[1]
  plot2[g,6] <- panel2$d79_p[1]}
for (i in 1:nrow(plot2)) {
  for (j in 1:ncol(plot2)) {
    if (isTRUE(as.numeric(plot2[i,j])>=0.05)) {plot2[i,j] <- ""}
    if (isTRUE(as.numeric(plot2[i,j])<0.05 & as.numeric(plot2[i,j])>=0.01)) {plot2[i,j] <- "*"}
    if (isTRUE(as.numeric(plot2[i,j])<0.01 & as.numeric(plot2[i,j])>=0.001)) {plot2[i,j] <- "**"}
    if (isTRUE(as.numeric(plot2[i,j])<0.001)) {plot2[i,j] <- "***"}}}

# generating a matrix for the drug annotation
annot <- matrix(ncol=13, nrow=1)
for (g in genes) {
  geneannot2 <- subset(geneannot, geneannot[,1]==g)
  colnames(annot) <- colnames(geneannot)
  annot <- rbind(annot, geneannot2)}
annot <- annot[-1,]
rownames(annot) <- annot[,1]
annot <- annot[,-c(1,3)]

# defining colour functions
col_fun3 <- c(Fast="#a1d76a", Slow="#b2182b")
col_fun4 <- c(yes="goldenrod2", no="grey96")
col_fun5 <- c(resistance="salmon1", sensitivity="#5ab4ac", na.value="grey96")
col_fun6 <- c('B cell development'="plum3", 'transcriptional regulator'="powderblue", 'tumor suppressor'="chocolate", 'oncogene'="seagreen")
col_fun = colorRamp2(c(0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2), c("grey96", "mistyrose", "pink", "lightpink2", "lightcoral", "brown3", "#b2182b", "indianred4"))

# heatmap generation
pdf("func_filt_mut_heatmap_fastmed29.pdf", width=12, height=5)
Heatmap(plot, col=col_fun, column_order = colnames(plot), cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf(plot2[i, j]), x, y, gp = gpar(fontsize = 10))}, cluster_rows=FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd = 2), column_split = c(rep("Day 15", 2), rep("Day 29", 2), rep("Day 79", 2)), column_title_gp = gpar(col = c("black", "black", "black")), top_annotation=HeatmapAnnotation("Response"=rep(c("Fast", "Slow"), 3), col=list("Response"=col_fun3), gp = gpar(col="white", lwd=2)), show_column_names =F,  heatmap_legend_param = list(title="Fraction of cases with a mutation", legend_name_gp = gpar(fontsize = 5)), right_annotation=rowAnnotation(na_col="grey96", "Driver gene"=annot[,1], "Gene function"=annot[,2], "Function"=annot[,3], "VCR"=annot[,4], "6-MP"=annot[,5], "L-ASP"=annot[,6], "ARAC"=annot[,7], "DNR"=annot[,8], "MAF"=annot[,9], "MTX"=annot[,10], "DEX"=annot[,11], annotation_name_rot = 45, gp = gpar(col="white", lwd=2), col=list("Driver gene"=col_fun4, "Drug target gene"=col_fun4, "Gene function"=col_fun6, "Function"=col_fun6, "VCR"=col_fun5, "6-MP"=col_fun5, "L-ASP"=col_fun5, "ARAC"=col_fun5, "DNR"=col_fun5, "MAF"=col_fun5, "MTX"=col_fun5, "DEX"=col_fun5)))
dev.off()

## day 29 fast vs intermediate + slow comparison

# table for heatmap generation
genes <- setdiff(tab[,1], notms)
plot <- matrix(ncol=6, nrow=length(genes))
rownames(plot) <- genes
colnames(plot) <- c("d15_fast", "d15_slow", "d29_fast_intermediate", "d29_slow", "d79_fast", "d79_slow")
# adding the fraction to the table
for (g in genes) {
  panel2 <- subset(tab, tab$genes==g)
  plot[g,1] <- panel2$d15_frac_fast[1]
  plot[g,2] <- panel2$d15_frac_slow[1]
  plot[g,3] <- panel2$d29ms_frac_fast[1]
  plot[g,4] <- panel2$d29ms_frac_slow[1]
  plot[g,5] <- panel2$d79_frac_fast[1]
  plot[g,6] <- panel2$d79_frac_slow[1]}

# table for describing the p-value levels
plot2 <- matrix(ncol=6, nrow=length(genes))
rownames(plot2) <- genes
colnames(plot2) <- c("pval_d15", "pval_d152", "pval_d79", "pval_d792", "pval_d79", "pval_d792")
for (g in genes) {
  panel2 <- subset(tab, tab$genes==g)
  plot2[g,1] <- panel2$d15_p[1]
  plot2[g,2] <- panel2$d15_p[1]
  plot2[g,3] <- panel2$d29ms_p[1]
  plot2[g,4] <- panel2$d29ms_p[1]
  plot2[g,5] <- panel2$d79_p[1]
  plot2[g,6] <- panel2$d79_p[1]}
for (i in 1:nrow(plot2)) {
  for (j in 1:ncol(plot2)) {
    if (isTRUE(as.numeric(plot2[i,j])>=0.05)) {plot2[i,j] <- ""}
    if (isTRUE(as.numeric(plot2[i,j])<0.05 & as.numeric(plot2[i,j])>=0.01)) {plot2[i,j] <- "*"}
    if (isTRUE(as.numeric(plot2[i,j])<0.01 & as.numeric(plot2[i,j])>=0.001)) {plot2[i,j] <- "**"}
    if (isTRUE(as.numeric(plot2[i,j])<0.001)) {plot2[i,j] <- "***"}}}

# matrix for the drug annotation
annot <- matrix(ncol=13, nrow=1)
for (g in genes) {
  geneannot2 <- subset(geneannot, geneannot[,1]==g)
  colnames(annot) <- colnames(geneannot)
  annot <- rbind(annot, geneannot2)}
annot <- annot[-1,]
rownames(annot) <- annot[,1]
annot <- annot[,-c(1,3)]

# defining colour functions
col_fun3 <- c(Fast="#a1d76a", Slow="#b2182b")
col_fun4 <- c(yes="goldenrod2", no="grey96")
col_fun5 <- c(resistance="salmon1", sensitivity="#5ab4ac", na.value="grey96")
col_fun6 <- c('B cell development'="plum3", 'transcriptional regulator'="powderblue", 'tumor suppressor'="chocolate", 'oncogene'="seagreen")
col_fun = colorRamp2(c(0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2), c("grey96", "mistyrose", "pink", "lightpink2", "lightcoral", "brown3", "#b2182b", "indianred4"))

# heatmap generation
pdf("func_filt_mut_heatmap_medslow29.pdf", width=12, height=5)
Heatmap(plot, col=col_fun, column_order = colnames(plot), cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf(plot2[i, j]), x, y, gp = gpar(fontsize = 10))}, cluster_rows=FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd = 2), column_split = c(rep("Day 15", 2), rep("Day 29", 2), rep("Day 79", 2)), column_title_gp = gpar(col = c("black", "black", "black")), top_annotation=HeatmapAnnotation("Response"=rep(c("Fast", "Slow"), 3), col=list("Response"=col_fun3), gp = gpar(col="white", lwd=2)), show_column_names =F,  heatmap_legend_param = list(title="Fraction of cases with a mutation", legend_name_gp = gpar(fontsize = 5)),  right_annotation=rowAnnotation(na_col="grey96", "Driver gene"=annot[,1], "Gene function"=annot[,2], "Function"=annot[,3], "VCR"=annot[,4], "6-MP"=annot[,5], "L-ASP"=annot[,6], "ARAC"=annot[,7], "DNR"=annot[,8], "MAF"=annot[,9], "MTX"=annot[,10], "DEX"=annot[,11], annotation_name_rot = 45, gp = gpar(col="white", lwd=2), col=list("Driver gene"=col_fun4, "Drug target gene"=col_fun4, "Gene function"=col_fun6, "Function"=col_fun6, "VCR"=col_fun5, "6-MP"=col_fun5, "L-ASP"=col_fun5, "ARAC"=col_fun5, "DNR"=col_fun5, "MAF"=col_fun5, "MTX"=col_fun5, "DEX"=col_fun5)))
dev.off()

## day 29 fast vs intermediate vs slow comparison

# table for heatmap generation
genes <- setdiff(tab[,1], notall)
plot <- matrix(ncol=7, nrow=length(genes))
rownames(plot) <- genes
colnames(plot) <- c("d15_fast", "d15_slow", "d29_fast", "d29_intermediate", "d29_slow", "d79_fast", "d79_slow")
# adding the fraction to the table
for (g in genes) {
  panel2 <- subset(tab, tab$genes==g)
  plot[g,1] <- panel2$d15_frac_fast[1]
  plot[g,2] <- panel2$d15_frac_slow[1]
  plot[g,3] <- panel2$d29all_frac_fast[1]
  plot[g,4] <- panel2$d29all_frac_med[1]
  plot[g,5] <- panel2$d29all_frac_slow[1]
  plot[g,6] <- panel2$d79_frac_fast[1]
  plot[g,7] <- panel2$d79_frac_slow[1]}

# table describing the p-value levels
plot2 <- matrix(ncol=7, nrow=length(genes))
rownames(plot2) <- genes
colnames(plot2) <- c("pval_d15", "pval_d152", "pval_d79", "pval_d79", "pval_d792", "pval_d79", "pval_d792")
for (g in genes) {
  panel2 <- subset(tab, tab$genes==g)
  plot2[g,1] <- panel2$d15_p[1]
  plot2[g,2] <- panel2$d15_p[1]
  plot2[g,3] <- panel2$d29all_p[1]
  plot2[g,4] <- panel2$d29all_p[1]
  plot2[g,5] <- panel2$d29all_p[1]
  plot2[g,6] <- panel2$d79_p[1]
  plot2[g,7] <- panel2$d79_p[1]}
for (i in 1:nrow(plot2)) {
  for (j in 1:ncol(plot2)) {
    if (isTRUE(as.numeric(plot2[i,j])>=0.05)) {plot2[i,j] <- ""}
    if (isTRUE(as.numeric(plot2[i,j])<0.05 & as.numeric(plot2[i,j])>=0.01)) {plot2[i,j] <- "*"}
    if (isTRUE(as.numeric(plot2[i,j])<0.01 & as.numeric(plot2[i,j])>=0.001)) {plot2[i,j] <- "**"}
    if (isTRUE(as.numeric(plot2[i,j])<0.001)) {plot2[i,j] <- "***"}}}

# matrix for the drug annotation
annot <- matrix(ncol=13, nrow=1)
for (g in genes) {
  geneannot2 <- subset(geneannot, geneannot[,1]==g)
  colnames(annot) <- colnames(geneannot)
  annot <- rbind(annot, geneannot2)}
annot <- annot[-1,]
rownames(annot) <- annot[,1]
annot <- annot[,-c(1,3)]

# defining colour functions
col_fun3 <- c(Fast="#a1d76a", Intermediate="#f1a340", Slow="#b2182b")
col_fun4 <- c(yes="goldenrod2", no="grey96")
col_fun5 <- c(resistance="salmon1", sensitivity="#5ab4ac", na.value="grey96")
col_fun = colorRamp2(c(0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2), c("grey96", "mistyrose", "pink", "lightpink2", "lightcoral", "brown3", "#b2182b", "indianred4"))

# heatmap generation
pdf("func_filt_mut_heatmap_all29.pdf", width=12, height=5)
Heatmap(plot, col=col_fun, column_order = colnames(plot), cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf(plot2[i, j]), x, y, gp = gpar(fontsize = 10))}, cluster_rows=FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "white", lwd = 2), column_split = c(rep("Day 15", 2), rep("Day 29", 3), rep("Day 79", 2)), column_title_gp = gpar(col = c("black", "black", "black")), top_annotation=HeatmapAnnotation("Response"=c("Fast", "Slow", "Fast", "Intermediate", "Slow", "Fast", "Slow"), col=list("Response"=col_fun3), gp = gpar(col="white", lwd=2)), show_column_names =F,  heatmap_legend_param = list(title="Fraction of cases with a mutation", legend_name_gp = gpar(fontsize = 5)),  right_annotation=rowAnnotation(na_col="grey96", "Driver gene"=annot[,1], "Gene function"=annot[,2], "Function"=annot[,3], "VCR"=annot[,4], "6-MP"=annot[,5], "L-ASP"=annot[,6], "ARAC"=annot[,7], "DNR"=annot[,8], "MAF"=annot[,9], "MTX"=annot[,10], "DEX"=annot[,11], annotation_name_rot = 45, gp = gpar(col="white", lwd=2), col=list("Driver gene"=col_fun4, "Drug target gene"=col_fun4, "Gene function"=col_fun6, "Function"=col_fun6, "VCR"=col_fun5, "6-MP"=col_fun5, "L-ASP"=col_fun5, "ARAC"=col_fun5, "DNR"=col_fun5, "MAF"=col_fun5, "MTX"=col_fun5, "DEX"=col_fun5)))
dev.off()


## oncoprints visualizing functionally filtered top genes

# reading in the functionally filtered variant results
tab <- read.table("func_filtered_snv_indel_results.tab", sep="\t", header=T)
fisher <- read.table("sig_func_filt_panel_wgs_snv_fishers_test_results.tab", header=T)
genes <- as.vector(fisher[,1])

# saving the genes of interest as a vector
genes15 <- subset(fisher, as.numeric(fisher$d15_p) < 0.05)
genes15 <- as.vector(genes15[,1])
genes29 <- subset(fisher, as.numeric(fisher$d29all_p) < 0.05 | as.numeric(fisher$d29fm_p) < 0.05 | as.numeric(fisher$d29ms_p) < 0.05)
genes29 <- as.vector(genes29[,1])
genes79 <- subset(fisher, as.numeric(fisher$d79_p) < 0.05)
genes79 <- as.vector(genes79[,1])

# subsetting the table to only include the genes of interest
tab <- subset(tab, tab$gene %in% genes)
# correcting the variant type annotations
for (i in 1:nrow(tab)) {
  tab$mut_type[i] <- gsub("Missense_Mutation", "Missense", tab$mut_type[i])
  tab$mut_type[i] <- gsub("Nonsense_Mutation", "Nonsense", tab$mut_type[i])
  tab$mut_type[i] <- gsub("Frame_Shift_Del", "Frameshift", tab$mut_type[i])
  tab$mut_type[i] <- gsub("Frame_Shift_Ins", "Frameshift", tab$mut_type[i])
  tab$mut_type[i] <- gsub("Splice_Site", "Splicesite", tab$mut_type[i])}

# response info
mrd <- read.table("sample_responses.txt", header=T)
# substituting NA values with "Unknown", and modifying everything else to start with a capital letter
for (i in 1:nrow(mrd)) {
  for (j in 2:ncol(mrd)) {
    if (is.na(mrd[i,j]) | mrd[i,j]=="unknown") {mrd[i,j] <- c("Unknown")}
    if (mrd[i,j]=="slow") {mrd[i,j] <- c("Slow")}
    if (mrd[i,j]=="fast") {mrd[i,j] <- c("Fast")}
    if (mrd[i,j]=="intermediate") {mrd[i,j] <- c("Intermediate")}
    if (mrd[i,j]=="positive") {mrd[i,j] <- c("Positive")}
    if (mrd[i,j]=="negative") {mrd[i,j] <- c("Negative")}}}

# ordering the table according to day 29 response
slow <- subset(mrd, mrd$d29_response=="Slow")
med <- subset(mrd, mrd$d29_response=="Intermediate")
fast <- subset(mrd, mrd$d29_response=="Fast")
uk <- subset(mrd, mrd$d29_response=="Unknown")
mrd <- rbind(slow, med, fast)
# saving the samples in this order as well
samples <- mrd[,1]

# collecting the hits to a table
toplot <- matrix(ncol=length(genes29), nrow=length(samples))
colnames(toplot) <- genes29
rownames(toplot) <- samples
for (s in samples) {
  for (g in genes29) {
    tab2 <- subset(tab, tab[,1]==s & tab[,5]==g)
    if (nrow(tab2)!=0) {toplot[s,g] <- tab2[1,3]}}}

# defining the oncoprint function
alter_fun = list(background = function(x, y, w, h) {grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "grey96", col = NA))},
                 Nonsense = function(x, y, w, h) {grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill="tomato2", col = NA))},
                 Frameshift = function(x, y, w, h) {grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill="olivedrab", col = NA))},
                 Splicesite = function(x, y, w, h) {grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill="#5ab4ac", col = NA))},
                 Missense = function(x, y, w, h) {grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill="#cab2d6", col = NA))
                 })

# determining the colour functions
col_fun1 <- c(Fast="#a1d76a", Intermediate="#f1a340", Slow="#b2182b", Unknown="grey90")
col_fun2 <- c(Negative="white", Positive="black", Unknown="grey90")

# oncoprint generation
pdf("manuscript_SNV_indel_oncoprint_d29_all_genes.pdf", width=15, height=3)
oncoPrint(t(toplot), column_order=samples, alter_fun=alter_fun, show_pct=FALSE, show_column_names=FALSE, top_annotation=HeatmapAnnotation(na_col="grey90", "Mid-induction response"=mrd[,2], "EOI response"=mrd[,3], "EOC MRD status"=mrd[,4], col=list("Mid-induction response"=col_fun1, "EOI response"=col_fun1, "EOC MRD status"=col_fun2)))
dev.off()

## genes arising from EOC analyses

# reading in the response info
mrd <- read.table("sample_responses.txt", header=T)
# substituting NA values with "Unknown", and modifying everything else to start with a capital letter
for (i in 1:nrow(mrd)) {
  for (j in 2:ncol(mrd)) {
    if (is.na(mrd[i,j]) | mrd[i,j]=="unknown") {mrd[i,j] <- c("Unknown")}
    if (mrd[i,j]=="slow") {mrd[i,j] <- c("Slow")}
    if (mrd[i,j]=="fast") {mrd[i,j] <- c("Fast")}
    if (mrd[i,j]=="intermediate") {mrd[i,j] <- c("Intermediate")}
    if (mrd[i,j]=="positive") {mrd[i,j] <- c("Positive")}
    if (mrd[i,j]=="negative") {mrd[i,j] <- c("Negative")}}}
# ordering the table according to day 79 response
slow <- subset(mrd, mrd$d79_response=="Positive")
fast <- subset(mrd, mrd$d79_response=="Negative")
uk <- subset(mrd, mrd$d79_response=="Unknown")
mrd <- rbind(slow, fast)
# saving the samples in this order as well
samples <- mrd[,1]

# collecting the hits to a table
toplot <- matrix(ncol=length(genes79), nrow=length(samples))
colnames(toplot) <- genes79
rownames(toplot) <- samples
for (s in samples) {
  for (g in genes79) {
    tab2 <- subset(tab, tab[,1]==s & tab[,5]==g)
    if (nrow(tab2)!=0) {toplot[s,g] <- tab2[1,3]}}}

# oncoprint generation
pdf("manuscript_SNV_indel_oncoprint_d79_all_genes.pdf", width=15, height=3)
oncoPrint(t(toplot), column_order=samples, alter_fun=alter_fun, show_pct=FALSE, show_column_names=FALSE, top_annotation=HeatmapAnnotation(na_col="grey90", "Mid-induction response"=mrd[,2], "EOI response"=mrd[,3], "EOC MRD status"=mrd[,4], col=list("Mid-induction response"=col_fun1, "EOI response"=col_fun1, "EOC MRD status"=col_fun2)))
dev.off()

## genes arising from mid-induction analyses

# reading in the response info
mrd <- read.table("sample_responses.txt", header=T)
# substituting NA values with "Unknown", and modifying everything else to start with a capital letter
for (i in 1:nrow(mrd)) {
  for (j in 2:ncol(mrd)) {
    if (is.na(mrd[i,j]) | mrd[i,j]=="unknown") {mrd[i,j] <- c("Unknown")}
    if (mrd[i,j]=="slow") {mrd[i,j] <- c("Slow")}
    if (mrd[i,j]=="fast") {mrd[i,j] <- c("Fast")}
    if (mrd[i,j]=="intermediate") {mrd[i,j] <- c("Intermediate")}
    if (mrd[i,j]=="positive") {mrd[i,j] <- c("Positive")}
    if (mrd[i,j]=="negative") {mrd[i,j] <- c("Negative")}}}

# ordering the table according to day 15 response
slow <- subset(mrd, mrd$d15_response=="Slow")
fast <- subset(mrd, mrd$d15_response=="Fast")
uk <- subset(mrd, mrd$d15_response=="Unknown")
mrd <- rbind(slow, fast)
# saving the samples in this order as well
samples <- mrd[,1]

# collecting the hits to a table
toplot <- matrix(ncol=length(genes15), nrow=length(samples))
colnames(toplot) <- genes15
rownames(toplot) <- samples
for (s in samples) {
  for (g in genes15) {
    tab2 <- subset(tab, tab[,1]==s & tab[,5]==g)
    if (nrow(tab2)!=0) {toplot[s,g] <- tab2[1,3]}}}

# oncoprint generation
pdf("manuscript_SNV_indel_oncoprint_d15_all_genes.pdf", width=15, height=2)
oncoPrint(t(toplot), column_order=samples, alter_fun=alter_fun, show_pct=FALSE, show_column_names=FALSE, top_annotation=HeatmapAnnotation(na_col="grey90", "Mid-induction response"=mrd[,2], "EOI response"=mrd[,3], "EOC MRD status"=mrd[,4], col=list("Mid-induction response"=col_fun1, "EOI response"=col_fun1, "EOC MRD status"=col_fun2)))
dev.off()
