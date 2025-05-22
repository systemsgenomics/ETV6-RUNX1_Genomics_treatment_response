# generate oncoprints of the IG rearrangements with functionality predictions
library(sjmisc)
library(ggplot2)
library(gridExtra)
library(circlize)
library(ComplexHeatmap)

# reading in the table with the comprehensive functionality predictions
tab <- read.table("collective_IGH_IGK_IGL_alterations_functionalities.txt", sep="\t", header=T)

# simplifying the productive/unproductive colunms
for (i in 1:nrow(tab)) {
  if (isTRUE(str_contains(tab$Functionality, "Productive"))) {tab$Functionality <- c("productive")}
  if (isTRUE(str_contains(tab$Functionality, "Unproductive"))) {tab$Functionality <- c("unproductive")}
  if (isFALSE(str_contains(tab$Functionality, "productive")) & isFALSE(str_contains(tab$Functionality, "unproductive"))) {tab$Functionality <- c("unknown")}
  if (isTRUE(str_contains(tab$IMGT_functionality, "Productive"))) {tab$IMGT_functionality <- c("productive")}
  if (isTRUE(str_contains(tab$IMGT_functionality, "Unproductive"))) {tab$IMGT_functionality <- c("unproductive")}
  if (isFALSE(str_contains(tab$IMGT_functionality, "productive")) & isFALSE(str_contains(tab$IMGT_functionality, "unproductive"))) {tab$IMGT_functionality <- c("unknown")}
  if (isTRUE(str_contains(tab$IgBlast_functionality, "Productive"))) {tab$IgBlast_functionality <- c("productive")}
  if (isTRUE(str_contains(tab$IgBlast_functionality, "Unproductive"))) {tab$IgBlast_functionality <- c("unproductive")}
  if (isFALSE(str_contains(tab$IgBlast_functionality, "productive")) & isFALSE(str_contains(tab$IgBlast_functionality, "unproductive"))) {tab$IgBlast_functionality <- c("unknown")}
}
# simplifying the table
tab <- cbind(tab$PatientID, tab$Analysis, tab$Annotation, tab$Mechanism, tab$Functionality, tab$IMGT_functionality, tab$IgBlast_functionality)

# then determining the response groups for the patients
# reading in the MRD data
mrd <- read.table("MRD_status_new.txt", header=T)
# and ordering it according to the flow day 29 values
mrd <- mrd[order(-mrd$day29_flow),]
# saving the samples in a vector in the MRD order
samples <- as.vector(mrd[,1])

# separating the IGH, IGK and IGL tables, and also making one with both light chains
igh <- subset(tab, tab[,2]=="IGH")
igk <- subset(tab, tab[,2]=="IGK")
igl <- subset(tab, tab[,2]=="IGL")
light <- subset(tab, tab[,2]=="IGK" | tab[,2]=="IGL")

# defining the alter_fun for the oncoprint generation
alter_fun = list(background = function(x, y, w, h) {grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))},
                 igcaller_productive = function(x, y, w, h) {grid.polygon(unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h), gp = gpar(fill="tomato1", col ="white"))},
                 igblast_productive = function(x, y, w, h) {grid.polygon(unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h), gp = gpar(fill="tomato3", col ="white"))},
                 igcaller_unproductive = function(x, y, w, h) {grid.polygon(unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h), gp = gpar(fill="mediumturquoise", col ="white"))},
                 igblast_unproductive = function(x, y, w, h) {grid.polygon(unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h), gp = gpar(fill="#5ab4ac", col ="white"))},
                 igcaller_unknown = function(x, y, w, h) {grid.polygon(unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h), gp = gpar(fill="grey58", col ="white"))},
                 igblast_unknown = function(x, y, w, h) {grid.polygon(unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h), gp = gpar(fill="grey58", col ="white"))})
# heatmap annotation colour functions
col_fun1 = colorRamp2(c(0,0.0001,0.01,0.1,0.8), c("#f7f7f7", "#d1e5f0", "#91bfdb", "#af8dc3", "#b2182b"))
col_fun2 = colorRamp2(c(0,0.0001,0.0002,0.001,0.03), c("#f7f7f7", "#d1e5f0", "#91bfdb", "#af8dc3", "#b2182b"))
col_fun3 <- c(fast="#a1d76a", intermediate="#f1a340", slow="#b2182b")
col_fun4 <- c(positive="black", unknown="grey", negative="#f7f7f7")

# day 29 MRD order and IGH rearrangements
# saving the gene combinations as a vector
genes <- unique(igh[,3])
# and generating a matrix which can be used to generate the plot
toplot <- matrix(ncol=length(samples), nrow=length(genes))
colnames(toplot) <- samples
rownames(toplot) <- genes
# and collecting the hits into the table, including the functionalities
for (g in genes) {
  igh2 <- subset(igh, igh[,3]==g)
  for (i in 1:nrow(igh2)) {
    toplot[g,igh2[i,1]] <- paste0("igcaller_", igh2[i,5], ",", "igblast_", igh2[i,7])}}
# drawing the oncoprint
pdf("IGH_rearrangements_with_functionalities_oncoprint_d29.pdf", width=11)
oncoPrint(toplot, column_order=samples, row_order=genes, alter_fun=alter_fun, show_pct=FALSE, show_column_names=TRUE, left_annotation=NULL, top_annotation=HeatmapAnnotation("day15_flow"=mrd[,2], "day29_flow"=mrd[,3], "day79_status"=mrd[,5], "response"=mrd[,4], col=list("day15_flow"=col_fun1, "day29_flow"=col_fun2, "response"=col_fun3, "day79_status"=col_fun4)))
dev.off()

# day 29 MRD order and IGK rearrangements
# saving the gene combinations as a vector
genes <- unique(igk[,3])
# and generating a matrix which can be used to generate the plot
toplot <- matrix(ncol=length(samples), nrow=length(genes))
colnames(toplot) <- samples
rownames(toplot) <- genes
# and collecting the hits into the table, including the functionalities
for (g in genes) {
  igk2 <- subset(igk, igk[,3]==g)
  for (i in 1:nrow(igk2)) {
    toplot[g,igk2[i,1]] <- paste0("igcaller_", igk2[i,5], ",", "igblast_", igk2[i,7])}}
# drawing the oncoprint
pdf("IGK_rearrangements_with_functionalities_oncoprint_d29.pdf", width=11)
oncoPrint(toplot, column_order=samples, row_order=genes, alter_fun=alter_fun, show_pct=FALSE, show_column_names=TRUE, left_annotation=NULL, top_annotation=HeatmapAnnotation("day15_flow"=mrd[,2], "day29_flow"=mrd[,3], "day79_status"=mrd[,5], "response"=mrd[,4], col=list("day15_flow"=col_fun1, "day29_flow"=col_fun2, "response"=col_fun3, "day79_status"=col_fun4)))
dev.off()

# day 29 MRD order and IGL rearrangements
# saving the gene combinations as a vector
genes <- unique(igl[,3])
# and generating a matrix which can be used to generate the plot
toplot <- matrix(ncol=length(samples), nrow=length(genes))
colnames(toplot) <- samples
rownames(toplot) <- genes
# and collecting the hits into the table, including the functionalities
for (g in genes) {
  igl2 <- subset(igl, igl[,3]==g)
  for (i in 1:nrow(igl2)) {
    toplot[g,igl2[i,1]] <- paste0("igcaller_", igl2[i,5], ",", "igblast_", igl2[i,7])}}
# drawing the oncoprint
pdf("IGL_rearrangements_with_functionalities_oncoprint_d29.pdf", width=11)
oncoPrint(toplot, column_order=samples, row_order=genes, alter_fun=alter_fun, show_pct=FALSE, show_column_names=TRUE, left_annotation=NULL, top_annotation=HeatmapAnnotation("day15_flow"=mrd[,2], "day29_flow"=mrd[,3], "day79_status"=mrd[,5], "response"=mrd[,4], col=list("day15_flow"=col_fun1, "day29_flow"=col_fun2, "response"=col_fun3, "day79_status"=col_fun4)))
dev.off()


# day 15 MRD order

# ordering the mrd table according to the flow day 15 values
mrd <- mrd[order(-mrd$day15_flow),]
# saving the samples in a vector in the MRD order
samples <- as.vector(mrd[,1])

# day 15 MRD order and IGH rearrangements
# saving the gene combinations as a vector
genes <- unique(igh[,3])
# and generating a matrix which can be used to generate the plot
toplot <- matrix(ncol=length(samples), nrow=length(genes))
colnames(toplot) <- samples
rownames(toplot) <- genes
# and collecting the hits into the table, including the functionalities
for (g in genes) {
  igh2 <- subset(igh, igh[,3]==g)
  for (i in 1:nrow(igh2)) {
    toplot[g,igh2[i,1]] <- paste0("igcaller_", igh2[i,5], ",", "igblast_", igh2[i,7])}}
# drawing the oncoprint
pdf("IGH_rearrangements_with_functionalities_oncoprint_d15.pdf", width=11)
oncoPrint(toplot, column_order=samples, row_order=genes, alter_fun=alter_fun, show_pct=FALSE, show_column_names=TRUE, left_annotation=NULL, top_annotation=HeatmapAnnotation("day15_flow"=mrd[,2], "day29_flow"=mrd[,3], "day79_status"=mrd[,5], "response"=mrd[,4], col=list("day15_flow"=col_fun1, "day29_flow"=col_fun2, "response"=col_fun3, "day79_status"=col_fun4)))
dev.off()

# day 15 MRD order and IGK rearrangements
# saving the gene combinations as a vector
genes <- unique(igk[,3])
# and generating a matrix which can be used to generate the plot
toplot <- matrix(ncol=length(samples), nrow=length(genes))
colnames(toplot) <- samples
rownames(toplot) <- genes
# and collecting the hits into the table, including the functionalities
for (g in genes) {
  igk2 <- subset(igk, igk[,3]==g)
  for (i in 1:nrow(igk2)) {
    toplot[g,igk2[i,1]] <- paste0("igcaller_", igk2[i,5], ",", "igblast_", igk2[i,7])}}
# drawing the oncoprint
pdf("IGK_rearrangements_with_functionalities_oncoprint_d15.pdf", width=11)
oncoPrint(toplot, column_order=samples, row_order=genes, alter_fun=alter_fun, show_pct=FALSE, show_column_names=TRUE, left_annotation=NULL, top_annotation=HeatmapAnnotation("day15_flow"=mrd[,2], "day29_flow"=mrd[,3], "day79_status"=mrd[,5], "response"=mrd[,4], col=list("day15_flow"=col_fun1, "day29_flow"=col_fun2, "response"=col_fun3, "day79_status"=col_fun4)))
dev.off()

# day 15 MRD order and IGL rearrangements
# saving the gene combinations as a vector
genes <- unique(igl[,3])
# and generating a matrix which can be used to generate the plot
toplot <- matrix(ncol=length(samples), nrow=length(genes))
colnames(toplot) <- samples
rownames(toplot) <- genes
# and collecting the hits into the table, including the functionalities
for (g in genes) {
  igl2 <- subset(igl, igl[,3]==g)
  for (i in 1:nrow(igl2)) {
    toplot[g,igl2[i,1]] <- paste0("igcaller_", igl2[i,5], ",", "igblast_", igl2[i,7])}}
# drawing the oncoprint
pdf("IGL_rearrangements_with_functionalities_oncoprint_d15.pdf", width=11)
oncoPrint(toplot, column_order=samples, row_order=genes, alter_fun=alter_fun, show_pct=FALSE, show_column_names=TRUE, left_annotation=NULL, top_annotation=HeatmapAnnotation("day15_flow"=mrd[,2], "day29_flow"=mrd[,3], "day79_status"=mrd[,5], "response"=mrd[,4], col=list("day15_flow"=col_fun1, "day29_flow"=col_fun2, "response"=col_fun3, "day79_status"=col_fun4)))
dev.off()


# day 79 MRD order

# ordering the mrd table according to the day 79 MRD values
samples <- append(positives, unknown)
samples <- append(samples, negatives)
# ordering the mrd table according to day 79 values
mrd2 <- matrix(ncol=5, nrow=33)
mrd2[,1] <- samples
for (i in 1:nrow(mrd2)) {
  mrd3 <- subset(mrd, mrd[,1]==mrd2[i,1])
  mrd2[i,2] <- mrd3[1,2]
  mrd2[i,3] <- mrd3[1,3]
  mrd2[i,4] <- mrd3[1,4]
  mrd2[i,5] <- mrd3[1,5]}
mrd <- mrd2
mrd[,2] <- as.numeric(mrd[,2])
mrd[,3] <- as.numeric(mrd[,3])

# day 79 MRD order and IGH rearrangements
# saving the gene combinations as a vector
genes <- unique(igh[,3])
# and generating a matrix which can be used to generate the plot
toplot <- matrix(ncol=length(samples), nrow=length(genes))
colnames(toplot) <- samples
rownames(toplot) <- genes
# and collecting the hits into the table, including the functionalities
for (g in genes) {
  igh2 <- subset(igh, igh[,3]==g)
  for (i in 1:nrow(igh2)) {
    toplot[g,igh2[i,1]] <- paste0("igcaller_", igh2[i,5], ",", "igblast_", igh2[i,7])}}
# drawing the oncoprint
pdf("IGH_rearrangements_with_functionalities_oncoprint_d79.pdf", width=11)
oncoPrint(toplot, column_order=samples, row_order=genes, alter_fun=alter_fun, show_pct=FALSE, show_column_names=TRUE, left_annotation=NULL, top_annotation=HeatmapAnnotation("day15_flow"=as.numeric(mrd[,2]), "day29_flow"=as.numeric(mrd[,3]), "day79_status"=mrd[,5], "response"=mrd[,4], col=list("day15_flow"=col_fun1, "day29_flow"=col_fun2, "response"=col_fun3, "day79_status"=col_fun4)))
dev.off()

# day 79 MRD order and IGK rearrangements
# saving the gene combinations as a vector
genes <- unique(igk[,3])
# and generating a matrix which can be used to generate the plot
toplot <- matrix(ncol=length(samples), nrow=length(genes))
colnames(toplot) <- samples
rownames(toplot) <- genes
# and collecting the hits into the table, including the functionalities
for (g in genes) {
  igk2 <- subset(igk, igk[,3]==g)
  for (i in 1:nrow(igk2)) {
    toplot[g,igk2[i,1]] <- paste0("igcaller_", igk2[i,5], ",", "igblast_", igk2[i,7])}}
# drawing the oncoprint
pdf("IGK_rearrangements_with_functionalities_oncoprint_d79.pdf", width=11)
oncoPrint(toplot, column_order=samples, row_order=genes, alter_fun=alter_fun, show_pct=FALSE, show_column_names=TRUE, left_annotation=NULL, top_annotation=HeatmapAnnotation("day15_flow"=as.numeric(mrd[,2]), "day29_flow"=as.numeric(mrd[,3]), "day79_status"=mrd[,5], "response"=mrd[,4], col=list("day15_flow"=col_fun1, "day29_flow"=col_fun2, "response"=col_fun3, "day79_status"=col_fun4)))
dev.off()

# day 79 MRD order and IGL rearrangements
# saving the gene combinations as a vector
genes <- unique(igl[,3])
# and generating a matrix which can be used to generate the plot
toplot <- matrix(ncol=length(samples), nrow=length(genes))
colnames(toplot) <- samples
rownames(toplot) <- genes
# and collecting the hits into the table, including the functionalities
for (g in genes) {
  igl2 <- subset(igl, igl[,3]==g)
  for (i in 1:nrow(igl2)) {
    toplot[g,igl2[i,1]] <- paste0("igcaller_", igl2[i,5], ",", "igblast_", igl2[i,7])}}
# drawing the oncoprint
pdf("IGL_rearrangements_with_functionalities_oncoprint_d79.pdf", width=11)
oncoPrint(toplot, column_order=samples, row_order=genes, alter_fun=alter_fun, show_pct=FALSE, show_column_names=TRUE, left_annotation=NULL, top_annotation=HeatmapAnnotation("day15_flow"=as.numeric(mrd[,2]), "day29_flow"=as.numeric(mrd[,3]), "day79_status"=mrd[,5], "response"=mrd[,4], col=list("day15_flow"=col_fun1, "day29_flow"=col_fun2, "response"=col_fun3, "day79_status"=col_fun4)))
dev.off()
