# Calculate the gene set and regulon scores

#read ETV6::RUNX1 related regulons (ER+ and normal cell comparison, Mehtonen el al. 2020 Genome Medicine)
er <- read.table("Fig5a_heatmap_regulons.txt")
er

#read table with all target genes (Mehtonen el al. 2020 Genome Medicine)
HCA_aucell_targets <- read.csv(gzfile("regulon_targets.csv.gz"))

#select only uniq genes
reg <- er[["V1"]]

for (i in reg) {
  targets <- HCA_aucell_targets[HCA_aucell_targets$regulon == i, ]
  targets_igt2 <- targets[targets$target %in% names(table(targets$target))[table(targets$target) > 2], ]
  targets_uniq <-unique(targets_igt2$target)
  write.table(targets_uniq, paste0(i, "_uniq_targets.txt"), row.names=F, col.names = F, quote=F)
}

#create vector including regulons with genes left
reg.new <- read.table("selected_regulons.txt", header = F)
reg.new <- reg.new[["V1"]]

#create the list including regulons and their genes
ll <- list()

for (r in 1:length(reg.new)) {
  genes <- read.table(paste0(reg.new[r], "_uniq_targets.txt"), header=F)
  genes <- genes[["V1"]]
  ll[[r]]=genes
}

names(ll) <- reg.new

#GSVA analysis
#ssGSVA
library(GSVA)

#read count matrix
cpm = read.table("cpm_matrix_combo.tab", sep="\t", header = T)

#match ensembl id to gene name
library(org.Hs.eg.db)
cpm$symbol = mapIds(org.Hs.eg.db, keys = rownames(cpm), keytype = "ENSEMBL", column = "SYMBOL")

library(dplyr)

#remove rows without symbol
cpm_symbol <- cpm %>% 
  filter(!is.na(symbol))

#remove rows with duplicate symbol
cpm_symbol <- cpm_symbol[!duplicated(cpm_symbol[,52]),]

#move symbols in columns to row names
rownames(cpm_symbol) <- cpm_symbol[,52]
cpm_symbol <- cpm_symbol[-52]

#change into matrix format
cpm_symbol <- as.matrix(cpm_symbol)

#give zeros a small value
cpm_symbol[cpm_symbol==0] = 0.001
cpm_symbol

#log2 change
cpm_symbol_log2 <- log2(cpm_symbol)

#summary of log2 cpms of all genes per case
sum <- summary(cpm_symbol_log2)

#calculate the gene set scores
gsva.es.reg <- gsva(cpm_symbol_log2,ll, method = "ssgsea", min.sz > 1, verbose = FALSE) 
dim(gsva.es.reg)
head(gsva.es.reg, 10)
gsva.df.reg <- data.frame(gsva.es.reg)
summary(gsva.df.reg)

#Correlation calculations, MRD vs. score
#IDs as EOI order
col_order_d29 <- mrd_flow_comb$ID
col_order_d29

#reorganize score matrix to EOI (day29) order
score.day29 <- as.matrix(gsva.df.reg)
score.day29 <- score.day29[, col_order_d29]
score.day29

cor29 = matrix(nrow = 47, ncol = 2)
cor29

#calculate Spearman's rank correlation coefficient
for (i in 1:nrow(score.day29)) {
  cor_29=cor.test(mrd_flow_comb[,4], score.day29[i,], method = "spearman")
  cor29[i,1]=cor_29$p.value
  cor29[i,2]=cor_29$estimate }
colnames(cor29)= c("p-value", "cor")
rownames(cor29)= row.names(score.day29)
write.table(cor29, "correlation.txt", sep="\t")



#Heatmap visualization

library(grid)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

#top annotation color annotations
col_fun1 = colorRamp2(c(0,0.0001,0.01,0.1,0.9), c("#F7F7F7", "#D1E5F0", "#91BFDB", "#AF8DC3", "#B2182B"))
col_fun2 = colorRamp2(c(0,0.0001,0.0002,0.001,0.03), c("#F7F7F7", "#D1E5F0", "#91BFDB", "#AF8DC3", "#B2182B"))
col_fun3 = c("fast"="#A1D76A", "intermediate" = "#F1A340", "slow"="#B2182B")
col_fun5 = colorRamp2(c(0, 0.01, 0.1 ,0.5, 0.7, 0.8, 0.9), c("#F7F7F7", "#D1E5F0", "#91BFDB", "#AF8DC3", "#c38d96", "#e69598","#B2182B"))
col_fun4 <- c(positive="black", unknown="grey", negative="#F7F7F7")

#top annotation with mrd and response data
ha_d29=HeatmapAnnotation("dg_blast"=mrd_d29[,2],"day15_flow"=mrd_d29[,3], "day29_flow"=mrd_d29[,4],"day79_status"=mrd_d29[,6], "response_d29"=mrd_d29[,5], 
                         col=list("dg_blast"=col_fun5, "day15_flow"=col_fun1, "day29_flow"=col_fun2, "response_d29"=col_fun3, "day79_status"=col_fun4, 
                                  row_names_gp = gpar(fontsize = 3)))

#row scaling
score.day29 =  t(scale(t(score.day29)))
summary(score.day29)

color_score= colorRamp2(c(-4, 0, 3 ), c("steelblue3", "#f1eef6", "#980043"))

reg_29_comb <- Heatmap(score.day29, name = "regulon_combo_score_d29", col = color_score,
                       show_column_names = T, column_names_gp = gpar(fontsize = 6), 
                       show_row_names = TRUE, row_names_gp = gpar(fontsize = 6),
                       cluster_rows = T,
                       cluster_columns = F,
                       top_annotation = ha_d29,
                       column_order = col_order_d29)


## Gene set scored for cell cycle associated E2F regulons ##

library(tidyverse)

#select the E2F regulons only
e2f <- HCA_aucell_targets %>%
  filter(str_detect(regulon, "E2F"))

e2f.regulons <-unique(e2f$regulon)
e2f.regulons

#select uniq genes only
for (e in e2f.regulons) {
  targets <- HCA_aucell_targets[HCA_aucell_targets$regulon == e, ]
  targets_igt2 <- targets[targets$target %in% names(table(targets$target))[table(targets$target) > 2], ]
  targets_uniq <-unique(targets_igt2$target)
  write.table(targets_uniq, paste0(e, "_uniq_targets.txt"), row.names=F, col.names = F, quote=F)
}

#create vector including regulon with genes left
e2f.reg <- read.table("~/patient_RNAseq/ER_manuska/regulon/e2f.reg.txt", header = F)
e2f.reg <- e2f.reg[["V1"]]

#create the list including regulons and their genes
ll.e2f <- list()

for (f in 1:length(e2f.reg)) {
  genes <- read.table(paste0(e2f.reg[f], "_uniq_targets.txt"), header=F)
  genes <- genes[["V1"]]
  ll.e2f[[f]]=genes
}

names(ll.e2f) <- e2f.reg

#calculate the gene set scores 
gsva.es.e2f <- gsva(cpm_symbol_log2,ll.e2f, method = "ssgsea", min.sz > 1, verbose = FALSE) 
dim(gsva.es.e2f)
head(gsva.es.e2f, 10)
gsva.df.e2f <- data.frame(gsva.es.e2f)
summary(gsva.df.e2f)


#Correlation calculations, MRD vs. score
#reorganize score matrix to EOI (day29) order
score.2ef <- as.matrix(gsva.df.e2f)
score.2ef_29 <- score.2ef[, col_order_d29]
score.2ef_29

cor29 = matrix(nrow = 10, ncol = 2)
cor29

#calculate Spearman's rank correlation coefficient
for (i in 1:nrow(score.2ef_29)) {
  cor_29=cor.test(mrd_flow_comb[,4], score.2ef_29[i,], method = "spearman")
  cor29[i,1]=cor_29$p.value
  cor29[i,2]=cor_29$estimate }
colnames(cor29)= c("p-value", "cor")
rownames(cor29)= row.names(score.2ef_29)
write.table(cor29, "correlation.txt", sep="\t")


#Heatmap visualization

#row scalling
score.2ef_29 =  t(scale(t(score.2ef_29)))
summary(score.2ef_29)

color_score1= colorRamp2(c(-3, 0, 3 ), c("steelblue3", "#f1eef6", "#980043"))

reg_29_ef <- Heatmap(score.2ef_29, name = "regulon_score_combo_d29", col = color_score1,
                     show_column_names = T, column_names_gp = gpar(fontsize = 6), 
                     show_row_names = TRUE, row_names_gp = gpar(fontsize = 6),
                     cluster_rows = T,
                     cluster_columns = F,
                     top_annotation = ha_d29,
                     column_order = col_order_d29)


## CELL CYCLE GENES, from Tirosh et al, 2015 to calculate cell cycle scores##

s.genes <- c("MCM5","PCNA","TYMS", "FEN1","MCM2" ,"MCM4","RRM1","UNG", "GINS2","MCM6" ,"CDCA7" ,"DTL",
             "PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP","RAD51AP1", "GMNN" ,"WDR76" ,"SLBP","CCNE2",
             "UBR7","POLD3","MSH2", "ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM",
             "CASP8AP2","USP1","CLSPN","POLA1" , "CHAF1B","BRIP1","E2F8")

g2m.genes <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5",
               "TPX2","TOP2A","NDC80", "CKS2", "NUF2","CKS1B","MKI67","TMPO","CENPF" ,"TACC3","FAM64A" ,"SMC4",
               "CCNB2", "CKAP2L","CKAP2", "AURKB" ,"BUB1" ,"KIF11","ANP32E","TUBB4B" , "GTSE1" ,"KIF20B","HJURP","CDCA3",
               "HN1" ,"CDC20","TTK" , "CDC25C","KIF2C","RANGAP1", "NCAPD2" ,"DLGAP5","CDCA2" , "CDCA8" ,"ECT2" ,"KIF23",
               "HMMR", "AURKA","PSRC1","ANLN","LBR", "CKAP5" , "CENPE" , "CTCF" ,"NEK2" ,"G2E3" , "GAS2L3","CBX5","CENPA")

#create a list containing a vector, a matrix and a list
list_data <- list(s.genes, g2m.genes)

#give names to the elements in the list
names(list_data) <- c("S Phase", "G2M Phase")

#calculate the gene set scores  
gsva.es.new <- gsva(cpm_symbol_log2,list_data, method = "ssgsea", min.sz > 1, verbose = FALSE) 
dim(gsva.es.new)
head(gsva.es.new, 10)
gsva.df.new <- data.frame(gsva.es.new)
summary(gsva.df.new)

score <- as.matrix(gsva.df.new)

#Correlation calculations, MRD vs. score
#reorganize score matrix to EOI (day29) order
score_29 <- score[, col_order_d29]
score_29

cor29 = matrix(nrow = 2, ncol = 2)
cor29

#calculate Spearman's rank correlation coefficient
for (i in 1:nrow(score_29)) {
  cor_29=cor.test(mrd_flow_comb[,4], score_29[i,], method = "spearman")
  cor29[i,1]=cor_29$p.value
  cor29[i,2]=cor_29$estimate }
colnames(cor29)= c("p-value", "cor")
rownames(cor29)= row.names(score_29)
write.table(cor29, "correlation.txt", sep="\t")

#Heatmap visualization

#row scalling 
score_29 =  t(scale(t(score_29)))
summary(score_29)

color_score= colorRamp2(c(-2, 0, 2 ), c("steelblue3", "#f1eef6", "#980043"))

sc_29 <- Heatmap(score_29, name = "cell_cycle_score_combo_d29", col = color_score,
                 show_column_names = T, column_names_gp = gpar(fontsize = 8), 
                 show_row_names = TRUE, row_names_gp = gpar(fontsize = 6),
                 cluster_rows = F,
                 cluster_columns = F,
                 column_order = col_order_d29)

#combine E2F and cell cycle heatmaps into same image
ht_list = reg_29_ef %v% sc_29
draw(ht_list)
