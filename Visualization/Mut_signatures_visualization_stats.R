#Mutation signature visualization and statistics 

# Load necessary library
library(grid)
library(ComplexHeatmap)
library(circlize)
library(dplyr)


#Create top annotation with MRD and treatment response data
#Read EOI (day29) mrd and response data
mrd=read.table("MRD_status.txt", header=T)
samples=as.vector(mrd[,1])

#color annotations
col_fun1 = colorRamp2(c(0,0.0001,0.01,0.1,0.8), c("#F7F7F7", "#D1E5F0", "#91BFDB", "#AF8DC3", "#B2182B"))
col_fun2 = colorRamp2(c(0,0.0001,0.0002,0.001,0.03), c("#F7F7F7", "#D1E5F0", "#91BFDB", "#AF8DC3", "#B2182B"))
col_fun3 = c("fast"="#A1D76A", "intermediate" = "#F1A340", "slow"="#B2182B")
col_fun4 <- c(positive="black", unknown="grey", negative="#F7F7F7")

#EOI (day29) top annotation
ha_29=HeatmapAnnotation("day15_flow"=mrd[,2], "day29_flow"=mrd[,3],"day79_status"=mrd[,5], "response"=mrd[,4], 
                     col=list("day15_flow"=col_fun1, "day29_flow"=col_fun2, "response"=col_fun3, "day79_status"=col_fun4))

#Read mid-induction (day15) mrd and response data
mrd_15=read.table("MRD_status_day15.txt", header=T)
samples_d15=as.vector(mrd_15[,1]) #tehdään vectori näytteistä

#color annotations
col_fun1 = colorRamp2(c(0,0.0001,0.01,0.1,0.8), c("#F7F7F7", "#D1E5F0", "#91BFDB", "#AF8DC3", "#B2182B"))
col_fun2 = colorRamp2(c(0,0.0001,0.0002,0.001,0.03), c("#F7F7F7", "#D1E5F0", "#91BFDB", "#AF8DC3", "#B2182B"))
col_fun3 = c("fast"="#A1D76A", "intermediate" = "#F1A340", "slow"="#B2182B")
col_fun5 = c("fast"="#A1D76A", "slow"="#B2182B")
col_fun4 <- c(positive="black", unknown="grey", negative="#F7F7F7")

#Mid-induction (day15) top annotation
ha_d15=HeatmapAnnotation("day15_flow"=mrd_15[,2], "day29_flow"=mrd_15[,3],"day79_status"=mrd_15[,6], "response"=mrd_15[,4], "response_d15"=mrd_15[,5],  
                            col=list("day15_flow"=col_fun1, "day29_flow"=col_fun2, "response"=col_fun3, "day79_status"=col_fun4, "response_d15"= col_fun5))


## SigprofilerAssignment COSMIC v3.4 ##

#Read signature contributions
data_v3=read.table("Signature_contributions_COSMIC_v3.4_mod.txt",header=T,na.strings="null",row.names=1)
data_v3=data.matrix(data_v3)
data_v3

#Reorganize data according to the sample order (EOI)
data3 = matrix(ncol=33, nrow=86)
colnames(data3) = samples
data3
signatures = as.vector(rownames(data_v3))
rownames(data3) = signatures
data3

for (i in signatures) {
  for (j in samples) {
    data3[i,j] = as.numeric(data_v3[i,j])}}

data3=as.matrix(data3)
data3

#color annotation for signature contributions
color = colorRamp2(c(0, 0.5, 1), c("snow", "brown3", "brown4"))

#Heatmap
Heatmap(data3, 
        name = "signature contribution", 
        top_annotation = ha_29, 
        column_order=samples, 
        column_names_gp = gpar(fontsize=9), 
        row_names_gp = gpar(fontsize=9), col = color)

#select only specific signatures
data3_sel <- data3[c(2:3, 19),]
data3_sel

Heatmap(data3_sel, 
        name = "signature contribution", 
        top_annotation = ha_29, 
        column_order=samples, 
        cluster_rows = F,
        column_names_gp = gpar(fontsize=9), 
        row_names_gp = gpar(fontsize=9), col = color)

#Reorganize data according to the sample order (mid-induction, day15)
data_d15_v3 = matrix(ncol=33, nrow=86)
colnames(data_d15_v3) = samples_d15
signatures = as.vector(rownames(data_v3))
rownames(data_d15_v3) = signatures
data_d15_v3

for (i in signatures) {
  for (j in samples_d15) {
    data_d15_v3[i,j] = as.numeric(data_v3[i,j])}}

data_d15_v3=as.matrix(data_d15_v3)
data_d15_v3

#color annotation for signature contributions
color = colorRamp2(c(0, 0.5, 1), c("snow", "brown3", "brown4"))

#Heatmap
Heatmap(data_d15_v3, 
        name = "signature contribution", 
        top_annotation = ha_d15, 
        column_order=samples_d15, 
        column_names_gp = gpar(fontsize=9), 
        row_names_gp = gpar(fontsize=5), col = color)


#Correlations
#EOI (day29)

cor29 = matrix(nrow = 86, ncol = 2)
cor29
#calculate Spearman's rank correlation coefficient
for (i in 1:nrow(data3_mod)) {
  cor_29=cor.test(mrd_mod[,3], data3_mod[i,], method = "spearman")
  cor29[i,1]=cor_29$p.value
  cor29[i,2]=cor_29$estimate }
colnames(cor29)= c("p-value", "cor")
write.table(cor29, "spearman_cor29_COSMIC_v3.4.txt", sep="\t")

#Mid-induction (day15)
cor15 = matrix(nrow = 86, ncol = 2)
cor15
#calculate Spearman's rank correlation coefficient
for (i in 1:nrow(data_d15_v3_mod)) {
  cor_15=cor.test(mrd_15_mod[,2], data_d15_v3_mod[i,], method = "spearman")
  cor15[i,1]=cor_15$p.value
  cor15[i,2]=cor_15$estimate }
colnames(cor15)= c("p-value", "cor")
write.table(cor15, "spearman_cor15_COSMIC_v3.4.txt", sep="\t")

### FOR WGS DATA ###
## COSMIC v3.4(SigProfilerAssignment) and v2(Musica) signatures ##

#Read signature contributions
data_comb_15=read.table("combined_day15_v2_v3.4.txt",header=T,na.strings="null",row.names=1)
data_comb_15=data.matrix(data_comb_15)

#Reorganize data according to the sample order (mid-induction, day15)
data_d15_v3 = matrix(ncol=32, nrow=8)
colnames(data_d15_v3) = samples_d15
data_d15_v3

signatures = as.vector(rownames(data_comb_15))
rownames(data_d15_v3) = signatures
data_d15_v3

for (i in signatures) {
  for (j in samples_d15) {
    data_d15_v3[i,j] = as.numeric(data_comb_15[i,j])}}

data_d15_v3=as.matrix(data_d15_v3)
data_d15_v3

#Color annotation for signature contributions
color = colorRamp2(c(0, 0.5, 1), c("snow", "brown3", "brown4"))

#Row annotation for v2 and v3.4 phases day15
#Define group labels
row_groups <- ifelse(grepl("^SBS", rownames(data_d15_v3)), "v3.4", "v2")

#Define colors
group_colors <- c("v3.4" = "#b6cc85", "v2" = "#769828")

#Create the annotation
row_anno <- rowAnnotation(Version = factor(row_groups, levels = c("v3.4", "v2")),
                          col = list(Version = group_colors))
#Heatmap
Heatmap(data_d15_v3, 
        name = "Signature contribution", 
        top_annotation = ha_d15,
        left_annotation = row_anno,
        column_order=samples_d15, 
        cluster_rows = F,
        column_names_gp = gpar(fontsize=9), 
        row_names_gp = gpar(fontsize=9), col = color)