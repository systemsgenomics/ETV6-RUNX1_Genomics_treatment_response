# filter TNscope variants based on VEP consequences
library(dplyr)

# reading in the variant table
tab <- read.table("collective_VEP_ann_moderate.tab")

# determining the column names
colnames(tab) <- c("patientID", "chrom", "pos", "ref", "alt", "tlod", "af_tumor", "af_normal", "afdp_tumor", "afdp_normal", "variant_type", "impact", "consequence", "gene")

# selecting only variants that have HIGH/MODERATE annotated impact
tab <- subset(tab, tab$impact=="HIGH" | tab$impact=="MODERATE")

# and then including only the most severe consequence if multiple are reported
for (i in 1:nrow(tab)){
  vector <- unlist(strsplit(tab$consequence[i], "&"))
  tab$consequence[i] <- vector[1]}
tab <- unique(tab)

# writing the table out
write.table(tab, "collective_high_moderate_impact_moderately_filtered_tnscope_variants.tab", col.names=F, row.names=F, quote=F, sep="\t")
