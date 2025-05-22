# Fisher's exact tests on combined WGS and panelseq SNV/indel results

# panelseq SNV/indel results
panel <- read.table("collective.SNV_indel.af0.05.vardict.csq_filt.tab", header=T)

# WGS SNV/indel results
snv <- read.table("collective_high_moderate_impact_filtered_tnscope_variants.tab", header=T)

## combined cohort analyses

# response info
resp <- read.table("sample_responses.txt", header=T)

# combining the data, and checking which genes are affected in both cohorts
genes <- intersect(snv[,15], panel[,20])

# simplifying both tables and combining the tables
snv <- subset(snv, snv$gene %in% genes)
panel <- subset(panel, panel$gene %in% genes)
panel <- panel[,c(1,6,18,20,21)]
snv <- snv[,c(3,12,14,15,1)]
colnames(snv) <- colnames(panel)
tab <- rbind(snv, panel)

# adding columns to the table describing the response groups and cohort
newcols <- matrix(ncol=4, nrow=nrow(tab))
for (i in 1:nrow(tab)) {
  resp2 <- subset(resp, resp[,1]==tab[i,5])
  newcols[i,1] <- resp2[1,2]
  newcols[i,2] <- resp2[1,3]
  newcols[i,3] <- resp2[1,4]
  if (isTRUE(str_contains(tab[i,5], "GE"))) {
    newcols[i,4] <- "WGS"}
  if (isFALSE(str_contains(tab[i,5], "GE"))) {
    newcols[i,4] <- "panel"}}
colnames(newcols) <- c("d15_response", "d29_response", "d79_response", "cohort")
tab <- cbind(tab, newcols)

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

# table into which the results can be collected
newcols <- matrix(ncol=24, nrow=length(genes))
colnames(newcols) <- c("d15_p", "d29all_p", "d29fm_p", "d29ms_p", "d79_p", "d15_top_prevalence", "d29all_top_prevalence", "d29fm_top_prevalence", "d29ms_top_prevalence", "d79_top_prevalence", "d15_frac_fast", "d15_frac_slow", "d29all_frac_fast", "d29all_frac_med", "d29all_frac_slow", "d29fm_frac_fast", "d29fm_frac_slow", "d29ms_frac_fast", "d29ms_frac_slow", "d79_frac_fast", "d79_frac_slow", "n_WGS", "n_panel", "sd_WGS_panel")
genes <- cbind(genes, newcols)
rownames(genes) <- genes[,1]
genes <- as.data.frame(genes)
# Fisher's exact tests for the genes affected by SNVs and indels
for (i in 1:nrow(genes)) {
  tab2 <- subset(tab, tab[,4]==genes[i,1])
  # d15
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,6]=="slow")
  fast <- subset(tab2, tab2[,6]=="fast")
  test[1,1] <- length(unique(slow[,5]))
  test[2,1] <- length(unique(fast[,5]))
  test[1,2] <- (nrow(slow15)-length(unique(slow[,5])))
  test[2,2] <- (nrow(fast15)-length(unique(fast[,5])))
  fisher <- fisher.test(test)
  genes$d15_p[i] <- round(fisher$p, digits=4)
  genes$d15_frac_slow[i] <- round(length(unique(slow[,5]))/nrow(slow15), digits=3)
  genes$d15_frac_fast[i] <- round(length(unique(fast[,5]))/nrow(fast15), digits=3)
  # d29 all
  test <- matrix(ncol=2, nrow=3)
  slow <- subset(tab2, tab2[,7]=="slow")
  med <- subset(tab2, tab2[,7]=="intermediate")
  fast <- subset(tab2, tab2[,7]=="fast")
  test[1,1] <- length(unique(slow[,5]))
  test[2,1] <- length(unique(med[,5]))
  test[3,1] <- length(unique(fast[,5]))
  test[1,2] <- (nrow(slow29)-length(unique(slow[,5])))
  test[2,2] <- (nrow(med29)-length(unique(med[,5])))
  test[3,2] <- (nrow(fast29)-length(unique(fast[,5])))
  fisher <- fisher.test(test)
  genes$d29all_p[i] <- round(fisher$p, digits=4)
  genes$d29all_frac_slow[i] <- round(length(unique(slow[,5]))/nrow(slow29), digits=3)
  genes$d29all_frac_med[i] <- round(length(unique(med[,5]))/nrow(med29), digits=3)
  genes$d29all_frac_fast[i] <- round(length(unique(fast[,5]))/nrow(fast29), digits=3)
  # d29 fast + intermediate vs slow
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,7]=="slow")
  fast <- subset(tab2, tab2[,7]=="fast" | tab2[,7]=="intermediate")
  test[1,1] <- length(unique(slow[,5]))
  test[2,1] <- length(unique(fast[,5]))
  test[1,2] <- (nrow(slow29)-length(unique(slow[,5])))
  test[2,2] <- (nrow(fastmed29)-length(unique(fast[,5])))
  fisher <- fisher.test(test)
  genes$d29fm_p[i] <- round(fisher$p, digits=4)
  genes$d29fm_frac_slow[i] <- round(length(unique(slow[,5]))/nrow(slow29), digits=3)
  genes$d29fm_frac_fast[i] <- round(length(unique(fast[,5]))/nrow(fastmed29), digits=3)
  # d29 fast vs intermediate + slow
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,7]=="slow" | tab2[,7]=="intermediate")
  fast <- subset(tab2, tab2[,7]=="fast")
  test[1,1] <- length(unique(slow[,5]))
  test[2,1] <- length(unique(fast[,5]))
  test[1,2] <- (nrow(medslow29)-length(unique(slow[,5])))
  test[2,2] <- (nrow(fast29)-length(unique(fast[,5])))
  fisher <- fisher.test(test)
  genes$d29ms_p[i] <- round(fisher$p, digits=4)
  genes$d29ms_frac_slow[i] <- round(length(unique(slow[,5]))/nrow(medslow29), digits=3)
  genes$d29ms_frac_fast[i] <- round(length(unique(fast[,5]))/nrow(fast29), digits=3)
  # d79
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,8]=="positive")
  fast <- subset(tab2, tab2[,8]=="negative")
  test[1,1] <- length(unique(slow[,5]))
  test[2,1] <- length(unique(fast[,5]))
  test[1,2] <- (nrow(slow79)-length(unique(slow[,5])))
  test[2,2] <- (nrow(fast79)-length(unique(fast[,5])))
  fisher <- fisher.test(test)
  genes$d79_p[i] <- round(fisher$p, digits=4)
  genes$d79_frac_slow[i] <- round(length(unique(slow[,5]))/nrow(slow79), digits=3)
  genes$d79_frac_fast[i] <- round(length(unique(fast[,5]))/nrow(fast79), digits=3)}

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

# adding info about the mutation numbers in each dataset, and the standard deviation
for (i in 1:nrow(genes)) {
  tab2 <- subset(tab, tab[,4]==genes[i,1])
  wgs <- subset(tab2, tab2$cohort=="WGS")
  panel <- subset(tab2, tab2$cohort=="panel")
  genes$n_WGS[i] <- length(unique(wgs[,5]))
  genes$n_panel[i] <- length(unique(panel[,5]))
  genes$sd_WGS_panel[i] <- sd(c((length(unique(wgs[,5]))/33), (length(unique(panel[,5]))/261)))}

# removing hits with the standard deviation > 0.15 between the WGS and panelseq
remove <- vector()
for (i in 1:nrow(genes)) {
  if (as.numeric(genes[i,25])>0.15) (remove <- append(remove, i))}
genes <- genes[-remove,]

# writing the table out
write.table(genes, "panel_wgs_snv_fishers_test_results.tab", col.names = T, row.names = F, sep="\t", quote=F)

# distinguishing the more significant genes using the p-value cutoff of 0.1
genes <- as.matrix(genes)
genes <- subset(genes, as.numeric(genes[,2]) < 0.1 | as.numeric(genes[,3]) < 0.1 | as.numeric(genes[,4]) < 0.1 | as.numeric(genes[,5]) < 0.1 | as.numeric(genes[,6]) < 0.1)

# removing the BCR and TCR genes from the comparison
remove <- vector()
for (i in 1:nrow(genes)) {
  if (isTRUE(str_contains(genes[i,1], "IG")) | isTRUE(str_contains(genes[i,1], "TRB"))) {remove <- append(remove, i)}}
genes <- genes[-remove,]

# writing the table out
write.table(genes, "sig_panel_wgs_snv_fishers_test_results.tab", col.names = T, row.names = F, sep="\t", quote=F)


## analyses only using panelseq data

# panelseq SNV/indel results
panel <- read.table("collective.SNV_indel.af0.05.vardict.csq_filt.tab", header=T)

# response info
resp <- read.table("sample_responses_panel.txt", header=T)

# checking which genes are affected in the data
genes <- unique(panel[,20])

# simplifying the table
tab <- panel[,c(1,6,18,20,21)]
tab <- unique(tab)

# adding columns to the table describing the response groups
newcols <- matrix(ncol=3, nrow=nrow(tab))
for (i in 1:nrow(tab)) {
  resp2 <- subset(resp, resp[,1]==tab[i,5])
  newcols[i,1] <- resp2[1,2]
  newcols[i,2] <- resp2[1,3]
  newcols[i,3] <- resp2[1,4]}
colnames(newcols) <- c("d15_response", "d29_response", "d79_response")
tab <- cbind(tab, newcols)

# generating tables of the samples in each responder group
fast15 <- subset(resp, resp$d15_response=="fast")
slow15 <- subset(resp, resp$d15_response=="slow")
fast29 <- subset(resp, resp$d29_response=="fast")
med29 <- subset(resp, resp$d29_response=="intermediate")
slow29 <- subset(resp, resp$d29_response=="slow")
medslow29 <- subset(resp, resp$d29_response=="slow" | resp$d29_response=="intermediate")
fastmed29 <- subset(resp, resp$d29_response=="fast" | resp$d29_response=="intermediate")
fast79 <- subset(resp, resp[,4]=="negative")
slow79 <- subset(resp, resp[,4]=="positive")

# generating a table into which the results can be collected
newcols <- matrix(ncol=22, nrow=length(genes))
colnames(newcols) <- c("d15_p", "d29all_p", "d29fm_p", "d29ms_p", "d79_p", "d15_top_prevalence", "d29all_top_prevalence", "d29fm_top_prevalence", "d29ms_top_prevalence", "d79_top_prevalence", "d15_frac_fast", "d15_frac_slow", "d29all_frac_fast", "d29all_frac_med", "d29all_frac_slow", "d29fm_frac_fast", "d29fm_frac_slow", "d29ms_frac_fast", "d29ms_frac_slow", "d79_frac_fast", "d79_frac_slow", "n_mut")
genes <- cbind(genes, newcols)
rownames(genes) <- genes[,1]
genes <- as.data.frame(genes)
for (i in 1:nrow(genes)) {
  tab2 <- subset(tab, tab[,4]==genes[i,1])
  # d15
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,6]=="slow")
  fast <- subset(tab2, tab2[,6]=="fast")
  test[1,1] <- length(unique(slow[,5]))
  test[2,1] <- length(unique(fast[,5]))
  test[1,2] <- (nrow(slow15)-length(unique(slow[,5])))
  test[2,2] <- (nrow(fast15)-length(unique(fast[,5])))
  fisher <- fisher.test(test)
  genes$d15_p[i] <- round(fisher$p, digits=4)
  genes$d15_frac_slow[i] <- round(length(unique(slow[,5]))/nrow(slow15), digits=3)
  genes$d15_frac_fast[i] <- round(length(unique(fast[,5]))/nrow(fast15), digits=3)
  # d29 all
  test <- matrix(ncol=2, nrow=3)
  slow <- subset(tab2, tab2[,7]=="slow")
  med <- subset(tab2, tab2[,7]=="intermediate")
  fast <- subset(tab2, tab2[,7]=="fast")
  test[1,1] <- length(unique(slow[,5]))
  test[2,1] <- length(unique(med[,5]))
  test[3,1] <- length(unique(fast[,5]))
  test[1,2] <- (nrow(slow29)-length(unique(slow[,5])))
  test[2,2] <- (nrow(med29)-length(unique(med[,5])))
  test[3,2] <- (nrow(fast29)-length(unique(fast[,5])))
  fisher <- fisher.test(test)
  genes$d29all_p[i] <- round(fisher$p, digits=4)
  genes$d29all_frac_slow[i] <- round(length(unique(slow[,5]))/nrow(slow29), digits=3)
  genes$d29all_frac_med[i] <- round(length(unique(med[,5]))/nrow(med29), digits=3)
  genes$d29all_frac_fast[i] <- round(length(unique(fast[,5]))/nrow(fast29), digits=3)
  # d29 fast + intermediate vs fast
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,7]=="slow")
  fast <- subset(tab2, tab2[,7]=="fast" | tab2[,7]=="intermediate")
  test[1,1] <- length(unique(slow[,5]))
  test[2,1] <- length(unique(fast[,5]))
  test[1,2] <- (nrow(slow29)-length(unique(slow[,5])))
  test[2,2] <- (nrow(fastmed29)-length(unique(fast[,5])))
  fisher <- fisher.test(test)
  genes$d29fm_p[i] <- round(fisher$p, digits=4)
  genes$d29fm_frac_slow[i] <- round(length(unique(slow[,5]))/nrow(slow29), digits=3)
  genes$d29fm_frac_fast[i] <- round(length(unique(fast[,5]))/nrow(fastmed29), digits=3)
  # d29 fast vs intermediate + fast
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,7]=="slow" | tab2[,7]=="intermediate")
  fast <- subset(tab2, tab2[,7]=="fast")
  test[1,1] <- length(unique(slow[,5]))
  test[2,1] <- length(unique(fast[,5]))
  test[1,2] <- (nrow(medslow29)-length(unique(slow[,5])))
  test[2,2] <- (nrow(fast29)-length(unique(fast[,5])))
  fisher <- fisher.test(test)
  genes$d29ms_p[i] <- round(fisher$p, digits=4)
  genes$d29ms_frac_slow[i] <- round(length(unique(slow[,5]))/nrow(medslow29), digits=3)
  genes$d29ms_frac_fast[i] <- round(length(unique(fast[,5]))/nrow(fast29), digits=3)
  # d79
  test <- matrix(ncol=2, nrow=2)
  slow <- subset(tab2, tab2[,8]=="positive")
  fast <- subset(tab2, tab2[,8]=="negative")
  test[1,1] <- length(unique(slow[,5]))
  test[2,1] <- length(unique(fast[,5]))
  test[1,2] <- (nrow(slow79)-length(unique(slow[,5])))
  test[2,2] <- (nrow(fast79)-length(unique(fast[,5])))
  fisher <- fisher.test(test)
  genes$d79_p[i] <- round(fisher$p, digits=4)
  genes$d79_frac_slow[i] <- round(length(unique(slow[,5]))/nrow(slow79), digits=3)
  genes$d79_frac_fast[i] <- round(length(unique(fast[,5]))/nrow(fast79), digits=3)}

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

# adding info about the total mutation numbers
for (i in 1:nrow(genes)) {
  tab2 <- subset(tab, tab[,4]==genes[i,1])
  genes$n_mut[i] <- length(unique(tab2[,5]))}

# removing genes with less than 3 hits altogether
genes <- as.matrix(genes)
genes <- subset(genes, as.numeric(genes[,23]) > 2)

# writing the table out
write.table(genes, "panel_only_snv_fishers_test_results.tab", col.names = T, row.names = F, sep="\t", quote=F)

# distinguishing the more significant genes using the p-value cutoff of 0.1
genes <- subset(genes, as.numeric(genes[,2]) < 0.1 | as.numeric(genes[,3]) < 0.1 | as.numeric(genes[,4]) < 0.1 | as.numeric(genes[,5]) < 0.1 | as.numeric(genes[,6]) < 0.1)
# removing the BCR and TCR genes from the comparison
remove <- vector()
for (i in 1:nrow(genes)) {
  if (isTRUE(str_contains(genes[i,1], "IG")) | isTRUE(str_contains(genes[i,1], "TRA")) | isTRUE(str_contains(genes[i,1], "TRB")) | isTRUE(str_contains(genes[i,1], "TRD")) | isTRUE(str_contains(genes[i,1], "TRG"))) {remove <- append(remove, i)}}
genes <- genes[-remove,]

# writing the table out
write.table(genes, "sig_panel_only_snv_fishers_test_results.tab", col.names = T, row.names = F, sep="\t", quote=F)
