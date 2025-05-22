# compute the empirical p-values for the genes of interest
# based on randomized CNV data

# list of genes of interest
goi <- read.table("cnv_genes_of_interest_cn_de_results.tab", header=T)

# deletion and amplification fisher's test results
del <- read.table("del_fishers_test_results.tab", header=T)
amp <- read.table("amp_fishers_test_results.tab", header=T)

# adding columns to the table that describe the chromosome and coordinate
newcols <- matrix(ncol=3, nrow=nrow(goi))
cnv <- read.table("collective_CNV_filt_annot_per_gene.tab", header=T)
check <- vector()
for (i in 1:nrow(goi)) {
  cnv2 <- subset(cnv, cnv$Gene_name==goi[i,1])
  if (length(unique(cnv2$Tx_start))==1 & length(unique(cnv2$Tx_end))==1) {
    newcols[i,1] <- cnv2[1,1]
    newcols[i,2] <- cnv2[1,9]
    newcols[i,3] <- cnv2[1,10]}
  else {
    newcols[i,1] <- cnv2[1,1]
    newcols[i,2] <- min(as.numeric(cnv2[,9]))
    newcols[i,3] <- max(as.numeric(cnv2[,10]))}}
goi <- cbind(goi[,1:2], newcols)

# adding columns describing from which analyses the gene arises from
newcols <- matrix(ncol=8, nrow=nrow(goi))
for (i in 1:nrow(goi)) {
  del2 <- subset(del, del$Gene_name==goi[i,1])
  if (isTRUE(as.numeric(del2[1,3]) < 0.1)) {newcols[i,1] <- "yes"}
  if (isTRUE(as.numeric(del2[1,5]) < 0.1)) {newcols[i,3] <- "yes"}
  if (isTRUE(as.numeric(del2[1,6]) < 0.1)) {newcols[i,5] <- "yes"}
  if (isTRUE(as.numeric(del2[1,7]) < 0.1)) {newcols[i,7] <- "yes"}
  amp2 <- subset(amp, amp$Gene_name==goi[i,1])
  if (isTRUE(as.numeric(amp2[1,3]) < 0.1)) {newcols[i,2] <- "yes"}
  if (isTRUE(as.numeric(amp2[1,5]) < 0.1)) {newcols[i,4] <- "yes"}
  if (isTRUE(as.numeric(amp2[1,6]) < 0.1)) {newcols[i,6] <- "yes"}
  if (isTRUE(as.numeric(amp2[1,7]) < 0.1)) {newcols[i,8] <- "yes"}}
colnames(newcols) <- c("d15_del_gene", "d15_amp_gene", "d29fm_del_gene", "d29fm_amp_gene", "d29ms_del_gene",	"d29ms_amp_gene",	"d79_del_gene",	"d79_amp_gene")
goi <- cbind(goi, newcols)

# response table
response <- read.table("sample_responses_cnv_cohort.txt", header=T)

# generating vectors of all cases in each response group
fast15 <- subset(response, response$d15_response=="fast")
fast15 <- as.vector(fast15[,1])
slow15 <- subset(response, response$d15_response=="slow")
slow15 <- as.vector(slow15[,1])
fast29 <- subset(response, response$d29_response=="fast")
fast29 <- as.vector(fast29[,1])
fastmed29 <- subset(response, response$d29_response=="fast" | response$d29_response=="intermediate")
fastmed29 <- as.vector(fastmed29[,1])
slow29 <- subset(response, response$d29_response=="slow")
slow29 <- as.vector(slow29[,1])
medslow29 <- subset(response, response$d29_response=="slow" | response$d29_response=="intermediate")
medslow29 <- as.vector(medslow29[,1])
fast79 <- subset(response, response$d79_response=="negative")
fast79 <- as.vector(fast79[,1])
slow79 <- subset(response, response$d79_response=="positive")
slow79 <- as.vector(slow79[,1])
samples15 <- subset(response, response$d15_response!="unknown")
samples15 <- as.vector(samples15[,1])
samples29 <- subset(response, response$d29_response!="unknown")
samples29 <- as.vector(samples29[,1])
samples79 <- subset(response, response$d79_response!="unknown")
samples79 <- as.vector(samples79[,1])

# saving the samples as a vector
samples <- as.vector(response[,1])

# WGS CNV results
cnv <- read.table("collective_coding_CNVs.tab", header=T)
# simplifying the table
getab <- cnv[,c(1:5,7,8)]

# panelseq and methylation array results
tab <- read.table("collective_CNV_filt_annot_per_gene.tab", header=T)
# simplifying the table
tab <- tab[,c(15,1:3,5,7,6)]

# combining the tables
colnames(getab) <- colnames(tab)
tab <- rbind(tab, getab)
# removing LOH events
tab <- subset(tab, tab[,5]!="LOH")
# correcting DUP into AMP
tab[,5] <- gsub("DUP", "AMP", tab[,5])
# including only coding consequences
tab <- subset(tab, tab[,7]=="loss" | tab[,7]=="nc_gene_loss" | tab[,7]=="coding_sequence_variant" | tab[,7]=="gain" | tab[,7]=="nc_gene_gain" | tab[,7]=="nc_gene_sequence_variant")

# subsetting the table to only include the genes of interest, and the present samples
tab <- subset(tab, tab[,6] %in% goi[,1])
tab <- unique(tab)
# adding columns to the table describing the response groups
newcols <- matrix(ncol=3, nrow=nrow(tab))
colnames(newcols) <- c("d15_response", "d29_response", "d79_response")
for (i in 1:nrow(tab)) {
  resp2 <- subset(response, response[,1]==tab[i,1])
  newcols[i,1] <- resp2[1,2]
  newcols[i,2] <- resp2[1,3]
  newcols[i,3] <- resp2[1,4]}
tab <- cbind(tab, newcols)

## starting with amplifications

# determining the genes to include
amps <- subset(goi, goi[,7]=="yes" | goi[,9]=="yes" | goi[,11]=="yes" | goi[,13]=="yes")
# checking the unique chromosomes
chr <- unique(amps[,3])

# generating a table into which the results can be included
res <- matrix(ncol=4, nrow=nrow(amps))
colnames(res) <- c("d15_epval", "d29fm_epval", "d29ms_epval", "d79_epval")
rownames(res) <- amps[,1]

# conducting the following steps for each chromosome and each gene, and each shuffle round
for (c in chr) {
  
  # determining the genes within the chromosomal region
  genes <- subset(amps, amps[,3]==c)
  genes <- as.vector(genes[,1])
  
  # conducting the following steps for each gene
  for (g in genes) {
    amps2 <- subset(amps, amps[,1]==g)
    # determining the gene start and end
    start <- amps2[1,4]
    end <- amps2[1,5]
    
    # starting with day 15 comparison
    if (isTRUE(amps2[1,7]=="yes")) {
      # conducting the following steps for each shuffle round, adding the response group difference to a collective table
      gene_dif <- matrix(ncol=1, nrow=10000)
      for (i in 1:10000) {
        # reading in the segment table for the chromosome for the present shuffle round
        if (c=="12") {seg <- read.table(paste0(c, "_amp/collective_", c, "p_AMPs_shuffle", i, ".bed"), header=F)}
        if (c=="21") {seg <- read.table(paste0(c, "_amp/collective_", c, "q_AMPs_shuffle", i, ".bed"), header=F)}
        if (c=="10") {seg <- read.table(paste0(c, "_amp/collective_", c, "_AMPs_shuffle", i, ".bed"), header=F)}
        # subsetting the table to include only the regions that overlap the present gene
        gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
        # randomizing the response groups 
        slow <- sample(samples15, length(slow15))
        fast <- setdiff(samples15, slow)
        # subsetting the table to slow and fast
        slow <- subset(gene_seg, gene_seg[,5] %in% slow)
        fast <- subset(gene_seg, gene_seg[,5] %in% fast)
        # adding the response group difference to the table
        slowfrac <- (length(unique(slow[,5]))/length(slow15))
        fastfrac <- (length(unique(fast[,5]))/length(fast15))
        gene_dif[i,1] <- (slowfrac-fastfrac)}
      # comparing to the true difference
      # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
      if (c=="12") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,4] <=15000000 & tab[,6]==g)}
      if (c=="21") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,3] >= 9400000 & tab[,6]==g)}
      if (c=="10") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,6]==g)}
      trueslow <- subset(tab2, tab2[,8]=="slow")
      truefast <- subset(tab2, tab2[,8]=="fast")
      trueslowfrac <- (length(unique(trueslow[,1]))/length(slow15))
      truefastfrac <- (length(unique(truefast[,1]))/length(fast15))
      truedif <- (trueslowfrac-truefastfrac)
      if (as.numeric(truedif) > 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
        res[g,"d15_epval"] <- (nrow(gene_dif_higher)/10000)}
      if (as.numeric(truedif) < 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
        res[g,"d15_epval"] <- (nrow(gene_dif_higher)/10000)}}
    
    # day 29 fastmed comparison
    if (isTRUE(amps2[1,9]=="yes")) {
      # conducting the following steps for each shuffle round, adding the response group difference to a collective table
      gene_dif <- matrix(ncol=1, nrow=10000)
      for (i in 1:10000) {
        # reading in the segment table for the chromosome for the present shuffle round
        if (c=="12") {seg <- read.table(paste0(c, "_amp/collective_", c, "p_AMPs_shuffle", i, ".bed"), header=F)}
        if (c=="21") {seg <- read.table(paste0(c, "_amp/collective_", c, "q_AMPs_shuffle", i, ".bed"), header=F)}
        if (c=="10") {seg <- read.table(paste0(c, "_amp/collective_", c, "_AMPs_shuffle", i, ".bed"), header=F)}
        # subsetting the table to include only the regions that overlap the present gene
        gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
        # randomizing the response groups 
        slow <- sample(samples29, length(slow29))
        fast <- setdiff(samples29, slow)
        # subsetting the table to slow and fast
        slow <- subset(gene_seg, gene_seg[,5] %in% slow)
        fast <- subset(gene_seg, gene_seg[,5] %in% fast)
        # adding the response group difference to the table
        slowfrac <- (length(unique(slow[,5]))/length(slow29))
        fastfrac <- (length(unique(fast[,5]))/length(fastmed29))
        gene_dif[i,1] <- (slowfrac-fastfrac)}
      # comparing to the true difference
      # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
      if (c=="12") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,4] <=15000000 & tab[,6]==g)}
      if (c=="21") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,3] >= 9400000 & tab[,6]==g)}
      if (c=="10") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,6]==g)}
      trueslow <- subset(tab2, tab2[,9]=="slow")
      truefast <- subset(tab2, tab2[,9]=="fast" | tab2[,9]=="intermediate")
      trueslowfrac <- (length(unique(trueslow[,1]))/length(slow29))
      truefastfrac <- (length(unique(truefast[,1]))/length(fastmed29))
      truedif <- (trueslowfrac-truefastfrac)
      if (as.numeric(truedif) > 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
        res[g,"d29fm_epval"] <- (nrow(gene_dif_higher)/10000)}
      if (as.numeric(truedif) < 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
        res[g,"d29fm_epval"] <- (nrow(gene_dif_higher)/10000)}}
    
    # day 29 medslow comparison
    if (isTRUE(amps2[1,11]=="yes")) {
      # conducting the following steps for each shuffle round, adding the response group difference to a collective table
      gene_dif <- matrix(ncol=1, nrow=10000)
      for (i in 1:10000) {
        # reading in the segment table for the chromosome for the present shuffle round
        if (c=="12") {seg <- read.table(paste0(c, "_amp/collective_", c, "p_AMPs_shuffle", i, ".bed"), header=F)}
        if (c=="21") {seg <- read.table(paste0(c, "_amp/collective_", c, "q_AMPs_shuffle", i, ".bed"), header=F)}
        if (c=="10") {seg <- read.table(paste0(c, "_amp/collective_", c, "_AMPs_shuffle", i, ".bed"), header=F)}
        # subsetting the table to include only the regions that overlap the present gene
        gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
        # randomizing the response groups 
        slow <- sample(samples29, length(medslow29))
        fast <- setdiff(samples29, slow)
        # subsetting the table to slow and fast
        slow <- subset(gene_seg, gene_seg[,5] %in% slow)
        fast <- subset(gene_seg, gene_seg[,5] %in% fast)
        # adding the response group difference to the table
        slowfrac <- (length(unique(slow[,5]))/length(medslow29))
        fastfrac <- (length(unique(fast[,5]))/length(fast29))
        gene_dif[i,1] <- (slowfrac-fastfrac)}
      # comparing to the true difference
      # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
      if (c=="12") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,4] <=15000000 & tab[,6]==g)}
      if (c=="21") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,3] >= 9400000 & tab[,6]==g)}
      if (c=="10") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,6]==g)}
      trueslow <- subset(tab2, tab2[,9]=="slow" | tab2[,9]=="intermediate")
      truefast <- subset(tab2, tab2[,9]=="fast")
      trueslowfrac <- (length(unique(trueslow[,1]))/length(medslow29))
      truefastfrac <- (length(unique(truefast[,1]))/length(fast29))
      truedif <- (trueslowfrac-truefastfrac)
      if (as.numeric(truedif) > 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
        res[g,"d29ms_epval"] <- (nrow(gene_dif_higher)/10000)}
      if (as.numeric(truedif) < 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
        res[g,"d29ms_epval"] <- (nrow(gene_dif_higher)/10000)}}
    
    # day 79 comparison
    if (isTRUE(amps2[1,13]=="yes")) {
      # conducting the following steps for each shuffle round, adding the response group difference to a collective table
      gene_dif <- matrix(ncol=1, nrow=10000)
      for (i in 1:10000) {
        # reading in the segment table for the chromosome for the present shuffle round
        if (c=="12") {seg <- read.table(paste0(c, "_amp/collective_", c, "p_AMPs_shuffle", i, ".bed"), header=F)}
        if (c=="21") {seg <- read.table(paste0(c, "_amp/collective_", c, "q_AMPs_shuffle", i, ".bed"), header=F)}
        if (c=="10") {seg <- read.table(paste0(c, "_amp/collective_", c, "_AMPs_shuffle", i, ".bed"), header=F)}
        # subsetting the table to include only the regions that overlap the present gene
        gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
        # randomizing the response groups 
        slow <- sample(samples79, length(slow79))
        fast <- setdiff(samples79, slow)
        # subsetting the table to slow and fast
        slow <- subset(gene_seg, gene_seg[,5] %in% slow)
        fast <- subset(gene_seg, gene_seg[,5] %in% fast)
        # adding the response group difference to the table
        slowfrac <- (length(unique(slow[,5]))/length(slow79))
        fastfrac <- (length(unique(fast[,5]))/length(fast79))
        gene_dif[i,1] <- (slowfrac-fastfrac)}
      # comparing to the true difference
      # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
      if (c=="12") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,4] <= 15000000 & tab[,6]==g)}
      if (c=="21") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,3] >= 9400000 & tab[,6]==g)}
      if (c=="10") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="AMP" & tab[,6]==g)}
      trueslow <- subset(tab2, tab2[,10]=="positive")
      truefast <- subset(tab2, tab2[,10]=="negative")
      trueslowfrac <- (length(unique(trueslow[,1]))/length(slow79))
      truefastfrac <- (length(unique(truefast[,1]))/length(fast79))
      truedif <- (trueslowfrac-truefastfrac)
      if (as.numeric(truedif) > 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
        res[g,"d79_epval"] <- (nrow(gene_dif_higher)/10000)}
      if (as.numeric(truedif) < 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
        res[g,"d79_epval"] <- (nrow(gene_dif_higher)/10000)}}}}

# writing the table out
res <- cbind(rownames(res), res)
colnames(res)[1] <- c("gene")
write.table(res, "collective_amp_epvals.tab", sep="\t", quote=F, col.names=T, row.names=F)

## continuing with deletions

# determining the genes to include
dels <- subset(goi, goi[,6]=="yes" | goi[,8]=="yes" | goi[,10]=="yes" | goi[,12]=="yes")
# checking the unique chromosomes
chr <- unique(dels[,3])

# generating a table into which the results can be included
res <- matrix(ncol=4, nrow=nrow(dels))
colnames(res) <- c("d15_epval", "d29fm_epval", "d29ms_epval", "d79_epval")
rownames(res) <- dels[,1]

# conducting the following steps for each chromosome and each gene, and each shuffle round
for (c in chr) {
  
  # determining the genes within the chromosomal region
  genes <- subset(dels, dels[,3]==c)
  genes <- as.vector(genes[,1])
  
  # conducting the following steps for each gene
  for (g in genes) {
    dels2 <- subset(dels, dels[,1]==g)
    # determining the gene start and end
    start <- dels2[1,4]
    end <- dels2[1,5]
    
    # starting with day 15 comparison
    if (isTRUE(dels2[1,6]=="yes")) {
      # conducting the following steps for each shuffle round, adding the response group difference to a collective table
      gene_dif <- matrix(ncol=1, nrow=10000)
      for (i in 1:10000) {
        # reading in the segment table for the chromosome for the present shuffle round
        if (c=="12") {seg <- read.table(paste0(c, "_del/collective_", c, "p_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="15") {seg <- read.table(paste0(c, "_del/collective_", c, "q_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="6") {seg <- read.table(paste0(c, "_del/collective_", c, "q_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="3") {seg <- read.table(paste0(c, "_del/collective_", c, "_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="9") {seg <- read.table(paste0(c, "_del/collective_", c, "_DELs_shuffle", i, ".bed"), header=F)}
        # subsetting the table to include only the regions that overlap the present gene
        gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
        # randomizing the response groups 
        slow <- sample(samples15, length(slow15))
        fast <- setdiff(samples15, slow)
        # subsetting the table to slow and fast
        slow <- subset(gene_seg, gene_seg[,5] %in% slow)
        fast <- subset(gene_seg, gene_seg[,5] %in% fast)
        # adding the response group difference to the table
        slowfrac <- (length(unique(slow[,5]))/length(slow15))
        fastfrac <- (length(unique(fast[,5]))/length(fast15))
        gene_dif[i,1] <- (slowfrac-fastfrac)}
      # comparing to the true difference
      # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
      if (c=="12") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,4] <= 35000000 & tab[,6]==g)}
      if (c=="15") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,3] >= 30000000 & tab[,6]==g)}
      if (c=="6") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,3] >= 60000000 & tab[,6]==g)}
      if (c=="3") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,6]==g)}
      if (c=="9") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,6]==g)}
      trueslow <- subset(tab2, tab2[,8]=="slow")
      truefast <- subset(tab2, tab2[,8]=="fast")
      trueslowfrac <- (length(unique(trueslow[,1]))/length(slow15))
      truefastfrac <- (length(unique(truefast[,1]))/length(fast15))
      truedif <- (trueslowfrac-truefastfrac)
      if (as.numeric(truedif) > 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
        res[g,"d15_epval"] <- (nrow(gene_dif_higher)/10000)}
      if (as.numeric(truedif) < 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
        res[g,"d15_epval"] <- (nrow(gene_dif_higher)/10000)}}
    
    # day 29 fastmed comparison
    if (isTRUE(dels2[1,8]=="yes")) {
      # conducting the following steps for each shuffle round, adding the response group difference to a collective table
      gene_dif <- matrix(ncol=1, nrow=10000)
      for (i in 1:10000) {
        # reading in the segment table for the chromosome for the present shuffle round
        if (c=="12") {seg <- read.table(paste0(c, "_del/collective_", c, "p_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="15") {seg <- read.table(paste0(c, "_del/collective_", c, "q_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="6") {seg <- read.table(paste0(c, "_del/collective_", c, "q_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="3") {seg <- read.table(paste0(c, "_del/collective_", c, "_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="9") {seg <- read.table(paste0(c, "_del/collective_", c, "_DELs_shuffle", i, ".bed"), header=F)}
        # subsetting the table to include only the regions that overlap the present gene
        gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
        # randomizing the response groups 
        slow <- sample(samples29, length(slow29))
        fast <- setdiff(samples29, slow)
        # subsetting the table to slow and fast
        slow <- subset(gene_seg, gene_seg[,5] %in% slow)
        fast <- subset(gene_seg, gene_seg[,5] %in% fast)
        # adding the response group difference to the table
        slowfrac <- (length(unique(slow[,5]))/length(slow29))
        fastfrac <- (length(unique(fast[,5]))/length(fastmed29))
        gene_dif[i,1] <- (slowfrac-fastfrac)}
      # comparing to the true difference
      # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
      if (c=="12") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,4] <= 35000000 & tab[,6]==g)}
      if (c=="15") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,3] >= 30000000 & tab[,6]==g)}
      if (c=="6") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,3] >= 60000000 & tab[,6]==g)}
      if (c=="3") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,6]==g)}
      if (c=="9") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,6]==g)}
      trueslow <- subset(tab2, tab2[,9]=="slow")
      truefast <- subset(tab2, tab2[,9]=="fast" | tab2[,9]=="intermediate")
      trueslowfrac <- (length(unique(trueslow[,1]))/length(slow29))
      truefastfrac <- (length(unique(truefast[,1]))/length(fastmed29))
      truedif <- (trueslowfrac-truefastfrac)
      if (as.numeric(truedif) > 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
        res[g,"d29fm_epval"] <- (nrow(gene_dif_higher)/10000)}
      if (as.numeric(truedif) < 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
        res[g,"d29fm_epval"] <- (nrow(gene_dif_higher)/10000)}}
    
    # day 29 medslow comparison
    if (isTRUE(dels2[1,10]=="yes")) {
      # conducting the following steps for each shuffle round, adding the response group difference to a collective table
      gene_dif <- matrix(ncol=1, nrow=10000)
      for (i in 1:10000) {
        # reading in the segment table for the chromosome for the present shuffle round
        if (c=="12") {seg <- read.table(paste0(c, "_del/collective_", c, "p_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="15") {seg <- read.table(paste0(c, "_del/collective_", c, "q_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="6") {seg <- read.table(paste0(c, "_del/collective_", c, "q_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="3") {seg <- read.table(paste0(c, "_del/collective_", c, "_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="9") {seg <- read.table(paste0(c, "_del/collective_", c, "_DELs_shuffle", i, ".bed"), header=F)}
        # subsetting the table to include only the regions that overlap the present gene
        gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
        # randomizing the response groups 
        slow <- sample(samples29, length(medslow29))
        fast <- setdiff(samples29, slow)
        # subsetting the table to slow and fast
        slow <- subset(gene_seg, gene_seg[,5] %in% slow)
        fast <- subset(gene_seg, gene_seg[,5] %in% fast)
        # adding the response group difference to the table
        slowfrac <- (length(unique(slow[,5]))/length(medslow29))
        fastfrac <- (length(unique(fast[,5]))/length(fast29))
        gene_dif[i,1] <- (slowfrac-fastfrac)}
      # comparing to the true difference
      # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
      if (c=="12") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,4] <= 35000000 & tab[,6]==g)}
      if (c=="15") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,3] >= 30000000 & tab[,6]==g)}
      if (c=="6") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,3] >= 60000000 & tab[,6]==g)}
      if (c=="3") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,6]==g)}
      if (c=="9") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,6]==g)}
      trueslow <- subset(tab2, tab2[,9]=="slow" | tab2[,9]=="intermediate")
      truefast <- subset(tab2, tab2[,9]=="fast")
      trueslowfrac <- (length(unique(trueslow[,1]))/length(medslow29))
      truefastfrac <- (length(unique(truefast[,1]))/length(fast29))
      truedif <- (trueslowfrac-truefastfrac)
      if (as.numeric(truedif) > 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
        res[g,"d29ms_epval"] <- (nrow(gene_dif_higher)/10000)}
      if (as.numeric(truedif) < 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
        res[g,"d29ms_epval"] <- (nrow(gene_dif_higher)/10000)}}
    
    # day 79 comparison
    if (isTRUE(dels2[1,12]=="yes")) {
      # conducting the following steps for each shuffle round, adding the response group difference to a collective table
      gene_dif <- matrix(ncol=1, nrow=10000)
      for (i in 1:10000) {
        # reading in the segment table for the chromosome for the present shuffle round
        if (c=="12") {seg <- read.table(paste0(c, "_del/collective_", c, "p_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="15") {seg <- read.table(paste0(c, "_del/collective_", c, "q_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="6") {seg <- read.table(paste0(c, "_del/collective_", c, "q_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="3") {seg <- read.table(paste0(c, "_del/collective_", c, "_DELs_shuffle", i, ".bed"), header=F)}
        if (c=="9") {seg <- read.table(paste0(c, "_del/collective_", c, "_DELs_shuffle", i, ".bed"), header=F)}
        # subsetting the table to include only the regions that overlap the present gene
        gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
        # randomizing the response groups 
        slow <- sample(samples79, length(slow79))
        fast <- setdiff(samples79, slow)
        # subsetting the table to slow and fast
        slow <- subset(gene_seg, gene_seg[,5] %in% slow)
        fast <- subset(gene_seg, gene_seg[,5] %in% fast)
        # adding the response group difference to the table
        slowfrac <- (length(unique(slow[,5]))/length(slow79))
        fastfrac <- (length(unique(fast[,5]))/length(fast79))
        gene_dif[i,1] <- (slowfrac-fastfrac)}
      # comparing to the true difference
      # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
      if (c=="12") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,4] <= 35000000 & tab[,6]==g)}
      if (c=="15") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,3] >= 30000000 & tab[,6]==g)}
      if (c=="6") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,3] >= 60000000 & tab[,6]==g)}
      if (c=="3") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,6]==g)}
      if (c=="9") {tab2 <- subset(tab, tab[,2]==c & tab[,5]=="DEL" & tab[,6]==g)}
      trueslow <- subset(tab2, tab2[,10]=="positive")
      truefast <- subset(tab2, tab2[,10]=="negative")
      trueslowfrac <- (length(unique(trueslow[,1]))/length(slow79))
      truefastfrac <- (length(unique(truefast[,1]))/length(fast79))
      truedif <- (trueslowfrac-truefastfrac)
      if (as.numeric(truedif) > 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
        res[g,"d79_epval"] <- (nrow(gene_dif_higher)/10000)}
      if (as.numeric(truedif) < 0) {
        gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
        res[g,"d79_epval"] <- (nrow(gene_dif_higher)/10000)}}}}

# writing the table out
res <- cbind(rownames(res), res)
colnames(res)[1] <- c("gene")
write.table(res, "collective_del_epvals.tab", sep="\t", quote=F, col.names=T, row.names=F)


## conducting the analysis for the control region, i.e., chr16 amplifications

# reading in the amplification fisher's test results
amp <- read.table("amp_fishers_test_results.tab", header=T)
amp <- subset(amp, amp[,2]=="16")

# adding columns to the table that describe the chromosome and coordinate
newcols <- matrix(ncol=3, nrow=nrow(amp))
cnv <- read.table("collective_CNV_filt_annot_per_gene.tab", header=T)
check <- vector()
for (i in 1:nrow(amp)) {
  cnv2 <- subset(cnv, cnv$Gene_name==amp[i,1])
  if (length(unique(cnv2$Tx_start))==1 & length(unique(cnv2$Tx_end))==1) {
    newcols[i,1] <- cnv2[1,1]
    newcols[i,2] <- cnv2[1,9]
    newcols[i,3] <- cnv2[1,10]}
  else {
    newcols[i,1] <- cnv2[1,1]
    newcols[i,2] <- min(as.numeric(cnv2[,9]))
    newcols[i,3] <- max(as.numeric(cnv2[,10]))}}
goi <- cbind(amp[,1:2], newcols)
# sampling 10 genes from the table to include in the tests
rows <- sample(1:nrow(goi), 10)
goi <- goi[rows,]

# generating a table into which the results can be included
res <- matrix(ncol=4, nrow=nrow(goi))
colnames(res) <- c("d15_epval", "d29fm_epval", "d29ms_epval", "d79_epval")
rownames(res) <- goi[,1]

# conducting the following steps for each chromosome and each gene, and each shuffle round
genes <- goi[,1]
for (g in genes) {
  
  amps2 <- subset(goi, goi[,1]==g)
  # determining the gene start and end
  start <- amps2[1,4]
  end <- amps2[1,5]
  
  # starting with day 15 comparison
  # conducting the following steps for each shuffle round, adding the response group difference to a collective table
  gene_dif <- matrix(ncol=1, nrow=10000)
  for (i in 1:10000) {
    # reading in the segment table for the chromosome for the present shuffle round
    seg <- read.table(paste0("16_amp/collective_16_AMPs_shuffle", i, ".bed"), header=F)
    # subsetting the table to include only the regions that overlap the present gene
    gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
    # randomizing the response groups 
    slow <- sample(samples15, length(slow15))
    fast <- setdiff(samples15, slow)
    # subsetting the table to slow and fast
    slow <- subset(gene_seg, gene_seg[,5] %in% slow)
    fast <- subset(gene_seg, gene_seg[,5] %in% fast)
    # adding the response group difference to the table
    slowfrac <- (length(unique(slow[,5]))/length(slow15))
    fastfrac <- (length(unique(fast[,5]))/length(fast15))
    gene_dif[i,1] <- (slowfrac-fastfrac)}
  # comparing to the true difference
  # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
  tab2 <- subset(tab, tab[,5]=="AMP" & tab[,6]==g)
  trueslow <- subset(tab2, tab2[,8]=="slow")
  truefast <- subset(tab2, tab2[,8]=="fast")
  trueslowfrac <- (length(unique(trueslow[,1]))/length(slow15))
  truefastfrac <- (length(unique(truefast[,1]))/length(fast15))
  truedif <- (trueslowfrac-truefastfrac)
  if (as.numeric(truedif) > 0) {
    gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
    res[g,"d15_epval"] <- (nrow(gene_dif_higher)/10000)}
  if (as.numeric(truedif) < 0) {
    gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
    res[g,"d15_epval"] <- (nrow(gene_dif_higher)/10000)}
  
  # day 29 fastmed comparison
  # conducting the following steps for each shuffle round, adding the response group difference to a collective table
  gene_dif <- matrix(ncol=1, nrow=10000)
  for (i in 1:10000) {
    # reading in the segment table for the chromosome for the present shuffle round
    seg <- read.table(paste0("16_amp/collective_16_AMPs_shuffle", i, ".bed"), header=F)
    # subsetting the table to include only the regions that overlap the present gene
    gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
    # randomizing the response groups 
    slow <- sample(samples29, length(slow29))
    fast <- setdiff(samples29, slow)
    # subsetting the table to slow and fast
    slow <- subset(gene_seg, gene_seg[,5] %in% slow)
    fast <- subset(gene_seg, gene_seg[,5] %in% fast)
    # adding the response group difference to the table
    slowfrac <- (length(unique(slow[,5]))/length(slow29))
    fastfrac <- (length(unique(fast[,5]))/length(fastmed29))
    gene_dif[i,1] <- (slowfrac-fastfrac)}
  # comparing to the true difference
  # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
  tab2 <- subset(tab, tab[,5]=="AMP" & tab[,6]==g)
  trueslow <- subset(tab2, tab2[,9]=="slow")
  truefast <- subset(tab2, tab2[,9]=="fast" | tab2[,9]=="intermediate")
  trueslowfrac <- (length(unique(trueslow[,1]))/length(slow29))
  truefastfrac <- (length(unique(truefast[,1]))/length(fastmed29))
  truedif <- (trueslowfrac-truefastfrac)
  if (as.numeric(truedif) > 0) {
    gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
    res[g,"d29fm_epval"] <- (nrow(gene_dif_higher)/10000)}
  if (as.numeric(truedif) < 0) {
    gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
    res[g,"d29fm_epval"] <- (nrow(gene_dif_higher)/10000)}
  
  # day 29 medslow comparison
  # conducting the following steps for each shuffle round, adding the response group difference to a collective table
  gene_dif <- matrix(ncol=1, nrow=10000)
  for (i in 1:10000) {
    # reading in the segment table for the chromosome for the present shuffle round
    seg <- read.table(paste0("16_amp/collective_16_AMPs_shuffle", i, ".bed"), header=F)
    # subsetting the table to include only the regions that overlap the present gene
    gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
    # randomizing the response groups 
    slow <- sample(samples29, length(medslow29))
    fast <- setdiff(samples29, slow)
    # subsetting the table to slow and fast
    slow <- subset(gene_seg, gene_seg[,5] %in% slow)
    fast <- subset(gene_seg, gene_seg[,5] %in% fast)
    # adding the response group difference to the table
    slowfrac <- (length(unique(slow[,5]))/length(medslow29))
    fastfrac <- (length(unique(fast[,5]))/length(fast29))
    gene_dif[i,1] <- (slowfrac-fastfrac)}
  # comparing to the true difference
  # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
  tab2 <- subset(tab, tab[,5]=="AMP" & tab[,6]==g)
  trueslow <- subset(tab2, tab2[,9]=="slow" | tab2[,9]=="intermediate")
  truefast <- subset(tab2, tab2[,9]=="fast")
  trueslowfrac <- (length(unique(trueslow[,1]))/length(medslow29))
  truefastfrac <- (length(unique(truefast[,1]))/length(fast29))
  truedif <- (trueslowfrac-truefastfrac)
  if (as.numeric(truedif) > 0) {
    gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
    res[g,"d29ms_epval"] <- (nrow(gene_dif_higher)/10000)}
  if (as.numeric(truedif) < 0) {
    gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
    res[g,"d29ms_epval"] <- (nrow(gene_dif_higher)/10000)}
  
  # day 79 comparison
  # conducting the following steps for each shuffle round, adding the response group difference to a collective table
  gene_dif <- matrix(ncol=1, nrow=10000)
  for (i in 1:10000) {
    # reading in the segment table for the chromosome for the present shuffle round
    seg <- read.table(paste0("16_amp/collective_16_AMPs_shuffle", i, ".bed"), header=F)
    # subsetting the table to include only the regions that overlap the present gene
    gene_seg <- subset(seg, seg[,2]<start & seg[,3]>end | seg[,2]<start & seg[,3]>start & seg[,3]<end | seg[,2]>start & seg[,2]<end & seg[,3]>end | seg[,2]>start & seg[,2]<end & seg[,3]<end)
    # randomizing the response groups 
    slow <- sample(samples79, length(slow79))
    fast <- setdiff(samples79, slow)
    # subsetting the table to slow and fast
    slow <- subset(gene_seg, gene_seg[,5] %in% slow)
    fast <- subset(gene_seg, gene_seg[,5] %in% fast)
    # adding the response group difference to the table
    slowfrac <- (length(unique(slow[,5]))/length(slow79))
    fastfrac <- (length(unique(fast[,5]))/length(fast79))
    gene_dif[i,1] <- (slowfrac-fastfrac)}
  # comparing to the true difference
  # subsetting the gene table to the present gene hits, using the same coordinate cutoffs as in the shuffles
  tab2 <- subset(tab, tab[,5]=="AMP" & tab[,6]==g)
  trueslow <- subset(tab2, tab2[,10]=="positive")
  truefast <- subset(tab2, tab2[,10]=="negative")
  trueslowfrac <- (length(unique(trueslow[,1]))/length(slow79))
  truefastfrac <- (length(unique(truefast[,1]))/length(fast79))
  truedif <- (trueslowfrac-truefastfrac)
  if (as.numeric(truedif) > 0) {
    gene_dif_higher <- subset(gene_dif, gene_dif[,1] >= truedif)
    res[g,"d79_epval"] <- (nrow(gene_dif_higher)/10000)}
  if (as.numeric(truedif) < 0) {
    gene_dif_higher <- subset(gene_dif, gene_dif[,1] <= truedif)
    res[g,"d79_epval"] <- (nrow(gene_dif_higher)/10000)}}

# writing the table out
res <- cbind(rownames(res), res)
colnames(res)[1] <- c("gene")
write.table(res, "collective_16amp_ctrl_epvals.tab", sep="\t", quote=F, col.names=T, row.names=F)
