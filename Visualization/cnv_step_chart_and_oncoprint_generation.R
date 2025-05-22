# generate the step charts describing the fraction of cases with deletions/amplifications
# in each chromosomal region per response group
library(ggplot2)
library(gridExtra)
library(scales)

# CNVs
tab <- read.table("collective_CNV_filt_annot_per_event.tab", sep="\t", header=T)
tab2 <- read.table("collective_CNV_annotated_only_full.tab")
# combining the tables
tab2 <- tab2[,1:5]
tab <- tab[,c(8,1:3,5)]
colnames(tab2) <- colnames(tab)
tab <- rbind(tab, tab2)
# removing LOH events
tab <- subset(tab, tab[,5]!="LOH")
# correcting DUP to AMP
tab[,5] <- gsub("DUP", "AMP", tab[,5])

# response table
resp <- read.table("sample_responses_cnv_cohort.txt", header=T)

# adding columns describing the response groups
newcols <- matrix(ncol=3, nrow=nrow(tab))
for (i in 1:nrow(tab)) {
  resp2 <- subset(resp, resp[,1]==tab[i,1])
  newcols[i,1] <- resp2[1,2]
  newcols[i,2] <- resp2[1,3]
  newcols[i,3] <- resp2[1,4]}
colnames(newcols) <- c("d15_response", "d29_response", "d79_response")
tab <- cbind(tab, newcols)

# generating vectors of all cases in each resp group
fast15 <- subset(resp, resp$d15_response=="fast")
fast15 <- as.vector(fast15[,1])
slow15 <- subset(resp, resp$d15_response=="slow")
slow15 <- as.vector(slow15[,1])
fast29 <- subset(resp, resp$d29_response=="fast")
fast29 <- as.vector(fast29[,1])
fastmed29 <- subset(resp, resp$d29_response=="fast" | resp$d29_response=="intermediate")
fastmed29 <- as.vector(fastmed29[,1])
slow29 <- subset(resp, resp$d29_response=="slow")
slow29 <- as.vector(slow29[,1])
medslow29 <- subset(resp, resp$d29_response=="slow" | resp$d29_response=="intermediate")
medslow29 <- as.vector(medslow29[,1])
fast79 <- subset(resp, resp$d79_response=="negative")
fast79 <- as.vector(fast79[,1])
slow79 <- subset(resp, resp$d79_response=="positive")
slow79 <- as.vector(slow79[,1])
samples15 <- subset(resp, resp$d15_response!="unknown")
samples15 <- as.vector(samples15[,1])
samples29 <- subset(resp, resp$d29_response!="unknown")
samples29 <- as.vector(samples29[,1])
samples79 <- subset(resp, resp$d79_response!="unknown")
samples79 <- as.vector(samples79[,1])

## CNV step chart generation

# conducting the following steps for each somatic chromosome
for (chr in c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)) {
  
  # subsetting the chromosome sizes table to only include the present chromosome
  size <- subset(sizes, sizes[,1]==paste0("chr", chr))
  
  # separating the table into deletions and amplifications affecting the present chromosome
  del <- subset(tab, tab[,2]==chr & tab[,5]=="DEL")
  amp <- subset(tab, tab[,2]==chr & tab[,5]=="AMP")
  
  # deletions, and day 15 comparison
  
  # focusing on the fast group
  fast <- subset(del, del[,6]=="fast")
  
  # conducting the next steps, if there are events reported
  if (nrow(fast)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(fast[,3])
    coords <- append(coords, fast[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the deletion fraction affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a deletion affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      delrows <- vector()
      for (r in 1:nrow(fast)) {
        if (as.numeric(fast[r,3])==as.numeric(cn_segments[i,2])) {delrows <- append(delrows, r)}
        if (as.numeric(fast[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(fast[r,4])>cn_segments[i,2]) {delrows <- append(delrows, r)}}
      del2 <- fast[delrows,]
      cn_segments[i,4] <- (length(unique(del2[,1]))/length(fast15))}
    
    # generating a data frame suitable for plot generation
    del_frac <- as.vector(cn_segments[,4])
    del_frac <- append(del_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    fast_df <- data.frame(del_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d15_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(fast)==0) {
    # generating a data frame suitable for plot generation
    del_frac <- c(0,0)
    coords <- append(0, size[1,2])
    fast_df <- data.frame(del_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d15_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # focusing on the slow group
  slow <- subset(del, del[,6]=="slow")
  
  # conducting the next steps, if there are events reported
  if (nrow(slow)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(slow[,3])
    coords <- append(coords, slow[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the Deletion fraction affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a deletion affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      delrows <- vector()
      for (r in 1:nrow(slow)) {
        if (as.numeric(slow[r,3])==as.numeric(cn_segments[i,2])) {delrows <- append(delrows, r)}
        if (as.numeric(slow[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(slow[r,4])>cn_segments[i,2]) {delrows <- append(delrows, r)}}
      del2 <- slow[delrows,]
      cn_segments[i,4] <- (length(unique(del2[,1]))/length(slow15))}
    
    # generating a data frame suitable for plot generation
    del_frac <- as.vector(cn_segments[,4])
    del_frac <- append(del_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    slow_df <- data.frame(del_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d15_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(slow)==0) {
    # generating a data frame suitable for plot generation
    del_frac <- c(0,0)
    coords <- append(0, size[1,2])
    slow_df <- data.frame(del_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d15_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # saving the plot
  p1 <- ggplot() + geom_step(data = fast_df, aes(x = coords, y = del_frac), colour="#a1d76a", size=0.3) + geom_step(data = slow_df, aes(x = coords, y = del_frac), colour="#b2182b", size=0.3) + ylim(0,0.1) + ggtitle(paste0("Chromosome ", chr, " deletion fractions by day 15 response"))  + xlab("Coordinates") + ylab("Deletion fraction") + scale_x_continuous(labels = label_comma()) + theme(panel.background=element_rect(fill="azure2", colour="azure2",size=0.5,linetype="solid"), text=element_text(size=5), plot.title=element_text(size=5))
  
  # deletions, and day 29 comparison
  
  # focusing on the fast/intermediate group
  fast <- subset(del, del[,7]=="fast" | del[,7]=="intermediate")
  
  # conducting the next steps, if there are events reported
  if (nrow(fast)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(fast[,3])
    coords <- append(coords, fast[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the Deletion fraction affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a deletion affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      delrows <- vector()
      for (r in 1:nrow(fast)) {
        if (as.numeric(fast[r,3])==as.numeric(cn_segments[i,2])) {delrows <- append(delrows, r)}
        if (as.numeric(fast[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(fast[r,4])>cn_segments[i,2]) {delrows <- append(delrows, r)}}
      del2 <- fast[delrows,]
      cn_segments[i,4] <- (length(unique(del2[,1]))/(length(fast29)+length(med29)))}
    # generating a data frame suitable for plot generation
    del_frac <- as.vector(cn_segments[,4])
    del_frac <- append(del_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    fast_df <- data.frame(del_frac, coords)
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d29_fastmed.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(fast)==0) {
    # generating a data frame suitable for plot generation
    del_frac <- c(0,0)
    coords <- append(0, size[1,2])
    fast_df <- data.frame(del_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d29_fastmed.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # focusing on the slow group
  slow <- subset(del, del[,7]=="slow")
  
  # conducting the next steps, if there are events reported
  if (nrow(slow)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(slow[,3])
    coords <- append(coords, slow[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the Deletion fraction affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a deletion affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      delrows <- vector()
      for (r in 1:nrow(slow)) {
        if (as.numeric(slow[r,3])==as.numeric(cn_segments[i,2])) {delrows <- append(delrows, r)}
        if (as.numeric(slow[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(slow[r,4])>cn_segments[i,2]) {delrows <- append(delrows, r)}}
      del2 <- slow[delrows,]
      cn_segments[i,4] <- (length(unique(del2[,1]))/length(slow29))}
    
    # generating a data frame suitable for plot generation
    del_frac <- as.vector(cn_segments[,4])
    del_frac <- append(del_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    slow_df <- data.frame(del_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d29_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(slow)==0) {
    # generating a data frame suitable for plot generation
    del_frac <- c(0,0)
    coords <- append(0, size[1,2])
    slow_df <- data.frame(del_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d29_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # saving the plot
  p2 <- ggplot() + geom_step(data = fast_df, aes(x = coords, y = del_frac), colour="#a1d76a", size=0.3) + geom_step(data = slow_df, aes(x = coords, y = del_frac), colour="#b2182b", size=0.3) + ylim(0,0.1) + ggtitle(paste0("Chromosome ", chr, " deletion fractions by day 29 response (fast+intermediate grouped)")) + xlab("Coordinates") + ylab("Deletion fraction") + scale_x_continuous(labels = label_comma()) + theme(panel.background=element_rect(fill="azure2", colour="azure2",size=0.5,linetype="solid"), text=element_text(size=5), plot.title=element_text(size=5))
  
  # focusing on the fast group only
  fast <- subset(del, del[,7]=="fast")
  
  # conducting the next steps, if there are events reported
  if (nrow(fast)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(fast[,3])
    coords <- append(coords, fast[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the Deletion fraction affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a deletion affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      delrows <- vector()
      for (r in 1:nrow(fast)) {
        if (as.numeric(fast[r,3])==as.numeric(cn_segments[i,2])) {delrows <- append(delrows, r)}
        if (as.numeric(fast[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(fast[r,4])>cn_segments[i,2]) {delrows <- append(delrows, r)}}
      del2 <- fast[delrows,]
      cn_segments[i,4] <- (length(unique(del2[,1]))/(length(fast29)))}
    
    # generating a data frame suitable for plot generation
    del_frac <- as.vector(cn_segments[,4])
    del_frac <- append(del_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    fast_df <- data.frame(del_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d29_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(fast)==0) {
    # generating a data frame suitable for plot generation
    del_frac <- c(0,0)
    coords <- append(0, size[1,2])
    fast_df <- data.frame(del_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d29_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # focusing on the intermediate/slow group
  slow <- subset(del, del[,7]=="slow" | del[,7]=="intermediate")
  
  # conducting the next steps, if there are events reported
  if (nrow(slow)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(slow[,3])
    coords <- append(coords, slow[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the Deletion fraction affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
   
    # determining the cases with a deletion affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      delrows <- vector()
      for (r in 1:nrow(slow)) {
        if (as.numeric(slow[r,3])==as.numeric(cn_segments[i,2])) {delrows <- append(delrows, r)}
        if (as.numeric(slow[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(slow[r,4])>cn_segments[i,2]) {delrows <- append(delrows, r)}}
      del2 <- slow[delrows,]
      cn_segments[i,4] <- (length(unique(del2[,1]))/(length(slow29)+length(med29)))}
    
    # generating a data frame suitable for plot generation
    del_frac <- as.vector(cn_segments[,4])
    del_frac <- append(del_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    slow_df <- data.frame(del_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d29_medslow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(slow)==0) {
    # generating a data frame suitable for plot generation
    del_frac <- c(0,0)
    coords <- append(0, size[1,2])
    slow_df <- data.frame(del_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d29_medslow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # saving the plot
  p3 <- ggplot() + geom_step(data = fast_df, aes(x = coords, y = del_frac), colour="#a1d76a", size=0.3) + geom_step(data = slow_df, aes(x = coords, y = del_frac), colour="#b2182b", size=0.3) + ylim(0,0.1) + ggtitle(paste0("Chromosome ", chr, " deletion fractions by day 29 response (slow+intermediate grouped)")) + xlab("Coordinates") + ylab("Deletion fraction") + scale_x_continuous(labels = label_comma()) + theme(panel.background=element_rect(fill="azure2", colour="azure2",size=0.5,linetype="solid"), text=element_text(size=5), plot.title=element_text(size=5))
  
  # deletions and day 79 comparison
  
  # focusing on the fast group
  fast <- subset(del, del[,8]=="negative")
  
  # conducting the next steps, if there are events reported
  if (nrow(fast)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(fast[,3])
    coords <- append(coords, fast[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the Deletion fraction affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a deletion affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      delrows <- vector()
      for (r in 1:nrow(fast)) {
        if (as.numeric(fast[r,3])==as.numeric(cn_segments[i,2])) {delrows <- append(delrows, r)}
        if (as.numeric(fast[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(fast[r,4])>cn_segments[i,2]) {delrows <- append(delrows, r)}}
      del2 <- fast[delrows,]
      cn_segments[i,4] <- (length(unique(del2[,1]))/length(fast79))}
    
    # generating a data frame suitable for plot generation
    del_frac <- as.vector(cn_segments[,4])
    del_frac <- append(del_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    fast_df <- data.frame(del_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d79_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(fast)==0) {
    # generating a data frame suitable for plot generation
    del_frac <- c(0,0)
    coords <- append(0, size[1,2])
    fast_df <- data.frame(del_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d79_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # focusing on the slow group
  slow <- subset(del, del[,8]=="positive")
  
  # conducting the next steps, if there are events reported
  if (nrow(slow)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(slow[,3])
    coords <- append(coords, slow[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the Deletion fraction affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a deletion affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      delrows <- vector()
      for (r in 1:nrow(slow)) {
        if (as.numeric(slow[r,3])==as.numeric(cn_segments[i,2])) {delrows <- append(delrows, r)}
        if (as.numeric(slow[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(slow[r,4])>cn_segments[i,2]) {delrows <- append(delrows, r)}}
      del2 <- slow[delrows,]
      cn_segments[i,4] <- (length(unique(del2[,1]))/length(slow79))}
    
    # generating a data frame suitable for plot generation
    del_frac <- as.vector(cn_segments[,4])
    del_frac <- append(del_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    slow_df <- data.frame(del_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d79_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(slow)==0) {
    # generating a data frame suitable for plot generation
    del_frac <- c(0,0)
    coords <- append(0, size[1,2])
    slow_df <- data.frame(del_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_del_fractions_d79_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # saving the plot
  p4 <- ggplot() + geom_step(data = fast_df, aes(x = coords, y = del_frac), colour="#a1d76a", size=0.3) + geom_step(data = slow_df, aes(x = coords, y = del_frac), colour="#b2182b", size=0.3) + ylim(0,0.1) + ggtitle(paste0("Chromosome ", chr, " deletion fractions by day 79 response")) + xlab("Coordinates") + ylab("Deletion fraction") + scale_x_continuous(labels = label_comma()) + theme(panel.background=element_rect(fill="azure2", colour="azure2",size=0.5,linetype="solid"), text=element_text(size=5), plot.title=element_text(size=5))
  
  # amplifications and day 15 comparison
  
  # focusing on the fast group
  fast <- subset(amp, amp[,6]=="fast")
  
  # conducting the next steps, if there are events reported
  if (nrow(fast)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(fast[,3])
    coords <- append(coords, fast[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the fraction of cases with a amplification affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a amplification affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      amprows <- vector()
      for (r in 1:nrow(fast)) {
        if (as.numeric(fast[r,3])==as.numeric(cn_segments[i,2])) {amprows <- append(amprows, r)}
        if (as.numeric(fast[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(fast[r,4])>cn_segments[i,2]) {amprows <- append(amprows, r)}}
      amp2 <- fast[amprows,]
      cn_segments[i,4] <- (length(unique(amp2[,1]))/length(fast15))}
    
    # generating a data frame suitable for plot generation
    amp_frac <- as.vector(cn_segments[,4])
    amp_frac <- append(amp_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    fast_df <- data.frame(amp_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d15_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(fast)==0) {
    # generating a data frame suitable for plot generation
    amp_frac <- c(0,0)
    coords <- append(0, size[1,2])
    fast_df <- data.frame(amp_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d15_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # focusing on the slow group
  slow <- subset(amp, amp[,6]=="slow")
  
  # conducting the next steps, if there are events reported
  if (nrow(slow)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(slow[,3])
    coords <- append(coords, slow[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the fraction of cases with a amplification affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a amplification affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      amprows <- vector()
      for (r in 1:nrow(slow)) {
        if (as.numeric(slow[r,3])==as.numeric(cn_segments[i,2])) {amprows <- append(amprows, r)}
        if (as.numeric(slow[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(slow[r,4])>cn_segments[i,2]) {amprows <- append(amprows, r)}}
      amp2 <- slow[amprows,]
      cn_segments[i,4] <- (length(unique(amp2[,1]))/length(slow15))}
    
    # generating a data frame suitable for plot generation
    amp_frac <- as.vector(cn_segments[,4])
    amp_frac <- append(amp_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    slow_df <- data.frame(amp_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d15_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(slow)==0) {
    # generating a data frame suitable for plot generation
    amp_frac <- c(0,0)
    coords <- append(0, size[1,2])
    slow_df <- data.frame(amp_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d15_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # saving the plot
  p5 <- ggplot() + geom_step(data = fast_df, aes(x = coords, y = amp_frac), colour="#a1d76a", size=0.3) + geom_step(data = slow_df, aes(x = coords, y = amp_frac), colour="#b2182b", size=0.3) + ylim(0,0.1) + ggtitle(paste0("Chromosome ", chr, " amplification fractions by day 15 response")) + xlab("Coordinates") + ylab("Amplification fraction")  + scale_x_continuous(labels = label_comma()) + theme(panel.background=element_rect(fill="seashell2", colour="seashell2",size=0.5,linetype="solid"), text=element_text(size=5), plot.title=element_text(size=5))
  
  # amplifications and day 29 comparison
  
  # focusing on the fast/intermediate group
  fast <- subset(amp, amp[,7]=="fast" | amp[,7]=="intermediate")
  
  # conducting the next steps, if there are events reported
  if (nrow(fast)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(fast[,3])
    coords <- append(coords, fast[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the fraction of cases with a amplification affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a amplification affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      amprows <- vector()
      for (r in 1:nrow(fast)) {
        if (as.numeric(fast[r,3])==as.numeric(cn_segments[i,2])) {amprows <- append(amprows, r)}
        if (as.numeric(fast[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(fast[r,4])>cn_segments[i,2]) {amprows <- append(amprows, r)}}
      amp2 <- fast[amprows,]
      cn_segments[i,4] <- (length(unique(amp2[,1]))/(length(fast29)+length(med29)))}
    
    # generating a data frame suitable for plot generation
    amp_frac <- as.vector(cn_segments[,4])
    amp_frac <- append(amp_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    fast_df <- data.frame(amp_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d29_fastmed.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(fast)==0) {
    # generating a data frame suitable for plot generation
    amp_frac <- c(0,0)
    coords <- append(0, size[1,2])
    fast_df <- data.frame(amp_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d29_fastmed.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # focusing on the slow group
  slow <- subset(amp, amp[,7]=="slow")
  
  # conducting the next steps, if there are events reported
  if (nrow(slow)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(slow[,3])
    coords <- append(coords, slow[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the fraction of cases with a amplification affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a amplification affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      amprows <- vector()
      for (r in 1:nrow(slow)) {
        if (as.numeric(slow[r,3])==as.numeric(cn_segments[i,2])) {amprows <- append(amprows, r)}
        if (as.numeric(slow[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(slow[r,4])>cn_segments[i,2]) {amprows <- append(amprows, r)}}
      amp2 <- slow[amprows,]
      cn_segments[i,4] <- (length(unique(amp2[,1]))/length(slow29))}
    
    # generating a data frame suitable for plot generation
    amp_frac <- as.vector(cn_segments[,4])
    amp_frac <- append(amp_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    slow_df <- data.frame(amp_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d29_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(slow)==0) {
    # generating a data frame suitable for plot generation
    amp_frac <- c(0,0)
    coords <- append(0, size[1,2])
    slow_df <- data.frame(amp_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d29_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # saving the plot
  p6 <- ggplot() + geom_step(data = fast_df, aes(x = coords, y = amp_frac), colour="#a1d76a", size=0.3) + geom_step(data = slow_df, aes(x = coords, y = amp_frac), colour="#b2182b", size=0.3) + ylim(0,0.1) + ggtitle(paste0("Chromosome ", chr, " amplification fractions by day 29 response (fast+intermediate grouped)")) + xlab("Coordinates") + ylab("Amplification fraction")  + scale_x_continuous(labels = label_comma()) + theme(panel.background=element_rect(fill="seashell2", colour="seashell2",size=0.5,linetype="solid"), text=element_text(size=5), plot.title=element_text(size=5))
  
  # then, focusing on the fast group only
  fast <- subset(amp, amp[,7]=="fast")
  
  # conducting the next steps, if there are events reported
  if (nrow(fast)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(fast[,3])
    coords <- append(coords, fast[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the fraction of cases with a amplification affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a amplification affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      amprows <- vector()
      for (r in 1:nrow(fast)) {
        if (as.numeric(fast[r,3])==as.numeric(cn_segments[i,2])) {amprows <- append(amprows, r)}
        if (as.numeric(fast[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(fast[r,4])>cn_segments[i,2]) {amprows <- append(amprows, r)}}
      amp2 <- fast[amprows,]
      cn_segments[i,4] <- (length(unique(amp2[,1]))/(length(fast29)))}
    
    # generating a data frame suitable for plot generation
    amp_frac <- as.vector(cn_segments[,4])
    amp_frac <- append(amp_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    fast_df <- data.frame(amp_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d29_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(fast)==0) {
    # generating a data frame suitable for plot generation
    amp_frac <- c(0,0)
    coords <- append(0, size[1,2])
    fast_df <- data.frame(amp_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d29_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # continuing with the intermediate/slow group
  slow <- subset(amp, amp[,7]=="slow" | amp[,7]=="intermediate")
  
  # conducting the next steps, if there are events reported
  if (nrow(slow)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(slow[,3])
    coords <- append(coords, slow[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the fraction of cases with a amplification affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a amplification affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      amprows <- vector()
      for (r in 1:nrow(slow)) {
        if (as.numeric(slow[r,3])==as.numeric(cn_segments[i,2])) {amprows <- append(amprows, r)}
        if (as.numeric(slow[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(slow[r,4])>cn_segments[i,2]) {amprows <- append(amprows, r)}}
      amp2 <- slow[amprows,]
      cn_segments[i,4] <- (length(unique(amp2[,1]))/(length(slow29)+length(med29)))}
    
    # generating a data frame suitable for plot generation
    amp_frac <- as.vector(cn_segments[,4])
    amp_frac <- append(amp_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    slow_df <- data.frame(amp_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d29_medslow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(slow)==0) {
    # generating a data frame suitable for plot generation
    amp_frac <- c(0,0)
    coords <- append(0, size[1,2])
    slow_df <- data.frame(amp_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d29_medslow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # saving the plot
  p7 <- ggplot() + geom_step(data = fast_df, aes(x = coords, y = amp_frac), colour="#a1d76a", size=0.3) + geom_step(data = slow_df, aes(x = coords, y = amp_frac), colour="#b2182b", size=0.3) + ylim(0,0.1) + ggtitle(paste0("Chromosome ", chr, " amplification fractions by day 29 response (slow+intermediate grouped)")) + xlab("Coordinates") + ylab("Amplification fraction")  + scale_x_continuous(labels = label_comma()) + theme(panel.background=element_rect(fill="seashell2", colour="seashell2",size=0.5,linetype="solid"), text=element_text(size=5), plot.title=element_text(size=5))
  
  # amplifications and day 79 comparison
  
  # focusing on the fast group
  fast <- subset(amp, amp[,8]=="negative")
  
  # conducting the next steps, if there are events reported
  if (nrow(fast)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(fast[,3])
    coords <- append(coords, fast[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the fraction of cases with a amplification affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a amplification affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      amprows <- vector()
      for (r in 1:nrow(fast)) {
        if (as.numeric(fast[r,3])==as.numeric(cn_segments[i,2])) {amprows <- append(amprows, r)}
        if (as.numeric(fast[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(fast[r,4])>cn_segments[i,2]) {amprows <- append(amprows, r)}}
      amp2 <- fast[amprows,]
      cn_segments[i,4] <- (length(unique(amp2[,1]))/length(fast79))}
    
    # generating a data frame suitable for plot generation
    amp_frac <- as.vector(cn_segments[,4])
    amp_frac <- append(amp_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    fast_df <- data.frame(amp_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d79_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(fast)==0) {
    # generating a data frame suitable for plot generation
    amp_frac <- c(0,0)
    coords <- append(0, size[1,2])
    fast_df <- data.frame(amp_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d79_fast.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # focusing on the slow group
  slow <- subset(amp, amp[,8]=="positive")
  
  # conducting the next steps, if there are events reported
  if (nrow(slow)!=0) {
    # determining each unique breakpoint, and sorting them
    coords <- as.vector(slow[,3])
    coords <- append(coords, slow[,4])
    coords <- unique(coords)
    coords <- sort(coords)
    
    # generating a new table for the copy number segments based on the unique breakpoints
    # which will describe the fraction of cases with a amplification affecting that segment
    cn_segments <- matrix(ncol=4, nrow=1)
    cn_segments[1,] <- c(chr,0,min(coords),0)
    newrows <- matrix(ncol=4, nrow=length(coords))
    newrows[,1] <- chr
    for (c in 2:length(coords)) {
      newrows[c,2] <- coords[c-1]
      newrows[c,3] <- coords[c]}
    newrows <- newrows[-1,]
    cn_segments <- rbind(cn_segments, newrows)
    newrow <- matrix(ncol=4, nrow=1)
    newrow[1,] <- c(chr,max(coords),size[1,2],0)
    cn_segments <- rbind(cn_segments, newrow)
    
    # determining the cases with a amplification affecting each segment
    # and adding the fraction of affected cases to the table
    for (i in 2:length(coords)) {
      amprows <- vector()
      for (r in 1:nrow(slow)) {
        if (as.numeric(slow[r,3])==as.numeric(cn_segments[i,2])) {amprows <- append(amprows, r)}
        if (as.numeric(slow[r,3])<as.numeric(cn_segments[i,2]) & as.numeric(slow[r,4])>cn_segments[i,2]) {amprows <- append(amprows, r)}}
      amp2 <- slow[amprows,]
      cn_segments[i,4] <- (length(unique(amp2[,1]))/length(slow79))}
    
    # generating a data frame suitable for plot generation
    amp_frac <- as.vector(cn_segments[,4])
    amp_frac <- append(amp_frac, 0)
    coords <- as.vector(cn_segments[,2])
    coords <- append(coords, size[1,2])
    slow_df <- data.frame(amp_frac, coords)
    
    # writing the cn_segments file out
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d79_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # conducting the next steps, if there are no events reported
  if (nrow(slow)==0) {
    # generating a data frame suitable for plot generation
    amp_frac <- c(0,0)
    coords <- append(0, size[1,2])
    slow_df <- data.frame(amp_frac, coords)
    # also generating a cn_segments file and writing it out
    cn_segments <- matrix(ncol=4, nrow=1, c(chr, 0, size[1,2], 0))
    colnames(cn_segments) <- c("chr", "start_seg", "end_seg", "amp_frac")
    write.table(cn_segments, paste0("cn_segment_files/chr", chr, "_amp_fractions_d79_slow.bed"), quote=F, sep="\t", col.names = T, row.names = F)}
  
  # saving the plot
  p8 <- ggplot() + geom_step(data = fast_df, aes(x = coords, y = amp_frac), colour="#a1d76a", size=0.3) + geom_step(data = slow_df, aes(x = coords, y = amp_frac), colour="#b2182b", size=0.3) + ylim(0,0.1) + ggtitle(paste0("Chromosome ", chr, " amplification fractions by day 79 response")) + xlab("Coordinates") + ylab("Amplification fraction")  + scale_x_continuous(labels = label_comma()) + theme(panel.background=element_rect(fill="seashell2", colour="seashell2",size=0.5,linetype="solid"), text=element_text(size=5), plot.title=element_text(size=5))
  
  # printing the plots
  pdf(paste0("chr", chr, "_collective_step_charts.pdf"))
  grid.arrange(p1,p5,p2,p6,p3,p7,p4,p8,ncol=2)
  dev.off()}

## oncoprint generation
library(ComplexHeatmap)

# response table
resp <- read.table("sample_responses_cnv_cohort.txt", header=T)

# WGS CNV results
nc <- read.table("/Volumes/groups/allseq/students/sanni_wrk/ETV6-RUNX1_cases/balsamic_results/data_files/lncRNA_genes/CNVs_affecting_nc_genes.tab")
cnv <- read.table("/Volumes/groups/allseq/students/sanni_wrk/ETV6-RUNX1_cases/balsamic_results/data_files/GEPARD_CNV/annotations/collective_coding_CNVs.tab", header=T)
# combining the tables
gene_type <- matrix(ncol=1, nrow=nrow(cnv))
gene_type[,1] <- c("protein_coding")
cnv <- cbind(cnv[,c(1,2,5,7)], gene_type, cnv[,8])
nc <- nc[,1:6]
colnames(cnv) <- colnames(nc)
getab <- rbind(cnv, nc)
getab <- getab[,c(1:4,6)]
# panelseq and methylation array results
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

# determining the genes of interest
# EOI genes affecting chr12
del <- read.table("d29_del_fishers_test_results_combined_de_info.tab", header=T)
del <- subset(del, del[,2]=="12")
amp <- read.table("d29_amp_fishers_test_results_combined_de_info.tab", header=T)
amp <- subset(amp, amp[,2]=="12")
genes <- append(del[,1], amp[,1])
genes <- unique(genes)

# reading in the file with the gene coordinates
coords <- read.table("collective_CNV_filt_annot_per_gene.tab", header=T)
coords <- unique(coords[,c(1,9,10,7)])
# ordering the table
coords <- coords[order(coords$Tx_start),]
coords <- subset(coords, coords[,1]=="12")
coords <- subset(coords, coords[,4] %in% genes)
genes <- unique(coords[,4])

# determining the EOI sample order
slow <- subset(resp, resp[,3]=="slow")
med <- subset(resp, resp[,3]=="intermediate")
fast <- subset(resp, resp[,3]=="fast")
resp29 <- rbind(slow, med, fast)
samples <- resp29[,1]

# generating a table for oncoprint generation
toplot <- matrix(ncol=length(samples), nrow=length(genes))
colnames(toplot) <- samples
rownames(toplot) <- genes
for (g in genes) {
  for (s in samples) {
    tab2 <- subset(tab, tab[,1]==s & tab[,4]==g)
    if (length(unique(tab2[,3]))==1) {
      toplot[g,s] <- tab2[1,3]}
    if (length(unique(tab2[,3]))==2) {
      toplot[g,s] <- paste0(unique(tab2[,3])[1], ",", unique(tab2[,3])[2])}}}

# correcting the names
for (g in genes) {
  for (s in samples) {
    toplot[g,s] <- gsub("DEL", "Deletion", toplot[g,s])
    toplot[g,s] <- gsub("AMP", "Gain", toplot[g,s] )}}
for (i in 1:nrow(resp29)) {
  for (j in 1:ncol(resp29)) {
    resp29[i,j] <- gsub("slow", "Slow", resp29[i,j])
    resp29[i,j] <- gsub("intermediate", "Intermediate", resp29[i,j])
    resp29[i,j] <- gsub("fast", "Fast", resp29[i,j])
    resp29[i,j] <- gsub("positive", "Positive", resp29[i,j])
    resp29[i,j] <- gsub("negative", "Negative", resp29[i,j])
    resp29[i,j] <- gsub("unknown", "Unknown", resp29[i,j])}}

# defining the alter_fun for generating the CNV/noncoding/intergenic SV oncoprint
alter_fun = list(background = function(x, y, w, h) {grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "grey95", col = NA))},
                 Deletion = function(x, y, w, h) {grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill="#5ab4ac", col = NA))},
                 Gain = function(x, y, w, h) {grid.rect(x, y, w-unit(0.5, "mm"), h*0.6, gp = gpar(fill="mediumorchid3", col = NA))})

# determining the colours for the top annotation
col_fun1 <- c(Fast="#a1d76a", Intermediate="#f1a340", Slow="#b2182b", Unknown="grey")
col_fun2 <- c(Positive="black", Unknown="grey", Negative="#f7f7f7")

# drawing the oncoprint
pdf("manuscript_CNV_oncoprint_EOI_chr12.pdf", width=16, height=16)
oncoPrint(toplot, column_order=samples, row_order=genes, alter_fun=alter_fun, show_pct=FALSE, show_column_names=FALSE, row_names_gp=gpar(fontsize=8), top_annotation=HeatmapAnnotation("Mid-induction response"=resp29[,2], "EOC MRD status"=resp29[,4], "EOI response"=resp29[,3], col=list("Mid-induction response"=col_fun1, "EOI response"=col_fun1, "EOC MRD status"=col_fun2)))
dev.off()
