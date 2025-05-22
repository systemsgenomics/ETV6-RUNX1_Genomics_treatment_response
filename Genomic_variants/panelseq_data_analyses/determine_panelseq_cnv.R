# determine the CNV deletions and amplifications from panelseq CNVkit output
library(dplyr)

# panelseq sample responses
panel <- read.table("MRD_panel_seq.txt", header=T)
panel <- panel[,c(1,6,5,7)]

# generating panelseq CNV files
for (s in panel[,1]) {
  
  # CNVkit output
  tab <- read.table(paste0("case", s, "/tumor.merged.cns"), sep="\t", header=T)
  
  # remove X and Y
  tab <- subset(tab, tab$chromosome %in% c(1:22))
  
  # determine deletions and amplifications
  # segments with log2 >= 0.3 determined to be amplifications (gains)
  # segments with log2 <= -0.4 determined to be deletions
  # otherwise neutral
  newcols <- matrix(ncol=2, nrow=nrow(tab))
  for (i in 1:nrow(tab)) {
    if (as.numeric(tab[i,5]) >= 0.3) {newcols[i,1] <- "AMP"}
    if (as.numeric(tab[i,5]) <= -0.4) {newcols[i,1] <- "DEL"}}
  newcols[,2] <- s
  tab <- cbind(tab[,1:3], newcols)
  
  # merging nearby events if they are less than 10000 bp apart
  remove <- vector()
  for (i in 1:nrow(tab)) {
    if (isTRUE(tab[i,4]=="DEL" & tab[i+1,4]=="DEL" & abs(as.numeric(tab[i,3])-as.numeric(tab[i+1,2])) < 10000)) {
      remove <- append(remove, i)
      remove <- append(remove, i+1)}}
  remove <- unique(remove)
  mergeddel <- tab[remove, ]
  if (nrow(mergeddel)!=0) {
    tab <- tab[-remove,]
    for (c in unique(mergeddel[,1])) {
      chr <- subset(mergeddel, mergeddel[,1]==c)
      bnd <- vector()
      for (i in 1:nrow(chr)) {
        if (isTRUE(abs(as.numeric(chr[i,3])-as.numeric(chr[i+1,2])) >= 10000)) {bnd <- append(bnd, i)}}
      if (length(bnd)==0) {
        newrow <- matrix(ncol=5, nrow=1)
        newrow[1,1] <- chr[1,1]
        newrow[1,2] <- min(as.numeric(chr[,2]))
        newrow[1,3] <- max(as.numeric(chr[,3]))
        newrow[1,4] <- chr[1,4]
        newrow[1,5] <- s
        colnames(newrow) <- colnames(tab)
        tab <- rbind(tab, newrow)}
      if (length(bnd)>0) {
        newrow <- matrix(ncol=5, nrow=(length(bnd)+1))
        newrow[,1] <- chr[1,1] 
        newrow[,4] <- chr[1,4]
        newrow[,5] <- s
        newrow[1,2] <- min(as.numeric(chr[1:bnd[1],2]))
        newrow[1,3] <- max(as.numeric(chr[1:bnd[1],3]))
        for (l in 1:length(bnd)) {
          if (l==length(bnd)) {
            newrow[l+1,2] <- min(as.numeric(chr[(bnd[l]+1):nrow(chr),2]))
            newrow[l+1,3] <- max(as.numeric(chr[(bnd[l]+1):nrow(chr),3]))}
          if (l != length(bnd)) {
            newrow[l+1,2] <- min(as.numeric(chr[(bnd[l]+1):(bnd[l+1]),2]))
            newrow[l+1,3] <- max(as.numeric(chr[(bnd[l]+1):(bnd[l+1]),3]))}}
        colnames(newrow) <- colnames(tab)
        tab <- rbind(tab, newrow)}}}
  remove <- vector()
  for (i in 1:nrow(tab)) {
    if (isTRUE(tab[i,4]=="AMP" & tab[i+1,4]=="AMP" & abs(as.numeric(tab[i,3])-as.numeric(tab[i+1,2])) < 10000)) {
      remove <- append(remove, i)
      remove <- append(remove, i+1)}}
  remove <- unique(remove)
  mergedamp <- tab[remove, ]
  if (nrow(mergedamp)!=0) {tab <- tab[-remove,]
  for (c in unique(mergedamp[,1])) {
    chr <- subset(mergedamp, mergedamp[,1]==c)
    bnd <- vector()
    for (i in 1:nrow(chr)) {
      if (isTRUE(abs(as.numeric(chr[i,3])-as.numeric(chr[i+1,2])) >= 10000)) {bnd <- append(bnd, i)}}
    if (length(bnd)==0) {
      newrow <- matrix(ncol=5, nrow=1)
      newrow[1,1] <- chr[1,1]
      newrow[1,2] <- min(as.numeric(chr[,2]))
      newrow[1,3] <- max(as.numeric(chr[,3]))
      newrow[1,4] <- chr[1,4]
      newrow[1,5] <- s
      colnames(newrow) <- colnames(tab)
      tab <- rbind(tab, newrow)}
    if (length(bnd)>0) {
      newrow <- matrix(ncol=5, nrow=(length(bnd)+1))
      newrow[,1] <- chr[1,1] 
      newrow[,4] <- chr[1,4]
      newrow[,5] <- s
      newrow[1,2] <- min(as.numeric(chr[1:bnd[1],2]))
      newrow[1,3] <- max(as.numeric(chr[1:bnd[1],3]))
      for (l in 1:length(bnd)) {
        if (l==length(bnd)) {
          newrow[l+1,2] <- min(as.numeric(chr[(bnd[l]+1):nrow(chr),2]))
          newrow[l+1,3] <- max(as.numeric(chr[(bnd[l]+1):nrow(chr),3]))}
        if (l != length(bnd)) {
          newrow[l+1,2] <- min(as.numeric(chr[(bnd[l]+1):(bnd[l+1]),2]))
          newrow[l+1,3] <- max(as.numeric(chr[(bnd[l]+1):(bnd[l+1]),3]))}}
      colnames(newrow) <- colnames(tab)
      tab <- rbind(tab, newrow)}}}
  
  # exclude neutral regions
  tab <- subset(tab, !is.na(tab[,4]))
  
  # remove events smaller than 1 Mb in size
  remove <- vector()
  for (i in 1:nrow(tab)) {
    if (diff(c(as.numeric(tab[i,2]), as.numeric(tab[i,3]))) < 1000000) {remove <- append(remove, i)}}
  if (length(remove)!=0) {tab <- tab[-remove,]}
  
  # writing the table out
  write.table(tab, paste0("bed_files/", s, ".CNVs.bed"), quote=F, sep="\t", col.names=F, row.names=F)
}
