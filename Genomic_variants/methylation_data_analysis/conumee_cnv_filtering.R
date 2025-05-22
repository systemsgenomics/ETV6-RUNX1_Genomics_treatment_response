# signal- and size based filtering of methylation data copy number calls (by conumee)

# reading in the file containing the case IDs and d29 response info
mrd <- read.table("d29_response.txt", header=T)

# saving the sample ids as vectors, modifying them to be in the same format as in the conumee file names
nopho <- subset(mrd, mrd[,3]=="NOPHO")
fmca <- subset(mrd, mrd[,3]=="FMCA")
nopho <- as.vector(nopho[,2])
fmca <- as.vector(fmca[,2])
for (l in 1:length(nopho)) {
  nopho[l] <- gsub("ALL_", "nopho.", nopho[l])}
for (l in 1:length(fmca)) {
  fmca[l] <- gsub("ALL_", "all.", fmca[l])}
# removing case 690 since it has no diagnosis data
fmca <- fmca[-13]

# reading in the segment files and filtering based on signal and size
for (s in nopho) {
  tab <- read.table(paste0(s, ".0.seg"), header=T)
  
  # accounting segments with median signal < -0.05 as deletions
  del <- subset(tab, tab$seg.median < -0.05)
  
  # accounting segments with median signal > 0.1 as duplications
  dup <- subset(tab, tab$seg.median > 0.1)
  
  # adding columns to each describing the variant type
  variant_type <- matrix(ncol=1, nrow=nrow(del), c("DEL"))
  del <- cbind(del, variant_type)
  variant_type <- matrix(ncol=1, nrow=nrow(dup), c("AMP"))
  dup <- cbind(dup, variant_type)
  
  # combining both tables
  tab <- rbind(del, dup)
  
  # adding a column describing the event size
  size <- matrix(ncol=1, nrow=nrow(tab))
  for (i in 1:nrow(tab)) {
    size[i,1] <- abs(tab$loc.start[i]-tab$loc.end[i])}
  tab <- cbind(tab[,c(1:4)], size, tab[,5:10])
  
  # removing the chr prefix
  for (i in 1:nrow(tab)) {
    vec <- unlist(strsplit(tab[i,2], "r"))
    tab[i,2] <- vec[2]}
  
  # separating CNVs based on size (1 Mb cutoff)
  cnv <- subset(tab, tab$size >= 1000000)

  # writing out as bed files
  s2 <- unlist(strsplit(s, "[.]"))
  s2 <- paste0("ALL_", s2[2])
  if (nrow(cnv)!=0) {
    cnv <- cnv[,c(2:5,11,1,9:10)]
    cnv$ID <- s2}
  write.table(cnv, paste0(s2, ".conumee_CNVs.bed"), col.names=F, row.names=F, quote=F, sep="\t")
}

for (s in fmca) {
  tab <- read.table(paste0(s, ".diag.seg"), header=T)
  # accounting segments with median signal < -0.05 as deletions
  del <- subset(tab, tab$seg.median < -0.05)
  
  # accounting segments with median signal > 0.1 as duplications
  dup <- subset(tab, tab$seg.median > 0.1)
  
  # adding columns to each describing the variant type
  variant_type <- matrix(ncol=1, nrow=nrow(del), c("DEL"))
  del <- cbind(del, variant_type)
  variant_type <- matrix(ncol=1, nrow=nrow(dup), c("AMP"))
  dup <- cbind(dup, variant_type)
  
  # combining both tables
  tab <- rbind(del, dup)
  
  # adding a column describing the event size
  size <- matrix(ncol=1, nrow=nrow(tab))
  for (i in 1:nrow(tab)) {
    size[i,1] <- abs(tab$loc.start[i]-tab$loc.end[i])}
  tab <- cbind(tab[,c(1:4)], size, tab[,5:10])
  
  # removing the chr prefix
  for (i in 1:nrow(tab)) {
    vec <- unlist(strsplit(tab[i,2], "r"))
    tab[i,2] <- vec[2]}
  
  # separating CNVs based on size (1 Mb cutoff)
  cnv <- subset(tab, tab$size >= 1000000)
  
  # writing out as bed files
  s2 <- unlist(strsplit(s, "[.]"))
  s2 <- paste0("ALL_", s2[2])
  if (nrow(cnv)!=0) {
    cnv <- cnv[,c(2:5,11,1,9:10)]
    cnv$ID <- s2}
  write.table(cnv, paste0(s2, ".conumee_CNVs.bed"), col.names=F, row.names=F, quote=F, sep="\t")
}
