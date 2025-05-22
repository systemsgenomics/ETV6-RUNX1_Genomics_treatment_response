# process filtered and annotated array CNVs

# case IDs
samples2 <- c("ALL_249", "ALL_269", "ALL_251", "ALL_284", "ALL_301", "ALL_330", "ALL_309", "ALL_305", "ALL_354", "ALL_326", "ALL_337", "ALL_322", "ALL_444", "ALL_364", "ALL_448", "ALL_440", "ALL_450", "ALL_381", "ALL_434", "ALL_446", "ALL_439", "ALL_430", "ALL_384", "ALL_386", "ALL_394", "ALL_400", "ALL_525", "ALL_458", "ALL_504", "ALL_510", "ALL_511", "ALL_512", "ALL_488", "ALL_493", "ALL_500", "ALL_540", "ALL_589", "ALL_581", "ALL_574", "ALL_586", "ALL_582", "ALL_548", "ALL_550", "ALL_551", "ALL_731", "ALL_732", "ALL_716", "ALL_701", "ALL_706", "ALL_906", "ALL_936", "ALL_949", "ALL_1027", "ALL_1051", "ALL_1060", "ALL_1069", "ALL_1071", "ALL_1085", "ALL_1098", "ALL_1123", "ALL_1158", "ALL_398", "ALL_755", "ALL_758", "ALL_695", "ALL_697", "ALL_698", "ALL_705", "ALL_709", "ALL_715", "ALL_881", "ALL_871", "ALL_883", "ALL_890", "ALL_893", "ALL_895", "ALL_922", "ALL_928", "ALL_933", "ALL_937", "ALL_939", "ALL_943", "ALL_963", "ALL_917", "ALL_958", "ALL_886", "ALL_887", "ALL_905", "ALL_916", "ALL_971", "ALL_980", "ALL_989", "ALL_759", "ALL_1014", "ALL_760", "ALL_1020", "ALL_1021", "ALL_1042", "ALL_995", "ALL_998", "ALL_1023", "ALL_1025", "ALL_1032", "ALL_1052", "ALL_981", "ALL_983", "ALL_1012", "ALL_1024", "ALL_970", "ALL_977", "ALL_997", "ALL_1115", "ALL_1104", "ALL_1119", "ALL_1120", "ALL_1112", "ALL_1126", "ALL_1134", "ALL_1133", "ALL_1144", "ALL_1172", "ALL_1182", "ALL_1180", "ALL_1181", "ALL_1186", "ALL_1184", "ALL_1191", "ALL_1192", "ALL_1143", "ALL_1146", "ALL_1147", "ALL_1163", "ALL_1168")

# conducting the following steps for the data of each sample
for (s in samples) {
  
  # CNV resultd
  tab <- read.csv(paste0(s, ".CNV.annotated.tsv"), sep="\t")
  # simplifying the table
  tab <- tab[,c(2:4,6:8,12,14:18,21,24,25,22)]
  # separating the full and split rows
  full <- subset(tab, tab[,7]=="full")
  split <- subset(tab, tab[,7]=="split")
  
  # simplifying the full table
  full <- full[,c(1:5,8,9,6)]
  # defining column names and writing the table out
  colnames(full) <- c("CNV_chrom", "CNV_start", "CNV_end", "CNV_length", "CNV_type", "Gene_name", "Gene_count", "sampleID")
  write.table(full, paste0(s, "_CNV_filt_annot_per_event.tab"), sep="\t", quote=F, col.names = T, row.names = F)
  
  # adding the consequence annotations to the gene-wise split table
  split <- split[,c(1:5,7:16,6)]
  consequence <- matrix(ncol=1, nrow=nrow(split))
  for (i in 1:nrow(split)) {
    # first annotating deletions which overlap 100% with the coding sequence with a loss consequence
    if (isTRUE(split[i,5]=="DEL" & split[i,12]==100)) {consequence[i,1] <- c("loss")}
    # annotating deletions which overlap 1-99% with the coding sequence with a coding_sequence_variant consequence
    if (isTRUE(split[i,5]=="DEL" & split[i,12]>0 & split[i,12]<100)) {consequence[i,1] <- c("coding_sequence_variant")}
    # annotating AMPlications which overlap 100% with the coding sequence with a gain consequence
    if (isTRUE(split[i,5]=="AMP" & split[i,12]==100)) {consequence[i,1] <- c("gain")}
    # annotating AMPlications which overlap 1-99% with the coding sequence with a coding_sequence_variant consequence
    if (isTRUE(split[i,5]=="AMP" & split[i,12]>0 & split[i,12]<100)) {consequence[i,1] <- c("coding_sequence_variant")}
    # annotating deletions which overlap 100% with the transcript of a non-coding gene with a nc_gene_loss consequence
    if (isTRUE(split[i,5]=="DEL" & split[i,12]==0 & split[i,13]=="txStart-txEnd")) {consequence[i,1] <- c("nc_gene_loss")}
    # annotating deletions which overlap 1-99% with the sequence of a non-coding gene with a nc_gene_sequence_variant consequence
    if (isTRUE(split[i,5]=="DEL" & as.numeric(split[i,2])>as.numeric(split[i,10]) & as.numeric(split[i,3])>as.numeric(split[i,11]) & as.numeric(split[i,2])<as.numeric(split[i,11]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR" | split[i,5]=="DEL" & as.numeric(split[i,2])<as.numeric(split[i,10]) & as.numeric(split[i,3])<as.numeric(split[i,11]) & as.numeric(split[i,3])>as.numeric(split[i,10]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR" | split[i,5]=="DEL" & as.numeric(split[i,2])>as.numeric(split[i,10]) & as.numeric(split[i,3])<as.numeric(split[i,11]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR")) {consequence[i,1] <- c("nc_gene_sequence_variant")}
    # annotating AMPlications which overlap 100% with the transcript of a non-coding gene with a nc_gene_gain consequence
    if (isTRUE(split[i,5]=="AMP" & split[i,12]==0 & split[i,13]=="txStart-txEnd")) {consequence[i,1] <- c("nc_gene_gain")}
    # annotating AMPlications which overlap 1-99% with the sequence of a non-coding gene with a nc_gene_sequence_variant consequence
    if (isTRUE(split[i,5]=="AMP" & as.numeric(split[i,2])>as.numeric(split[i,10]) & as.numeric(split[i,3])>as.numeric(split[i,11]) & as.numeric(split[i,2])<as.numeric(split[i,11]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR" | split[i,5]=="AMP" & as.numeric(split[i,2])<as.numeric(split[i,10]) & as.numeric(split[i,3])<as.numeric(split[i,11]) & as.numeric(split[i,3])>as.numeric(split[i,10]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR" | split[i,5]=="AMP" & as.numeric(split[i,2])>as.numeric(split[i,10]) & as.numeric(split[i,3])<as.numeric(split[i,11]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR")) {consequence[i,1] <- c("nc_gene_sequence_variant")}
    # annotating deletions and AMPlications affecting an intronic/UTR regions accordingly
    if (isTRUE(split[i,5]=="DEL" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="CDS")) {consequence[i,1] <- c("intronic")}
    if (isTRUE(split[i,5]=="DEL" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="3'UTR")) {consequence[i,1] <- c("3_UTR_variant")}
    if (isTRUE(split[i,5]=="DEL" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="5'UTR")) {consequence[i,1] <- c("5_UTR_variant")}
    if (isTRUE(split[i,5]=="AMP" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="CDS" & split[i,15]=="no")) {consequence[i,1] <- c("intronic")}
    if (isTRUE(split[i,5]=="AMP" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="CDS" & split[i,15]=="yes")) {consequence[i,1] <- c("coding_sequence_variant")}
    if (isTRUE(split[i,5]=="AMP" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="3'UTR")) {consequence[i,1] <- c("3_UTR_variant")}
    if (isTRUE(split[i,5]=="AMP" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="5'UTR")) {consequence[i,1] <- c("5_UTR_variant")}}
  split <- cbind(split[,c(1:5)], consequence, split[,c(7,9:12,15,13,14,16)])
  # defining column names and writing the table out
  colnames(split) <- c("CNV_chrom", "CNV_start", "CNV_end", "CNV_length", "CNV_type", "Consequence", "Gene_name", "Tx", "Tx_start", "Tx_end", "Overlapped_CDS_percent", "Frameshift","Location", "Location2", "sampleID")
  write.table(split, paste0(s, "_CNV_filt_annot_per_gene.tab"), sep="\t", quote=F, col.names = T, row.names = F)}

# generating collective tables for per event and per gene tables
event <- matrix(ncol=8, nrow=1)
gene <- matrix(ncol=15, nrow=1)
for (s in samples) {
  full <- read.table(paste0(s, "_CNV_filt_annot_per_event.tab"), sep="\t", header=T,quote = "")
  split <- read.table(paste0(s, "_CNV_filt_annot_per_gene.tab"), sep="\t", header=T,quote = "")
  colnames(event) <- colnames(full)
  colnames(gene) <- colnames(split)
  event <- rbind(event, full)
  gene <- rbind(gene, split)}
event <- event[-1,]
gene <- gene[-1,]
# writing the tables out
write.table(event, "collective_CNV_filt_annot_per_event.tab", sep="\t", quote=F, col.names=T, row.names = F)
write.table(gene, "collective_CNV_filt_annot_per_gene.tab", sep="\t", quote=F, col.names=T, row.names = F)
