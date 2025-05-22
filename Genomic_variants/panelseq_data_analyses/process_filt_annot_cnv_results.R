# process filtered and annotated panelseq CNVs

# case IDs
samples <- c("08-307", "10-047", "12-098", "10-095", "08-222", "08-186", "10219", "11-244", "13-185", "08-204", "10230", "11-245", "13-187", "08-209", "10264", "11-246", "13-191", "08-212", "10267", "11-255", "13-196", "08-219", "10274", "11-257", "13-197", "08-271", "10279", "11-263", "13-214", "08-286", "10289", "11-290", "13-223", "08-297", "10323", "11-299", "13-226", "08-304", "10339", "11-313", "13-236", "09-004", "10345", "11-317", "13-246", "09-019", "10346", "11-328", "13-270", "09-024", "10348", "11013", "13-271", "09-036", "10349", "11016", "13-285", "09-043", "10355", "11019", "14-023", "09-051", "10366", "11025", "14-024", "09-088", "10369", "11042", "14-044", "09-119", "10372", "11253", "14-046", "09-120", "10376", "12-014", "14-055", "09-121", "10377", "12-017", "14-070", "09-143", "10379", "12-033", "14-074", "09-162", "10399", "12-038", "14-087", "09-175", "10400", "12-040", "14-091", "09-185", "10416", "12-063", "14-103", "09-191", "10421", "12-126", "14-128", "09-215", "10432", "12-128", "14-144", "09-296", "10446", "12-130",  "14-161", "09-308", "10452", "12-153", "14-202", "10-024", "10460", "12-154", "14-211", "10-025", "10461", "12-168", "14-220", "10-055", "10469", "12-179", "14-223", "10-059", "10470", "12-180", "14-239", "10-064", "10480", "12-200", "14-250", "10-078", "10520", "12-213", "14-251", "10-113", "10532", "12-215", "14-254", "10-126", "10534", "12-216", "14-258", "10-127", "10541",  "12-222", "14-283", "10-175", "10611", "12-225", "14-284", "10-188", "10634", "12-227", "15-001", "10-191", "10660", "12-234", "15-010", "10-192", "10666", "12-236", "15-015", "10-198", "10685", "12-245", "15-019", "10-210", "10704", "12-256", "15-048", "10-217", "10746", "12-258", "15-106", "10-223", "10747", "12-260", "15-123", "10-232", "10759", "12-261", "15-139", "10-234", "10806", "12-264", "15-145", "10-235", "10812", "13-007", "15-148", "10-263", "10817", "13-011", "15-149", "10-291", "10855", "13-016", "15-150", "10-308", "10882", "13-018", "15-155", "10-331", "10888", "13-050", "15-176", "10-333", "10932", "13-061", "15-182", "10034", "10933", "13-069", "15-196", "10035", "10954", "13-078", "15-208", "10037", "10968", "13-092", "15-223", "10044", "10969", "13-101", "15-224", "10060", "11-065", "13-111", "15-228", "10082", "11-067", "13-112", "15-254", "10092", "11-098", "13-115", "15-264", "10104", "11-165", "13-128", "16-004", "10127", "11-171", "13-131", "10178", "11-189", "13-137", "10202",  "11-212", "13-155", "10213", "11-236", "13-171", "10216", "11-237", "13-184")

# conducting the following steps for the data of each sample
for (s in samples) {
  
  # annotated results
  tab <- read.csv(paste0(s, ".CNV.annotated.tsv"), sep="\t")
  # simplifying the table
  tab <- tab[,c(2:4,6,7,9,11:15,18,21,22,19)]
  # separating the full and split rows
  full <- subset(tab, tab[,6]=="full")
  split <- subset(tab, tab[,6]=="split")
  
  # adding the event size to the full table
  size <- matrix(ncol=1, nrow=nrow(full))
  for (i in 1:nrow(full)) {
    size[i,1] <- diff(c(as.numeric(full[i,2]), as.numeric(full[i,3])))}
  # simplifying the full table
  full <- cbind(full[,1:3], size, full[,c(4,7,8,5)])
  # defining column names and writing the table out
  colnames(full) <- c("CNV_chrom", "CNV_start", "CNV_end", "CNV_length", "CNV_type", "Gene_name", "Gene_count", "sampleID")
  write.table(full, paste0(s, "_CNV_filt_annot_per_event.tab"), sep="\t", quote=F, col.names = T, row.names = F)
  
  # adding the consequence annotations to the gene-wise split table
  consequence <- matrix(ncol=1, nrow=nrow(split))
  for (i in 1:nrow(split)) {
    # first annotating deletions which overlap 100% with the coding sequence with a loss consequence
    if (isTRUE(split[i,4]=="DEL" & split[i,12]==100)) {consequence[i,1] <- c("loss")}
    # annotating deletions which overlap 1-99% with the coding sequence with a coding_sequence_variant consequence
    if (isTRUE(split[i,4]=="DEL" & split[i,12]>0 & split[i,12]<100)) {consequence[i,1] <- c("coding_sequence_variant")}
    # annotating AMPlications which overlap 100% with the coding sequence with a gain consequence
    if (isTRUE(split[i,4]=="AMP" & split[i,12]==100)) {consequence[i,1] <- c("gain")}
    # annotating AMPlications which overlap 1-99% with the coding sequence with a coding_sequence_variant consequence
    if (isTRUE(split[i,4]=="AMP" & split[i,12]>0 & split[i,12]<100)) {consequence[i,1] <- c("coding_sequence_variant")}
    # annotating deletions which overlap 100% with the transcript of a non-coding gene with a nc_gene_loss consequence
    if (isTRUE(split[i,4]=="DEL" & split[i,12]==0 & split[i,13]=="txStart-txEnd")) {consequence[i,1] <- c("nc_gene_loss")}
    # annotating deletions which overlap 1-99% with the sequence of a non-coding gene with a nc_gene_sequence_variant consequence
    if (isTRUE(split[i,4]=="DEL" & as.numeric(split[i,2])>as.numeric(split[i,10]) & as.numeric(split[i,3])>as.numeric(split[i,11]) & as.numeric(split[i,2])<as.numeric(split[i,11]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR" | split[i,4]=="DEL" & as.numeric(split[i,2])<as.numeric(split[i,10]) & as.numeric(split[i,3])<as.numeric(split[i,11]) & as.numeric(split[i,3])>as.numeric(split[i,10]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR" | split[i,4]=="DEL" & as.numeric(split[i,2])>as.numeric(split[i,10]) & as.numeric(split[i,3])<as.numeric(split[i,11]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR")) {consequence[i,1] <- c("nc_gene_sequence_variant")}
    # annotating AMPlications which overlap 100% with the transcript of a non-coding gene with a nc_gene_gain consequence
    if (isTRUE(split[i,4]=="AMP" & split[i,12]==0 & split[i,13]=="txStart-txEnd")) {consequence[i,1] <- c("nc_gene_gain")}
    # annotating AMPlications which overlap 1-99% with the sequence of a non-coding gene with a nc_gene_sequence_variant consequence
    if (isTRUE(split[i,4]=="AMP" & as.numeric(split[i,2])>as.numeric(split[i,10]) & as.numeric(split[i,3])>as.numeric(split[i,11]) & as.numeric(split[i,2])<as.numeric(split[i,11]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR" | split[i,4]=="AMP" & as.numeric(split[i,2])<as.numeric(split[i,10]) & as.numeric(split[i,3])<as.numeric(split[i,11]) & as.numeric(split[i,3])>as.numeric(split[i,10]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR" | split[i,4]=="AMP" & as.numeric(split[i,2])>as.numeric(split[i,10]) & as.numeric(split[i,3])<as.numeric(split[i,11]) & split[i,13]!="txStart-txEnd" & split[i,12]==0 & split[i,14]=="UTR")) {consequence[i,1] <- c("nc_gene_sequence_variant")}
    # annotating deletions and AMPlications affecting an intronic/UTR regions accordingly
    if (isTRUE(split[i,4]=="DEL" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="CDS")) {consequence[i,1] <- c("intronic")}
    if (isTRUE(split[i,4]=="DEL" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="3'UTR")) {consequence[i,1] <- c("3_UTR_variant")}
    if (isTRUE(split[i,4]=="DEL" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="5'UTR")) {consequence[i,1] <- c("5_UTR_variant")}
    if (isTRUE(split[i,4]=="AMP" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="CDS" & split[i,15]=="no")) {consequence[i,1] <- c("intronic")}
    if (isTRUE(split[i,4]=="AMP" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="CDS" & split[i,15]=="yes")) {consequence[i,1] <- c("coding_sequence_variant")}
    if (isTRUE(split[i,4]=="AMP" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="3'UTR")) {consequence[i,1] <- c("3_UTR_variant")}
    if (isTRUE(split[i,4]=="AMP" & split[i,12]==0 & split[i,13]!="txStart-txEnd" & split[i,14]=="5'UTR")) {consequence[i,1] <- c("5_UTR_variant")}}
  
  # adding the event size to the table
  size <- matrix(ncol=1, nrow=nrow(split))
  for (i in 1:nrow(split)) {
    size[i,1] <- diff(c(as.numeric(split[i,2]), as.numeric(split[i,3])))}
  split <- cbind(split[,c(1:3)], size, split[,4], consequence, split[,c(7, 9:12,15,13,14,5)])
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
