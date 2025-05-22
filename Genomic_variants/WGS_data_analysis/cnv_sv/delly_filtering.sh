#!/bin/bash
# filter Delly results with the blacklist regions and based on read support
# and separate results into indels, SVs and CNVs based on size
SAMPLES="GE0310 GE1312 GE3309 GE3317 GE3320 GE8315 GE8321 GE8322 GE8329 GE9336 GE9344 GE0319 GE0335 GE1329 GE2328 GE5312 GE5338 GE5340 GE6339 GE6345 GE6347 GE7301 GE7304 GE7305 GE7318 GE7324 GE9321 GE9342 GE9347 GE9348 GEK20315 GET19344 GET20301"

# enabling bcftools and bedtools
module load bcftools
module load bedtools

# bcf to vcf
for i in ${SAMPLES}; do bcftools view ${i}.delly.geno.bcf > ${i}.delly.vcf; done

# read support filtering
# filtering based on paired read support, which is in INFO/PE column
# including variants with >= 5 reads supporting them
for i in ${SAMPLES}; do bcftools filter ${i}.delly.vcf --include 'INFO/PE >= 5' --soft-filter 'low_paired_read_support' --mode '+' > ${i}.filtering1.vcf; done

# including variants with >= 5 split reads supporting them (based on INFO/SR column)
for i in ${SAMPLES}; do bcftools filter ${i}.filtering1.vcf --include 'INFO/SR >= 5' --soft-filter 'low_split_read_support' --mode '+' > ${i}.filtering2.vcf; done

# variants that pass all filters
for i in ${SAMPLES}; do bcftools view -f PASS ${i}.filtering2.vcf > ${i}.filtering3.vcf; done

# extracting the breakends from the output to their own file, and the other types to their own
for i in ${SAMPLES}; do bcftools filter ${i}.filtering3.vcf --exclude 'INFO/SVTYPE="BND"' --soft-filter 'breakend' --mode '+' > ${i}.filtering4.vcf; done
for i in ${SAMPLES}; do bcftools view -f PASS ${i}.filtering4.vcf > ${i}.wo_breakends.vcf; done
for i in ${SAMPLES}; do bcftools view -f breakend ${i}.filtering4.vcf > ${i}.breakends.vcf; done

# non-bnd variants of which at least 50% locate onto a blacklist region are filtered out with bedtools
for i in ${SAMPLES}; do bedtools subtract -a ${i}.wo_breakends.vcf -b hg19-blacklist.v2.mod.bed -header -f 0.5 > ${i}.delly_wo_breakends_wo_blacklist.vcf; done

# breakend file is processed separately
# first defining the breakends that will be included in the results, since read support filters are passed
# and at least the 1st breakend does not locate into a blacklist region
for i in ${SAMPLES}; do bedtools subtract -a ${i}.breakends.vcf -b hg19-blacklist.v2.mod.bed -f 0.5 > ${i}.delly_breakends_wo_blacklist.vcf; done
# breakends locating onto the blacklist regions are saved into a separate file
for i in ${SAMPLES}; do bedtools intersect -a ${i}.breakends.vcf -b hg19-blacklist.v2.mod.bed -header > ${i}.delly_breakends_in_blacklist.vcf; done
# second positions are extracted and checked which of them locate outside blacklist regions and thgus should be 'saved' into he filtered result
for i in ${SAMPLES}; do bcftools query -f '%INFO/CHR2\t%INFO/POS2\t%INFO/POS2' ${i}.delly_breakends_in_blacklist.vcf > ${i}.delly_breakends2_in_blacklist.bed; done
# checking which of these 2nd breakends do not locate into blacklist regions, and thus the event should be saved
for i in ${SAMPLES}; do bedtools subtract -a ${i}.delly_breakends2_in_blacklist.bed -b hg19-blacklist.v2.mod.bed -wa > ${i}.breakends_to_save.bed; done
for i in ${SAMPLES}; do cut -f 2 ${i}.breakends_to_save.bed > ${i}.breakend_coords_to_save.txt; done
for i in ${SAMPLES}; do while read p; do grep "$p" ${i}.delly_breakends_in_blacklist.vcf > ${i}.rows_to_be_saved.txt ; done < ${i}.breakend_coords_to_save.txt; done
# adding the filtered breakends to the file with the other variant types
for i in ${SAMPLES}; do cat ${i}.delly_wo_breakends_wo_blacklist.vcf ${i}.delly_breakends_wo_blacklist.vcf ${i}.rows_to_be_saved.txt > ${i}.delly.filtered.vcf; done   

# size filtering
# classifying deletions/insertions less than 150 bp as indels
# and deletions/duplications larger than 1Mb as CNVs
# and all translocations and inversions, as well as insertions and and deletions <150bp and >1Mb, and duplications >1Mb as SVs
for i in ${SAMPLES}; do grep '^#' ${i}.delly.filtered.vcf > ${i}.header.vcf; done
for i in ${SAMPLES}; do awk -v OFS='\t' '{ split($8,a,";SVMETHOD=EMBL.DELLYv0.8.7;END="); split(a[1],b,";SVTYPE="); split(a[2],c,";PE="); if (c[1]-$2 >= 1000000 && b[2]=="DEL" || c[1]-$2 >= 1000000 && b[2] =="DUP") { print }}' ${i}.delly.filtered.vcf > ${i}.CNV.noheader.vcf; done
for i in ${SAMPLES}; do cat ${i}.header.vcf ${i}.CNV.noheader.vcf > ${i}.dellyCNV.filtered.vcf; done   
for i in ${SAMPLES}; do awk -v OFS='\t' '{ split($8,a,";SVMETHOD=EMBL.DELLYv0.8.7;END="); split(a[1],b,";SVTYPE="); split(a[2],c,";PE="); if (c[1]-$2 < 1000000 && c[1]-$2 > 150 && b[2]=="DEL" || c[1]-$2 < 1000000 && c[1]-$2 > 150 && b[2]=="INS" || c[1]-$2 < 1000000 && b[2]=="DUP" || b[2]=="BND" || b[2]=="INV") { print }}' ${i}.delly.filtered.vcf > ${i}.SV.noheader.vcf; done
for i in ${SAMPLES}; do uniq ${i}.SV.noheader.vcf > ${i}.SV.noheader2.vcf; done
for i in ${SAMPLES}; do cat ${i}.header.vcf ${i}.SV.noheader2.vcf > ${i}.dellySV.filtered.vcf; done
for i in ${SAMPLES}; do awk -v OFS='\t' '{ split($8,a,";SVMETHOD=EMBL.DELLYv0.8.7;END="); split(a[1],b,";SVTYPE="); split(a[2],c,";PE="); if (c[1]-$2 < 150 && b[2]=="DEL" || c[1]-$2 < 150 && b[2]=="INS") { print }}' ${i}.delly.filtered.vcf > ${i}.indel.noheader.vcf; done
for i in ${SAMPLES}; do cat ${i}.header.vcf ${i}.indel.noheader.vcf > ${i}.dellyIndels.filtered.vcf; done  

# removing unnecessary intermediate files
rm *header*
rm *filtering*
rm *wo_blacklist*
rm *in_blacklist*
rm *breakend*