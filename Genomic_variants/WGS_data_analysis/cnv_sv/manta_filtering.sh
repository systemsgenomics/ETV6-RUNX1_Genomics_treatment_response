#!/bin/bash
# filter Manta results with the blacklist regions and based on read support
# and separate results into indels, SVs and CNVs based on size
SAMPLES="GE0310 GE1312 GE3309 GE3317 GE3320 GE8315 GE8321 GE8322 GE8329 GE9336 GE9344 GE0319 GE0335 GE1329 GE2328 GE5312 GE5338 GE5340 GE6339 GE6345 GE6347 GE7301 GE7304 GE7305 GE7318 GE7324 GE9321 GE9342 GE9347 GE9348 GEK20315 GET19344 GET20301"

# enabling bcftools and bedtools
module load bcftools
module load bedtools

# variants that pass default Manta filters
for i in ${SAMPLES}; do bcftools view -f PASS SV.somatic.${i}.manta.vcf > SV.somatic.${i}.manta.pass.vcf; done

# variants of which at least 50% locate onto a blacklist region are filtered out with bedtools
for i in ${SAMPLES}; do bedtools subtract -a SV.somatic.${i}.manta.pass.vcf -b hg19-blacklist.v2.mod.bed -header -f 0.5 > ${i}.manta_wo_blacklist.vcf; done
# variants locating onto blacklist regions saved into a separate file
for i in ${SAMPLES}; do bedtools intersect -a SV.somatic.${i}.manta.pass.vcf -b hg19-blacklist.v2.mod.bed -header -f 0.5 > ${i}.manta_in_blacklist.vcf; done

# filter based on read support
# including variants with at least 5 paired and 5 split supporting reads
for i in ${SAMPLES}; do bcftools filter ${i}.manta_wo_blacklist.vcf --include 'FORMAT/PR[1:1]  >= 5' --soft-filter 'low_paired_read_support' --mode '+' > ${i}.manta.filtering1.vcf; done
for i in ${SAMPLES}; do bcftools filter ${i}.manta.filtering1.vcf --include 'FORMAT/SR[1:1]  >= 5' --soft-filter 'low_split_read_support' --mode '+' > ${i}.manta.filtering2.vcf; done
for i in ${SAMPLES}; do bcftools view -f PASS ${i}.manta.filtering2.vcf > ${i}.manta.pass.vcf; done
# saving the variants filtered out based on read support as a separate file
for i in ${SAMPLES}; do bcftools view -f low_paired_read_support,low_split_read_support ${i}.manta.filtering2.vcf > ${i}.manta.notpass.vcf; done

# saving variants locating onto balcklist regions and7or filtered out based on read support, if the mate breakend passes all of the filters
for i in ${SAMPLES}; do bcftools query -f '%INFO/MATEID\n' ${i}.manta.pass.vcf > ${i}.mateid.txt; done   
for i in ${SAMPLES}; do awk '$1!="."' ${i}.mateid.txt > ${i}.mateid.final.txt; done
for i in ${SAMPLES}; do while read p; do grep "$p" ${i}.manta_in_blacklist.vcf > ${i}.rows_to_be_saved.txt ; done < ${i}.mateid.final.txt; done
for i in ${SAMPLES}; do while read p; do grep "$p" ${i}.manta.notpass.vcf > ${i}.rows_to_be_saved2.txt ; done < ${i}.mateid.final.txt; done
for i in ${SAMPLES}; do cat ${i}.manta.pass.vcf ${i}.rows_to_be_saved.txt ${i}.rows_to_be_saved2.txt > ${i}.manta.filtered.vcf; done

# size filtering
# classifying deletions/insertions less than 150 bp as indels
# and deletions/duplications larger than 1Mb as CNVs
# and all translocations and inversions, as well as insertions and and deletions <150bp and >1Mb, and duplications >1Mb as SVs
for i in ${SAMPLES}; do grep '^#' ${i}.manta.filtered.vcf > ${i}.header.vcf; done
for i in ${SAMPLES}; do awk -v OFS='\t' '{ split($8,a,"END="); split(a[2],b,";SVTYPE"); if (b[1]-$2 >= 1000000 && $5=="<DEL>" || b[1]-$2 >= 1000000 && $5=="<DUP:TANDEM>") { print }}' ${i}.manta.filtered.vcf > ${i}.CNV.noheader.vcf; done
for i in ${SAMPLES}; do cat ${i}.header.vcf ${i}.CNV.noheader.vcf > ${i}.mantaCNV.filtered.vcf; done   
for i in ${SAMPLES}; do awk -v OFS='\t' '{ split($8,a,"END="); split(a[2],b,";SVTYPE="); split(b[2],c,";SVLEN"); if (b[1]-$2 < 1000000 && b[1]-$2 > 150 && c[1]=="DEL" || b[1]-$2 < 1000000 && b[1]-$2 > 150 && c[1]=="INS" || b[1]-$2 < 1000000 && c[1]=="DUP") { print }}' ${i}.manta.filtered.vcf > ${i}.SV1.noheader.vcf; done
for i in ${SAMPLES}; do awk -v OFS='\t' '{ split($8,a,"SVTYPE="); split(a[2],b,";MATEID"); if (b[1]=="BND") { print }}' ${i}.manta.filtered.vcf > ${i}.SV2.noheader.vcf; done
for i in ${SAMPLES}; do awk -v OFS='\t' '{ split($8,a,"SVTYPE="); split(a[2],b,";SVLEN"); if (b[1]=="INV") { print }}' ${i}.manta.filtered.vcf > ${i}.SV3.noheader.vcf; done
for i in ${SAMPLES}; do cat ${i}.header.vcf ${i}.SV1.noheader.vcf ${i}.SV2.noheader.vcf ${i}.SV3.noheader.vcf > ${i}.mantaSV.filtered.vcf; done
for i in ${SAMPLES}; do awk -v OFS='\t' '{ split($8,a,"END="); split(a[2],b,";SVTYPE="); split(b[2],c,";SVLEN"); if (b[1]-$2 < 150 && c[1]=="DEL" || b[1]-$2 < 150 && c[1]=="INS") { print }}' ${i}.manta.filtered.vcf > ${i}.indel.noheader.vcf; done
for i in ${SAMPLES}; do cat ${i}.header.vcf ${i}.indel.noheader.vcf > ${i}.mantaIndels.filtered.vcf; done  

# removing unnecessary intermediate files
rm *SV1*
rm *SV2*
rm *SV3*
rm *noheader*
rm *header*
rm *mateid*
rm *rows_to_be_saved*
rm *pass*
rm *filtering*
rm *wo_blacklist*
rm *in_blacklist*
