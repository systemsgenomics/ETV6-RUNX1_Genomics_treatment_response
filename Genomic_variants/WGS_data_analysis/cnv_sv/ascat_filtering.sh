#!/bin/bash
# filter ASCAT results with the blacklist regions and based on size

# enabling bedtools and bcftools
module load bedtools
module load bcftools

# remove vcf header
for i in GE0310 GE0319 GE0335 GE1312 GE1329 GE2328 GE3309 GE3317 GE3320 GE8315 GE8321 GE8322 GE8329 GE9336 GE9344 GE9347 GE9348 GEK20315 GET20301; do grep -v "^#" ./${i}/TUMOR.copynumber.caveman.vcf > ${i}.copynumber.caveman.noheader.vcf; done
for i in GE5312 GE5338 GE5340 GE6339 GE6345 GE6347 GE7301 GE7304 GE7305 GE7318 GE7324 GE9321 GE9342 GET19344; do grep -v "^#" ./${i}/${i}_tumor_R.copynumber.caveman.vcf > ${i}.copynumber.caveman.noheader.vcf; done

# exclude regions with normal copy number marked as .:2:1
SAMPLES="GE0310 GE0319 GE0335 GE1312 GE1329 GE2328 GE3309 GE3317 GE3320 GE5312 GE5338 GE5340 GE6339 GE6345 GE6347 GE7301 GE7304 GE7305 GE7318 GE7324 GE8315 GE8321 GE8322 GE8329 GE9321 GE9336 GE9342 GE9344 GE9347 GE9348 GEK20315 GET19344 GET20301"
for i in ${SAMPLES}; do awk -v OFS='\t' '{if ($11 !="./.:2:1") { print}}' ${i}.copynumber.caveman.noheader.vcf > ${i}.copynumber.caveman.noheader.tumor.vcf; done

# separate CNVs (size over 1 Mb)
for i in ${SAMPLES}; do awk -v OFS='\t' '{ split($8,a,"END=");  if (a[2]-$2 >= 1000000) { print }}' ${i}.copynumber.caveman.noheader.tumor.vcf > ${i}.ascat.CNV.vcf; done 

# separate SVs (size less than 1 Mb)
for i in ${SAMPLES}; do awk -v OFS='\t' '{ split($8,a,"END=");  if (a[2]-$2 < 1000000) { print }}' ${i}.copynumber.caveman.noheader.tumor.vcf > ${i}.ascat.SV.vcf; done

# extract vcf headers
for i in GE0310 GE0319 GE0335 GE1312 GE1329 GE2328 GE3309 GE3317 GE3320 GE8315 GE8321 GE8322 GE8329 GE9336 GE9344 GE9347 GE9348 GEK20315 GET20301; do bcftools view --header-only ./${i}/TUMOR.copynumber.caveman.vcf > ${i}.header.vcf; done
for i in GE5312 GE5338 GE5340 GE6339 GE6345 GE6347 GE7301 GE7304 GE7305 GE7318 GE7324 GE9321 GE9342 GET19344; do bcftools view --header-only ./${i}/${i}_tumor_R.copynumber.caveman.vcf > ${i}.header.vcf; done

# add vcf headers
for i in ${SAMPLES}; do cat ${i}.header.vcf ${i}.ascat.CNV.vcf > ${i}.ascat.CNV.header.vcf; done
for i in ${SAMPLES}; do cat ${i}.header.vcf ${i}.ascat.SV.vcf > ${i}.ascat.SV.header.vcf; done

# from GE9348 sample we will use both 1st and refitted ascat resuls because those present different clones
grep -v "^#" ./ascat-refit/GE9348_refit/TUMOR.copynumber.caveman.vcf > GE9348_refit.copynumber.caveman.noheader.vcf
awk -v OFS='\t' '{if ($11 !="./.:2:1") { print}}' GE9348_refit.copynumber.caveman.noheader.vcf > GE9348_refit.copynumber.caveman.noheader.tumor.vcf
awk -v OFS='\t' '{ split($8,a,"END=");  if (a[2]-$2 >= 1000000) { print }}' GE9348_refit.copynumber.caveman.noheader.tumor.vcf > GE9348_refit.ascat.CNV.vcf
awk -v OFS='\t' '{ split($8,a,"END=");  if (a[2]-$2 < 1000000) { print }}' GE9348_refit.copynumber.caveman.noheader.tumor.vcf > GE9348_refit.ascat.SV.vcf
bcftools view --header-only ./ascat-refit/GE9348_refit/TUMOR.copynumber.caveman.vcf > GE9348_refit.header.vcf
cat GE9348_refit.header.vcf GE9348_refit.ascat.CNV.vcf > GE9348_refit.ascat.CNV.header.vcf
cat GE9348_refit.header.vcf GE9348_refit.ascat.SV.vcf > GE9348_refit.ascat.SV.header.vcf

# same steps for refitted results for GE5338 GEK20315 GE7301 GE3309
SAMPLES_refit="GE5338_refit GEK20315_refit GE7301_refit GE3309_refit"
for i in ${SAMPLES_refit}; do grep -v "^#" ./ascat-refit/${i}/TUMOR.copynumber.caveman.vcf > ${i}.copynumber.caveman.noheader.vcf; done
for i in ${SAMPLES_refit}; do awk -v OFS='\t' '{if ($11 !="./.:2:1") { print}}' ${i}.copynumber.caveman.noheader.vcf > ${i}.copynumber.caveman.noheader.tumor.vcf; done
for i in ${SAMPLES_refit}; do awk -v OFS='\t' '{ split($8,a,"END=");  if (a[2]-$2 >= 1000000) { print }}' ${i}.copynumber.caveman.noheader.tumor.vcf > ${i}.ascat.CNV.vcf; done
for i in ${SAMPLES_refit}; do awk -v OFS='\t' '{ split($8,a,"END=");  if (a[2]-$2 < 1000000) { print }}' ${i}.copynumber.caveman.noheader.tumor.vcf > ${i}.ascat.SV.vcf; done
for i in ${SAMPLES_refit}; do bcftools view --header-only ./ascat-refit/${i}/TUMOR.copynumber.caveman.vcf > ${i}.header.vcf; done
for i in ${SAMPLES_refit}; do cat ${i}.header.vcf ${i}.ascat.CNV.vcf > ${i}.ascat.CNV.header.vcf; done
for i in ${SAMPLES_refit}; do cat ${i}.header.vcf ${i}.ascat.SV.vcf > ${i}.ascat.SV.header.vcf; done

# filter results with blacklist
# by discarding events if they overlap 50% or more with a blacklist region
SAMPLES="GE0310 GE0319 GE0335 GE1312 GE1329 GE2328 GE3309 GE3317 GE3320 GE5312 GE5338 GE5340 GE6339 GE6345 GE6347 GE7301 GE7304 GE7305 GE7318 GE7324 GE8315 GE8321 GE8322 GE8329 GE9321 GE9336 GE9342 GE9344 GE9347 GE9348 GEK20315 GET19344 GET20301 GE9348_refit"
for i in ${SAMPLES}; do bedtools subtract -a ${i}.ascat.CNV.header.vcf -b hg19-blacklist.v2.mod.bed -f 0.5 > ./blacklist_filtered/${i}.ascat.CNV.bl.filtered.vcf; done
for i in ${SAMPLES}; do cat .${i}.header.vcf ./blacklist_filtered/${i}.ascat.CNV.bl.filtered.vcf > ./blacklist_filtered/${i}.ascatCNV.filtered.vcf; done
for i in ${SAMPLES_refit}; do bedtools subtract -a ${i}.ascat.CNV.header.vcf -b ../hg19-blacklist.v2.mod.bed -f 0.5 > ./blacklist_filtered/${i}.ascat.CNV.bl.filtered.vcf; done
for i in ${SAMPLES_refit}; do cat .${i}.header.vcf ./blacklist_filtered/${i}.ascat.CNV.bl.filtered.vcf > ./blacklist_filtered/${i}.ascatCNV.filtered.vcf; done

# extracting rows which contain LOH events, i.e., 2:0 copy number
SAMPLES="GE7304 GE0310 GE0319 GE0335 GE1312 GE1329 GE2328 GE3317 GE3320 GE8315 GE8321 GE8322 GE8329 GE9336 GE9344 GE9347 GE5312 GE5340 GE6339 GE6345 GE6347 GE7305 GE7318 GE7324 GE9321 GE9342 GE9348 GET20301"
SAMPLES2="GE3309 GE5338 GE7301 GEK20315"
for i in ${SAMPLES}; do grep "2:0" ${i}.ascat.CNV.bl.filtered.270422.vcf > ${i}.loh_events.tab; done
for i in ${SAMPLES2}; do grep "2:0" ${i}_refit.ascat.CNV.bl.filtered.270422.vcf > ${i}.loh_events.tab; done
grep "2:0" GET19344_new.ascat.CNV.bl.filtered.270422.vcf > GET19344.loh_events.tab
grep "2:0" GE9348_refit.ascat.CNV.bl.filtered.270422.vcf > GE9348_refit.loh_events.tab

# adding the sample ID into a column
for i in ${SAMPLES}; do awk -v OFS='\t' -v SID=${i} '{print $1,$2,$8,$10,$11,SID}' ${i}.loh_events.tab > ${i}.loh_events.wid.tab; done
for i in ${SAMPLES2}; do awk -v OFS='\t' -v SID=${i} '{print $1,$2,$8,$10,$11,SID}' ${i}.loh_events.tab > ${i}.loh_events.wid.tab; done
awk -v OFS='\t' -v SID="GET19344" '{print $1,$2,$8,$10,$11,SID}' GET19344.loh_events.tab > GET19344.loh_events.wid.tab
awk -v OFS='\t' -v SID="GE9348.B" '{print $1,$2,$8,$10,$11,SID}' GE9348_refit.loh_events.tab > GE9348.B.loh_events.wid.tab
# combining the tables
cat *wid* > collective_loh_events.tab
uniq collective_loh_events.tab > collective_loh_events2.tab
mv collective_loh_events2.tab collective_loh_events.tab
rm *wid*
