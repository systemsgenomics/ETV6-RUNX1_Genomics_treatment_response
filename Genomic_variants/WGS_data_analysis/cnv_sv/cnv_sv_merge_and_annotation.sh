# CNV and SV merge and annotation with AnnotSV

# patient IDs
SAMPLES="GE0310 GE1312 GE3309 GE3317 GE3320 GE8315 GE8321 GE8322 GE8329 GE9336 GE9344 GE0319 GE0335 GE1329 GE2328 GE5312 GE5338 GE5340 GE6339 GE6345 GE6347 GE7301 GE7304 GE7305 GE7318 GE7324 GE9321 GE9342 GE9347 GE9348 GEK20315 GET19344 GET20301"

# merging CNV events (> 1 Mb) from Ascat, Delly, and Manta
# prioritizing Manta breakpoints over Delly over Ascat
for i in ${SAMPLES}; do svdb --merge --no_intra --bnd_distance 5000 --overlap 0.80 --vcf ${i}.mantaCNV.filtered.vcf:manta ${i}.dellyCNV.filtered.vcf:delly ${i}.ascatCNV.filtered.vcf:ascat --priority manta,delly,ascat > ${i}.merged_CNVs.vcf; done

# annotating merged CNV events with AnnotSV
for i in ${SAMPLES}; do AnnotSV -genomeBuild GRCh37 -SVinputFile ${i}.merged_CNVs.vcf -outputFile ${i}.CNV.merged.annotated.tsv; done

# merging SV events (< 1 Mb) from Delly, and Manta
# prioritizing Manta breakpoints over Delly
for i in ${SAMPLES}; do svdb --merge --no_intra --bnd_distance 5000 --overlap 0.80 --vcf ${i}.mantaSV.filtered.vcf:manta ${i}.dellySV.filtered.vcf:delly --priority manta,delly > ${i}.merged_SVs.vcf; done

# annotating merged SV events with AnnotSV
for i in ${SAMPLES}; do AnnotSV -genomeBuild GRCh37 -SVinputFile ${i}.merged_SVs.vcf -outputFile ${i}.SV.merged.annotated.tsv; done

# processing the annotated CNV events

# extracting the relevant columns
SAMPLES1="GE0335 GE3309 GE5312 GE5338 GE5340 GE6339 GE6345 GE7304 GE7305 GE7324 GE8315 GE9321 GE9344"
for i in ${SAMPLES1}; do cut -f 2,3,4,6,7,16,18,20,21,22,25,28 ${i}.CNV.merged.annotated.tsv > ${i}.CNV.merged.intermediate.tsv; done
SAMPLES2="GEK20315 GE9336 GE3320"
for i in ${SAMPLES2}; do cut -f 2,3,4,6,7,17,19,21,22,23,26,29 ${i}.CNV.merged.annotated.tsv > ${i}.CNV.merged.intermediate.tsv; done
SAMPLES3="GE0310 GE1312 GE2328 GE8321 GE8329 GE0319 GE8322 GE1329 GE6347 GE7301 GE7318 GE9342 GE9347 GET19344 GET20301"
for i in ${SAMPLES3}; do cut -f 2,3,4,6,7,18,20,22,23,24,27,30 ${i}.CNV.merged.annotated.tsv > ${i}.CNV.merged.intermediate.tsv; done
cut -f 2,3,4,6,7,19,21,23,24,25,28,31 GE9348.CNV.merged.annotated.tsv > GE9348.CNV.merged.intermediate.tsv
# relpacing empty values with .
for i in ${SAMPLES}; do awk -F'\t' '{ for (i = 1; i <= NF; ++i) sub(/^$/, ".", $i)}1' OFS='\t' ${i}.CNV.merged.intermediate.tsv > ${i}.CNV.merged.intermediate2.tsv; done
# adding the sample ID as one column
for i in ${SAMPLES}; do awk -v OFS='\t' -v SID=${i} '{print SID,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' ${i}.CNV.merged.intermediate2.tsv > ${i}.CNV.merged.simple.tsv; done
# removing the header
for i in ${SAMPLES}; do tail -n +2 ${i}.CNV.merged.simple.tsv > ${i}.CNV.merged.intermediate3.tsv; done
# collecting the info from all samples into one table
for i in ${SAMPLES}; do cat ${i}.CNV.merged.intermediate3.tsv >> collective_CNV_annotated.tab; done
# removing unnecessary intermediate files
rm *intermediate*

# processing the annotated SV events

# extracting the relevant columns
for i in ${SAMPLES}; do cut -f 2,3,4,5,6,7,17,19,21,22,23,26,29 ${i}.SV.merged.annotated.tsv > ${i}.SV.merged.intermediate.tsv; done
# relpacing empty values with .
for i in ${SAMPLES}; do awk -F'\t' '{ for (i = 1; i <= NF; ++i) sub(/^$/, ".", $i)}1' OFS='\t' ${i}.SV.merged.intermediate.tsv > ${i}.SV.merged.intermediate2.tsv; done
# adding the sample ID as one column
for i in ${SAMPLES}; do awk -v OFS='\t' -v SID=${i} '{print SID,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' ${i}.SV.merged.intermediate2.tsv > ${i}.SV.merged.simple.tsv; done
# removing the header
for i in ${SAMPLES}; do tail -n +2 ${i}.SV.merged.simple.tsv > ${i}.SV.merged.intermediate3.tsv; done
# collecting the info from all samples into one table
for i in ${SAMPLES}; do cat ${i}.SV.merged.intermediate3.tsv >> collective_SV_annotated.tab; done
# removing unnecessary intermediate files
rm *intermediate*
