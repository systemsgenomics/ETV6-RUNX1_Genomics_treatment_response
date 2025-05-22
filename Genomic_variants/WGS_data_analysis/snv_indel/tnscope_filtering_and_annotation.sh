#!/bin/bash
# TNscope SNV/indel result filtering

# enabling bcftools, bedtools, and nextflow
module load bcftools
module load bedtools
module load nextflow

SAMPLES="GE0310 GE1312 GE3309 GE3317 GE3320 GE8315 GE8321 GE8322 GE8329 GE9336 GE9344 GE0319 GE0335 GE1329 GE2328 GE5312 GE5338 GE5340 GE6339 GE6345 GE6347 GE7301 GE7304 GE7305 GE7318 GE7324 GE9321 GE9342 GE9347 GE9348 GEK20315 GET19344 GET20301"

# filtering based on tumor depth (>= 10)
for i in ${SAMPLES}; do bcftools filter SNV.somatic.${i}.tnscope.pass.vcf --include 'SUM(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= FORMAT/AFDP || SUM(FORMAT/AD[1:0]+FORMAT/AD[1:1]) >= 10' --soft-filter 'balsamic_low_tumor_dp' --mode '+' > ${i}.filtering1.vcf; done
# allelic depth >= 3
for i in ${SAMPLES}; do bcftools filter ${i}.filtering1.vcf --include 'FORMAT/AD[0:1] >= 3' --soft-filter 'balsamic_low_tumor_ad' --mode '+' > ${i}.filtering2.vcf; done
# allele fraction >= 0.05 & allele fration != 1
for i in ${SAMPLES}; do bcftools filter ${i}.filtering2.vcf --include 'FORMAT/AF[0] >= 0.05' --soft-filter 'balsamic_low_af' --mode '+' > ${i}.filtering3.vcf; done
for i in ${SAMPLES}; do bcftools filter ${i}.filtering3.vcf --include 'FORMAT/AF[0] <  1' --soft-filter 'balsamic_af_one' --mode '+' > ${i}.filtering4.vcf; done
# QUAL value >= 40 and/or PV2 value <= 0.05
for i in ${SAMPLES}; do bcftools filter ${i}.filtering4.vcf --include 'QUAL >= 40 | INFO/PV2 <= 0.05' --soft-filter 'low_qual_or_high_pv2' --mode '+' > ${i}.filtering5.vcf; done
# SOR <= 3
for i in ${SAMPLES}; do bcftools filter ${i}.filtering5.vcf --include 'INFO/SOR <= 3' --soft-filter 'high_sor' --mode '+' > ${i}.filtering6.vcf; done

# including the variants that pass all the aforementioned filters in final output
for i in ${SAMPLES}; do bcftools view -f PASS ${i}.filtering6.vcf > SNV.somatic.${i}.tnscope.filtering.vcf; done

# excluding the variants locating onto blacklist regions
for i in ${SAMPLES}; do bedtools subtract -a SNV.somatic.${i}.tnscope.filtering.vcf -b hg19_blacklist.bed -header > SNV.somatic.${i}.tnscope.filtered.vcf; done

# VEP annotation with nf-core/sarek
nextflow run nf-core/sarek -r 3.3.2 -profile singularity --genome hg19 --vep_cache s3://annotation-cache/vep_cache/110_GRCh37 --max_memory '16.GB' --outdir results --fasta human_g1k_v37.fasta --fasta_fai human_g1k_v37.fasta.fai --input samples.csv --step annotate --tools vep

# generating tables from the results
for i in ${SAMPLES}; do bcftools +split-vep -d -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/TLOD\t%VARIANT_CLASS\t%IMPACT\t%Consequence\t%SYMBOL[\t%AF][\t%AFDP]\n' results/annotation/${i}/SNV.somatic.${i}.tnscope.filtered.vcf.gz > SNV.somatic.${i}.tnscope.filtered.tab; done

# collective table with variants from all patient samples
for i in ${SAMPLES}; do awk -v OFS='\t' -v SID=${i} '{print SID,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' SNV.somatic.${i}.tnscope.filtered.tab > ${i}.VEP_annot.tab; done
for i in ${SAMPLES}; do cat ${i}.VEP_annot.tab >> collective_VEP_ann_moderate.tab; done

# removing unnecessary intermediate files
rm *filtering*
rm *VEP_annot*