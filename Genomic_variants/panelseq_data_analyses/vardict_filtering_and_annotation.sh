#!/bin/bash
# filter the panelseq Vardict SNV/indel results
SAMPLES="08-204 08-209 08-212 08-219 08-222 08-271 08-286 08-297 08-304 08-307 09-004 09-019 09-024 09-036 09-043 09-051 09-088 09-119 09-120 09-121 09-143 09-162 09-175 09-185 09-191 09-215 09-296 09-308 10-024 10-025 10-026 10-047 10-055 10-059 10-064 10-078 10-095 10-113 10-126 10-127 10-175 10-188 10-191 10-192 10-198 10-210 10-217 10-223 10-232 10-234 10-235 10-263 10-291 10-308 10-331 10-333 10034  10035  10037  10044  10060  10082  10092  10104  10127  10178  10202  10213  10216  10219  10230  10264  10267  10274  10279  10289  10323  10339  10345  10346  10348  10349  10355  10366  10369  10372  10376  10377  10379  10399  10400  10416  10421  10432  10446  10452  10460  10461  10469  10470  10480  10520  10524  10532  10534  10541  10611  10634  10660  10666  10685  10704  10746  10747  10759  10806  10812  10817  10855  10882  10888  10932  10933  10954  10968  10969  11-065 11-067 11-098 11-165 11-171 11-189 11-212 11-236 11-237 11-244 11-245 11-246 11-255 11-257 11-263 11-290 11-299 11-313 11-317 11-328 11013  11016  11019  11025  11042  11253  12-014 12-017 12-033 12-038 12-040 12-063 12-098 12-126 12-128 12-130 12-153 12-154 12-168 12-179 12-180 12-200 12-213 12-215 12-216 12-222 12-225 12-227 12-234 12-236 12-245 12-256 12-258 12-260 12-261 12-264 13-007 13-011 13-016 13-018 13-050 13-061 13-069 13-078 13-092 13-101 13-111 13-112 13-115 13-128 13-131 13-137 13-155 13-171 13-184 13-185 13-187 13-191 13-196 13-197 13-214 13-223 13-226 13-236 13-246 13-270 13-271 13-285 14-023 14-024 14-044 14-046 14-055 14-070 14-074 14-087 14-091 14-103 14-128 14-144 14-161 14-202 14-211 14-220 14-223 14-239 14-250 14-251 14-254 14-258 14-283 14-284 15-001 15-010 15-015 15-019 15-048 15-106 15-123 15-139 15-145 15-148 15-149 15-150 15-155 15-176 15-182 15-208 15-223 15-224 15-228 15-254 15-264 16-004 08-186 15-196"

# enabling bedtools, bcftools and nextflow
module load bedtools
module load bcftools
module load nextflow

# filtering the variants
# including variants found by VarDict passing default filters, and which have allele fraction >= 0.05
for s in ${SAMPLES}; do bcftools filter vcf/SNV.somatic.case${s}.merged.clinical.filtered.pass.vcf --exclude 'INFO/AF < 0.05' --soft-filter 'low_af' --mode '+' > ${s}.intermediate.vcf; done
for s in ${SAMPLES}; do bcftools view -f PASS ${s}.intermediate.vcf > vcf/${s}.filtered.vcf; done
for s in ${SAMPLES}; do bcftools filter --include 'INFO/FOUND_IN[*] ~ "vardict"' vcf/${s}.filtered.vcf > vcf/${s}.vardict_filtered.vcf; done
rm *intermediate*

# VEP annotations with nf-core/sarek
for s in ${SAMPLES}; do bcftools view vcf/${s}.vardict_filtered.vcf -Oz -o vcf/SNV.somatic.case${s}.vardict.vcf.gz; done
nextflow run nf-core/sarek -r 3.3.2  -profile singularity --genome hg19 --vep_cache s3://annotation-cache/vep_cache/110_GRCh37 --max_memory '16.GB' --outdir vep_annot --input samples.csv --step annotate --tools vep

# generating collective tables with the af-filtered variants
for s in ${SAMPLES}; do bcftools +split-vep -d -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/TYPE\t%QUAL\t%INFO/DP\t%INFO/AF\t%INFO/HIAF\t%INFO/ADJAF\t%INFO/VD\t%INFO/HICNT\t%INFO/PMEAN\t%INFO/MQ\t%INFO/VARBIAS\t%INFO/SBF\t%Consequence\t%IMPACT\t%SYMBOL\n' vep_annot/annotation/${s}/SNV.somatic.case${s}.vardict_VEP.ann.vcf.gz > ${s}.af0.05.res.remove.tab; done
for s in ${SAMPLES}; do awk -v OFS='\t' -v SID=${s} '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,SID}' ${s}.af0.05.res.remove.tab > ${s}.af0.05.res.remove2.tab; done
for s in ${SAMPLES}; do cat ${s}.af0.05.res.remove2.tab >> collective.SNV_indel.af0.05.vardict.filtered.tab; done
rm *remove*

# processing the VarDict results of the panel-of-normals samples
SAMPLES_PON="A69 A95 E04 E19 E33 E62 E64 E90 E93 I502 I503 I88 I91 I92 KE101 KE109 KE120 KE123 KI110 KI126 KI130 KI131 KI132 KI137 KI139 E11 E15 E16 E18 E21 E27 E39 E40 E42 E43 E66 I54 I61"
for s in ${SAMPLES_PON}; do bcftools +split-vep -d -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/TYPE\t%QUAL\t%INFO/DP\t%INFO/AF\t%INFO/HIAF\t%INFO/ADJAF\t%INFO/VD\t%INFO/HICNT\t%INFO/PMEAN\t%INFO/MQ\t%INFO/VARBIAS\t%INFO/SBF\t%Consequence\t%IMPACT\t%SYMBOL\n' pon/SNV.somatic.${s}.vardict.clinical.filtered.pass.vcf > ${s}.res.remove.tab; done
for s in ${SAMPLES_PON}; do awk -v OFS='\t' -v SID=${s} '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,SID}' ${s}.res.remove.tab > ${s}.res.remove2.tab; done
for s in ${SAMPLES_PON}; do cat ${s}.res.remove2.tab >> collective_pon_SNV_indel.vardict.tab; done
rm *remove*
