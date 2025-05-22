# blacklist filtering of the conumee CNV and SV results
# and AnnotSV annotation

# case IDs
SAMPLES="ALL_1085 ALL_294 ALL_346 ALL_439 ALL_440 ALL_488 ALL_498 ALL_504 ALL_525 ALL_540 ALL_541 ALL_586 ALL_589 ALL_701 ALL_705 ALL_706 ALL_720 ALL_731 ALL_732 ALL_755 ALL_758 ALL_760 ALL_249 ALL_251 ALL_284 ALL_309 ALL_322 ALL_326 ALL_330 ALL_337 ALL_354 ALL_361 ALL_364 ALL_381 ALL_394 ALL_430 ALL_431 ALL_434 ALL_444 ALL_448 ALL_458 ALL_500 ALL_510 ALL_511 ALL_548 ALL_550 ALL_582 ALL_695 ALL_697 ALL_709 ALL_715 ALL_759 ALL_716 ALL_269 ALL_301 ALL_305 ALL_384 ALL_386 ALL_398 ALL_400 ALL_446 ALL_450 ALL_493 ALL_512 ALL_551 ALL_574 ALL_581 ALL_698 ALL_871 ALL_883 ALL_886 ALL_887 ALL_893 ALL_895 ALL_916 ALL_922 ALL_928 ALL_933 ALL_949 ALL_958 ALL_963 ALL_1012 ALL_1013 ALL_1024 ALL_1027 ALL_1112 ALL_1119 ALL_1134 ALL_1144 ALL_1158 ALL_1180 ALL_905 ALL_937 ALL_980 ALL_981 ALL_995 ALL_997 ALL_998 ALL_1014 ALL_1020 ALL_1023 ALL_1025 ALL_1032 ALL_1042 ALL_1052 ALL_1071 ALL_1098 ALL_1104 ALL_1115 ALL_1120 ALL_1123 ALL_1126 ALL_1133 ALL_1143 ALL_1146 ALL_1147 ALL_1168 ALL_1172 ALL_1181 ALL_1182 ALL_1186 ALL_1192 ALL_917 ALL_971 ALL_989 ALL_906 ALL_1051 ALL_881 ALL_890 ALL_915 ALL_936 ALL_939 ALL_943 ALL_972 ALL_977 ALL_983 ALL_1021 ALL_1060 ALL_1069 ALL_1163 ALL_1184 ALL_1191 ALL_970"

# enabling bedtools
module load bedtools

# excluding CNVs that overlap 50% or more with the blacklist regions
for i in ${SAMPLES}; do bedtools subtract -a ${i}.conumee_CNVs.bed -b hg19-blacklist.v2.mod.bed -f 0.5 -A > ${i}.conumee_CNVs.filtered.bed; done

# excluding SVs that overlap 50% or more with the blacklist regions
for i in ${SAMPLES}; do bedtools subtract -a ${i}.conumee_SVs.bed -b hg19-blacklist.v2.mod.bed -f 0.5 -A > ${i}.conumee_SVs.filtered.bed; done

# generating collective tables from each variant type
for i in ${SAMPLES}; do cat ${i}.conumee_CNVs.filtered.bed >> collective_conumee_CNVs.bed; done
for i in ${SAMPLES}; do cat ${i}.conumee_SVs.filtered.bed >> collective_conumee_SVs.bed; done

# AnnotSV annotations for CNVs
for s in ${SAMPLES}; do AnnotSV/bin/AnnotSV -SVinputFile ${s}.conumee_CNVs.filtered.bed -outputFile ${s}.CNV.annotated -svtBEDcol 5 -REreport 1 -genomeBuild GRCh37; done
