#!/bin/bash
# bedtools shuffle to randomize the CNV coordinates for empirical p-value determination

# enabling bedtools
module load bedtools

# patient samples with CNV data
SAMPLES="15-015 12-168 13-131 10264 GE0310   GE0319   GE0335   GE1312   GE1329   GE2328   GE3309   GE3317   GE3320   GE5312   GE5338   GE5340   GE6339  GE6345   GE6347   GE7301   GE7304   GE7305   GE7318   GE7324   GE8315   GE8321   GE8322   GE8329   GE9321   GE9336  GE9342   GE9344   GE9347   GE9348   GEK20315 GET19344 GET20301 08-209   09-051   09-143   09-308   10-198   11-065  11-246   12-260   13-196   14-223   10037    16-004   10339    10379    10460    10532    10432    09-175   10524   10806    10377    10349    10355    11-328   10611    11-237   13-285   10-232   11013    09-296   10202    11-257  09-019   10969    10759    10888    11-313   08-271   14-074   11042    10541    10421    11-171   15-019   14-103  10480    14-144   15-015   10-191   10-188   09-088   11-299   14-283   09-024   13-223   10452    11-255   09-162  09-120   09-043   10932    12-033   12-153   12-236   10666    10267    11253    13-092   10289    10376    08-186  09-004   13-007   14-024   10216    10034    10399    10369    13-271   15-148   14-055   12-213   12-168   13-111  10219    15-224   15-176   11025    13-112   10-055   13-226   10634    15-264   10-308   08-212   15-123   12-040  12-215   13-011   10704    13-270   10230    10416    11019    15-182   12-245   12-017   11-212  09-121  14-070   15-145   13-246   10-059   11-098   10-234   13-128   10-291   13-236   10-235   10469    10-333   13-214  13-187   10968    12-225   12-126   15-010   10855    08-304   10372    15-106   10-113   13-069   12-038   10-223  10213    10954    10274    14-254   11-317   13-185   14-091   12-180   10-263   12-014   10092    12-216   15-223  13-131   09-185   14-046   11-290   10-192   12-154   10279    10812    10-175   12-264   11-263   12-256   11-067  10264    10-127   10520    14-161   15-254   12-063   10345    14-202   10-126   10470    13-197   14-211   11-245  10104    11-236   15-001   10685    14-284   12-200   10035    08-219   15-149   10-025   13-115   13-184   15-228  10366    13-171   13-018   10400    13-137   11-244   10746    13-050   15-155   10534    15-139   14-087   12-227  10882    10446    12-128   12-234   15-150   10-210   12-130   13-078   14-044   09-215   10127    11016    14-251  14-023   10346    08-297   08-286   08-204   09-036   10323    10060    10-331   10-024   10044    10-064   15-196  12-261   14-258   10-217   11-165  10178    11-189   14-239   10660    12-179   15-208   09-119   10747   12-222 13-061   13-016   10-026   13-191   15-048   13-155   10082    10348    14-128   14-220   10817   12-258   10461    10-078   13-101   09-191   14-250   10933    ALL_249  ALL_269  ALL_251  ALL_284  ALL_294  ALL_301 ALL_330  ALL_346  ALL_309  ALL_305  ALL_354  ALL_326  ALL_337  ALL_322  ALL_361  ALL_444  ALL_364  ALL_448  ALL_440 ALL_450  ALL_381  ALL_434  ALL_446  ALL_439  ALL_431  ALL_430  ALL_384  ALL_386  ALL_394  ALL_400  ALL_525  ALL_458 ALL_504  ALL_510  ALL_511  ALL_512  ALL_488  ALL_493  ALL_498  ALL_500  ALL_541  ALL_540  ALL_589  ALL_581  ALL_574 ALL_586  ALL_582  ALL_548  ALL_550  ALL_551  ALL_731  ALL_732  ALL_690  ALL_716  ALL_701  ALL_706  ALL_720  ALL_906 ALL_936  ALL_949  ALL_1027 ALL_1051 ALL_1060 ALL_1069 ALL_1071 ALL_1085 ALL_1098 ALL_1123 ALL_1158 ALL_398  ALL_755 ALL_758"

## chromosome 3 deletions

# subsetting the deletion data to chr3 deletions
for i in ${SAMPLES}; do awk '$1==3' input/${i}.DELs.bed > input/${i}.3_DELs.bed; done

# shuffling coordinates for 10000 times
for i in ${SAMPLES}; do for s in {1..10000}; do bedtools shuffle -i input/${i}.3_DELs.bed -g input/genome.bed -chrom -noOverlapping > shuffled_output/${i}.3_DELs_shuffle$s.bed; done; done

# collective files for each shuffle round
for s in {1..10000}; do cat shuffled_output/*.3_DELs_shuffle$s.bed >> shuffled_output/collective_3_DELs_shuffle$s.bed; done

## chromosome 6q deletions

# subsetting the deletion data to chr6 q arm deletions
for i in ${SAMPLES}; do awk '$1==6' input/${i}.DELs.bed > input/${i}.6_DELs.bed; done
for i in ${SAMPLES}; do awk '$2 >= 60000000' input/${i}.6_DELs.bed > input/${i}.6q_DELs.bed; done

# shuffling coordinates for 10000 times
for i in ${SAMPLES}; do for s in {1..10000}; do bedtools shuffle -i input/${i}.6q_DELs.bed -g input/genome.bed -chrom -noOverlapping -excl input/exclude.bed > shuffled_output/${i}.6q_DELs_shuffle$s.bed; done; done

# collective files for each shuffle round
for s in {1..10000}; do cat shuffled_output/*.6q_DELs_shuffle$s.bed >> shuffled_output/collective_6q_DELs_shuffle$s.bed; done

## chromosome 9 deletions

# subsetting the deletion data to chr9 deletions
for i in ${SAMPLES}; do awk '$1==9' input/${i}.DELs.bed > input/${i}.9_DELs.bed; done

# shuffling coordinates for 10000 times
for i in ${SAMPLES}; do for s in {1..10000}; do bedtools shuffle -i input/${i}.9_DELs.bed -g input/genome.bed -chrom -noOverlapping > shuffled_output/${i}.9_DELs_shuffle$s.bed; done; done

# collective files for each shuffle round
for s in {1..10000}; do cat shuffled_output/*.9_DELs_shuffle$s.bed >> shuffled_output/collective_9_DELs_shuffle$s.bed; done

## chromosome 10 amplifications

# subsetting the amplification data to chr9 amplifications
for i in ${SAMPLES}; do awk '$1==10' input/${i}.AMPs.bed > input/${i}.10_AMPs.bed; done

# shuffling coordinates for 10000 times
for i in ${SAMPLES}; do for s in {1..10000}; do bedtools shuffle -i input/${i}.10_AMPs.bed -g input/genome.bed -chrom -noOverlapping > shuffled_output/${i}.10_AMPs_shuffle$s.bed; done; done

# collective files for each shuffle round
for s in {1..10000}; do cat shuffled_output/*.10_AMPs_shuffle$s.bed >> shuffled_output/collective_10_AMPs_shuffle$s.bed; done

## chromosome 11 deletions

# subsetting the deletion data to chr11 deletions
for i in ${SAMPLES}; do awk '$1==11' input/${i}.DELs.bed > input/${i}.11_DELs.bed; done
for i in ${SAMPLES}; do awk '$2 >= 55000000' input/${i}.11_DELs.bed > input/${i}.11q_DELs.bed; done

# shuffling coordinates for 10000 times
for i in ${SAMPLES}; do for s in {1..10000}; do bedtools shuffle -i input/${i}.11q_DELs.bed -g input/genome.bed -chrom -noOverlapping -excl input/exclude.bed > shuffled_output/${i}.11q_DELs_shuffle$s.bed; done; done

# collective files for each shuffle round
for s in {1..10000}; do cat shuffled_output/*.11q_DELs_shuffle$s.bed >> shuffled_output/collective_11q_DELs_shuffle$s.bed; done

## chromosome 12 amplifications

# subsetting the amplification data to chr12 amplifications
for i in ${SAMPLES}; do awk '$1==12' input/${i}.AMPs.bed > input/${i}.12_AMPs.bed; done
for i in ${SAMPLES}; do awk '$3 <= 15000000' input/${i}.12_AMPs.bed > input/${i}.12p_AMPs.bed; done

# collective file for all amplifications
for i in ${SAMPLES}; do cat input/${i}.12p_AMPs.bed >> input/collective_12q_AMPs.bed; done

# copying the file 10000 times to make mock shuffle files
for s in {1..10000}; do cp input/collective_12q_AMPs.bed shuffled_output/collective_12p_AMPs_shuffle$s.bed; done

## chromsome 12 p arm deletions

# subsetting the deletion data to chr12p deletions
for i in ${SAMPLES}; do awk '$1==12' input/${i}.DELs.bed > input/${i}.12_DELs.bed; done
for i in ${SAMPLES}; do awk '$3 <= 35000000' input/${i}.12_DELs.bed > input/${i}.12p_DELs.bed; done

# shuffling coordinates for 10000 times
for i in ${SAMPLES}; do for s in {1..10000}; do bedtools shuffle -i input/${i}.12p_DELs.bed -g input/genome.bed -chrom -noOverlapping > shuffled_output/${i}.12p_DELs_shuffle$s.bed; done; done

# collective files for each shuffle round
for s in {1..10000}; do cat shuffled_output/*12p_DELs_shuffle$s.bed >> shuffled_output/collective_12p_DELs_shuffle$s.bed; done

## chromsome 15 deletions

# subsetting the deletion data to chr15 deletions
for i in ${SAMPLES}; do awk '$1==15' input/${i}.DELs.bed > input/${i}.15_DELs.bed; done
for i in ${SAMPLES}; do awk '$2 >= 30000000' input/${i}.15_DELs.bed > input/${i}.15q_DELs.bed; done

# shuffling coordinates for 10000 times
for i in ${SAMPLES}; do for s in {1..10000}; do bedtools shuffle -i input/${i}.15q_DELs.bed -g input/genome.bed -chrom -noOverlapping -excl input/exclude.bed > shuffled_output/${i}.15q_DELs_shuffle$s.bed; done; done

# collective files for each shuffle round
for s in {1..10000}; do cat shuffled_output/*.15q_DELs_shuffle$s.bed >> shuffled_output/collective_15q_DELs_shuffle$s.bed; done

## chromsome 16 amplifications

# subsetting the amplification data to chr16 amplifications
for i in ${SAMPLES}; do awk '$1==16' input/${i}.AMPs.bed > input/${i}.16_AMPs.bed; done

# shuffling coordinates for 10000 times
for i in ${SAMPLES}; do for s in {1..10000}; do bedtools shuffle -i input/${i}.16_AMPs.bed -g input/genome.bed -chrom -noOverlapping > shuffled_output/${i}.16_AMPs_shuffle$s.bed; done; done

# collective files for each shuffle round
for s in {1..10000}; do cat shuffled_output/*.16_AMPs_shuffle$s.bed >> shuffled_output/collective_16_AMPs_shuffle$s.bed; done

## chromsome 21 q arm amplifications

# subsetting the amplification data to chr21q amplifications
for i in ${SAMPLES}; do awk '$1==21' input/${i}.AMPs.bed > input/${i}.21_AMPs.bed; done
for i in ${SAMPLES}; do awk '$2 >= 9400000' input/${i}.21_AMPs.bed > input/${i}.21q_AMPs.bed; done

# shuffling coordinates for 10000 times
for i in ${SAMPLES}; do for s in {1..10000}; do bedtools shuffle -i input/${i}.21q_AMPs.bed -g input/genome.bed -chrom -noOverlapping -excl input/exclude.bed > shuffled_output/${i}.21q_AMPs_shuffle$s.bed; done; done

# collective files for each shuffle round
for s in {1..10000}; do cat shuffled_output/*.21q_AMPs_shuffle$s.bed >> shuffled_output/collective_21q_AMPs_shuffle$s.bed; done
