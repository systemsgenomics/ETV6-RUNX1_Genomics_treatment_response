#!/usr/bin/bash
# IgCaller analysis for WGS data
SAMPLES="GE0310 GE0319 GE0335 GE1312 GE1329 GE2328 GE3309 GE3317 GE3320 GE5312 GE5338 GE5340 GE6339 GE6345 GE6347 GE7301 GE7304 GE7305 GE7318 GE7324 GE8315 GE8321 GE8322 GE8329 GE9321 GE9336 GE9342 GE9344 GE9347 GE9348 GEK20315 GET19344 GET20301"

# enabling samtools
module load samtools

# running the IgCaller analysis for each patient sample
# utilizing hg19 reference genome, tumor normal data and the blast percentage estimates (-p parameter, replace by correct blast% estimate)
for s in ${SAMPLES}; do python3 IgCaller/IgCaller.py -I IgCaller/IgCaller_reference_files/ -V hg19 -C ensembl -T bam_files/${s}.tumor.bam -N bam_files/${s}.normal.bam -R reference_files/genome.fa -seq wgs -p 0.8; done
