#!/bin/bash
# RNAseq data analysis with nf-core/rnaseq

# enabling R
module load r/3.6.1

# running analysis
# using hg38 reference geneome
NXF_VER="21.02.0-edge" nextflow run nf-core/rnaseq -r 3.0 \
        --input sampleSheet_C101HW18062568-1.csv \
        --fasta GRCh38.103/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        --gtf GRCh38.103/Homo_sapiens.GRCh38.103.gtf \
        --star_index genome/index/star \
        --max_cpus 16 \
        --max_memory '70.GB' \
        --save_align_intermeds \
        --save_unaligned \
        --outdir data/RNAseq \
        -profile singularity
