#!/bin/bash
# run CellRanger

cellranger-6.0.2//cellranger count --id=ID --transcriptome=refdata-cellranger-GRCh38-2020-A/GRCh38 --feature-ref=ab_barcodes.csv --libraries=libraries.csv --localcores=24 --localmem=200 --disable-ui

