# methylation data preprocessing with methylumi

# 1. Apply methylumi preprocessing to matched case-control samples,
# 2. save QC plots and
# 3. save each stage both as methylumi file AND converted to minfi-style RGChannelSet.

library('methylumi')
library('minfi')
library('ggplot2')

PIC.W <- 1920
PIC.H <- 1080

PREPROC.DIR <-'preprocessed/lumi' 
PICS.DIR <- sprintf( '%s/pics', PREPROC.DIR )
CONVERTED.DIR <- sprintf( '%s/as-rgchannelset', PREPROC.DIR )

stopifnot( dir.exists( PREPROC.DIR ) )
stopifnot( dir.exists( PICS.DIR ) )
stopifnot( dir.exists( CONVERTED.DIR ) )

my.png <- function( basename ) {
	png( sprintf( "%s/%s.png", PICS.DIR, basename ), width=PIC.W, height=PIC.H, type="cairo", antialias="subpixel")
}

# Docs refer to "barcode" when they really just mean filename prefix (that might happen to contain Illumina barcodes).
#manifest <- read.table('manifest.tab',header=F, col.names=c('id','type','barcode'))
#nor	105	0	GSM1204169_5975801020_R03C01
#nor	106	0	GSM1204170_5975801020_R05C01
manifest <- read.table('manifest.tab',header=F,
	col.names=c('state','id','qual','barcode'),
	colClasses=c('factor','character','factor','character') )

# Docs erroneasouly state methylumIDAT expects a column called 'barcodes' (plural).
data.raw <- methylumIDAT( pdat=manifest, idatPath='data')

# Following is necessary apparently for side effects; actual prefix name assigned does not matter.
sampleNames(data.raw) <- with( manifest, paste(state,id,qual,sep=".") )

# PROCESSING ###############################################################

data.bgcorr <- methylumi.bgcorr(data.raw)
save( data.bgcorr, file=sprintf( "%s/bgcorr.RData", PREPROC.DIR ) )

data.normed <- normalizeMethyLumiSet(data.bgcorr)
save( data.normed, file=sprintf( "%s/normed.RData", PREPROC.DIR ) )

# methylumi-provided QC plots
my.png('raw')
qc.probe.plot( data.raw,    by.type=T )
dev.off()

my.png('bgcorr')
qc.probe.plot( data.bgcorr, by.type=T )
dev.off()

my.png('normed')
qc.probe.plot( data.normed, by.type=T )
dev.off()

# Turn each of preceding stages into minfi-style RGChannelSet, too, and save.
# Following is optional.

#data.stripped <- stripOOB( data.normed )
# data.mapped <- mapToGenome( data.preprocessed )

rgs.data.raw    <- as(data.raw,   "RGChannelSet")
rgs.data.bgcorr <- as(data.bgcorr,"RGChannelSet")
rgs.data.normed <- as(data.normed,"RGChannelSet")

save( rgs.data.raw,   file=sprintf( "%s/raw.RData",    CONVERTED.DIR ) )
save( rgs.data.bgcorr,file=sprintf( "%s/bgcorr.RData", CONVERTED.DIR ) )
save( rgs.data.normed,file=sprintf( "%s/normed.RData", CONVERTED.DIR ) )

