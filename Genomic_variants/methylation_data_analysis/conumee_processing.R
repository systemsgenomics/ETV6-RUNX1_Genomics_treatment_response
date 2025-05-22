# processing methylation data with conumee

library('conumee')
library('methylumi')

OUT.DIR <- 'conumee-analysis'
PIC.DIR <- sprintf( '%s/plots', OUT.DIR )

stopifnot( dir.exists(PIC.DIR) )

load('preprocessed/lumi/normed.RData')

SAMPLE.NAMES <- sampleNames(data.normed)

# Docs refer to "barcode" when they really just mean filename prefix (that might happen to contain Illumina barcodes).
#nor	105	0	GSM1204169_5975801020_R03C01
#nor	106	0	GSM1204170_5975801020_R05C01
manifest <- read.table('manifest.tab',header=F,
	col.names=c('state','id','qual','barcode'),
	colClasses=c('factor','character','factor','character') )

stopifnot( nrow(manifest) == 882 )

rois <- GRanges(
	seqnames=c("chr9","chr12"),
	IRanges( start=c(15000000,20000000), end=c(30000000,30000000) ),
	name=c('CDKN2A','KRAS') )

anno <- CNV.create_anno(
	array_type="450k",
	exclude_regions=NULL,
	detail_regions=rois)

mset <- as(data.normed,"MethylSet")
con.data <- CNV.load(mset)
# con.data is only indexable by integer index.
con.data.controls <- grep( '^nor', SAMPLE.NAMES )

my.pic <- function( id ) {
	jpeg( sprintf('%s/%s.jpg', PIC.DIR, id, type ),
		width=1920,height=1080,type="cairo",antialias="subpixel")
}

for( i in 1:length(SAMPLE.NAMES) ) {

	sample.id <- SAMPLE.NAMES[i]
	title <- sprintf( "Sample %s", sample.id )

	x <- CNV.fit( con.data[ i ], con.data[con.data.controls], anno )
	x <- CNV.bin(x)
	x <- CNV.detail(x)
	x <- CNV.segment(x)

	my.pic( sample.id )
		CNV.genomeplot(x, main=title )
	dev.off()

	save( x, file=sprintf("%s/%s.RData", OUT.DIR, sample.id ) )
	CNV.write( x, what="segments", file=sprintf("%s/%s.seg",      OUT.DIR, sample.id ) )
	CNV.write( x, what="bins",     file=sprintf("%s/%s-bins.igv", OUT.DIR, sample.id ) )
	CNV.write( x, what="detail",   file=sprintf("%s/%s-det.txt",  OUT.DIR, sample.id ) )
	CNV.write( x, what="probes",   file=sprintf("%s/%s-probes.igv", OUT.DIR, sample.id ) )
	#break; # for testing
}

