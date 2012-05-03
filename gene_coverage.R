# Rscript for coverage statistics for a list of genes in a sample
# Rscript coverage.R '<full path to bam file>' geneA geneB geneX minCoverage

.libPaths('/home/fro061/R/library')
suppressMessages(require(Rsamtools,quiet=TRUE))
suppressMessages(require(multicore,quiet=TRUE))
suppressMessages(require(GenomicFeatures,quiet=TRUE))
suppressMessages(require(GenomicRanges,quiet=TRUE))
suppressMessages(require(biomaRt,quiet=TRUE))

VERSION <- '1.0-Apr2012'

createBamIndex <- function
### creates index for bam file if it does not exist
(bamFile
### path to bam file
) {
    baiFile = paste(bam, '.bai', sep='')
    if (!file.exists(baiFile)) {
        indexBam(bamFile)
    } 
}

getSpecificRegion <- function
### function to return a GRanges object for specific region of a bam file
### should reduce memory consumption because the bam file is not read entirely
(chr, 
### chromosome
chrStart,
### start position (bp)
chrEnd,
### end position (bp)
bamFile
### bam file
) {
    param <- ScanBamParam(what = c("rname", "strand","pos", "qwidth"),
                        which = GRanges(chr,IRanges(chrStart, chrEnd)),
                        flag = scanBamFlag(isUnmappedQuery = FALSE)
                        )                    

    x <- scanBam(bamFile, param = param)[[1]]
    ranges = GRanges(seqnames=Rle(x$rname), ranges=IRanges(x$pos, width=x$qwidth))
    seqnames(ranges) <- sub("^(\\d+)","chr\\1",seqnames(ranges))
    ranges
}


args <- commandArgs(trailingOnly = TRUE)
bam <- args[1]
genes <- args[2:length(args)]
# genes <- c(args[2:(length(args)-1)])
# minCoverage <- as.integer(tail(args,1))

# create index file if it does not exist
createBamIndex(bam)

#load ensembl mart
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes_info <- getBM(attributes = c("hgnc_symbol","chromosome_name","start_position","end_position"), filters = "hgnc_symbol", values=genes, mart=human)

#convert genes to ccds ids
#genes_to_ccds[1] gene ids genes_to_ccds[2] ccds ids
# genes_to_ccds <- synergizer(authority="ensembl", species="Homo sapiens", domain="hgnc_symbol", range="ccds", ids=genes)

# txdb <- makeTranscriptDbFromUCSC(genome='hg19',tablename='ccdsGene')
# saveFeatures(txdb,'/export/astrakanfs/mpesj/reference/ccdsGene.hg19.mar2012.sqlite')
# stop('finished')
FEATURES <- '/export/astrakanfs/mpesj/reference/ccdsGene.hg19.mar2012.sqlite'
txdb = loadFeatures(FEATURES)

cat('Script version', VERSION,'\n')
cat('Database', FEATURES,'\n\n')
cat(bam,'\n')
cat('Candidate Genes: ')
cat(genes)
cat('\n\n')
cat("Individual gene coverage analysis")
cat('\n')

#iterate over each gene
for (i in 1:length(genes)) {
	info <- subset(genes_info, hgnc_symbol == genes[i])
	if (dim(info)[1] > 0) {
		#focus on the gene region
		chr <- info$chromosome_name
		start <- info$start_position
		end <- info$end_position
		bamRegion <- getSpecificRegion(chr,start,end,bam)
		region <- GRanges(chr,IRanges(start, end))
		seqnames(region) <- sub("^(\\d+)","chr\\1",seqnames(region))

		hg19.transcripts <- transcriptsByOverlaps(txdb,region)
		hg19.exons <- exonsByOverlaps(txdb,region)

		size.transcripts <- sum(width(hg19.transcripts))
		number.transcripts <- length(hg19.transcripts)
		number.exons <- length(hg19.exons)

		#intersect the transcript range with the actual reads reported, to calculate coverage
		coverage.transcripts <- Views(coverage(bamRegion)[paste('chr',chr,sep='')],as(hg19.transcripts,"RangesList")[paste('chr',chr,sep='')])
		coverage.exons <- Views(coverage(bamRegion)[paste('chr',chr,sep='')],as(hg19.exons,"RangesList")[paste('chr',chr,sep='')])

		#The mean is the weighted mean
		coverage.exon.means <- viewMeans(coverage.exons)[[1]]%*%width(coverage.exons)[[1]]/sum(width(coverage.exons)[[1]])

		# add coverage information to the GRanges object as metadata
		elementMetadata(hg19.transcripts)$mean_coverage <- as.vector(viewMeans(coverage.transcripts))
		elementMetadata(hg19.exons)$mean_coverage <- as.vector(viewMeans(coverage.exons))
		
		cat(genes[i]," - ",number.transcripts," transcripts, ",number.exons," exons, ", formatC(coverage.exon.means,digits=2,format='f'),'X mean exonic coverage',sep='')
		cat("\n")
		cat("List of exons:")
		cat("\n")
		cat(paste('chromosome','start','end', 'mean_coverage\n',sep='\t'))
		for (j in 1:number.exons) {
			meanCoverage <- formatC(elementMetadata(hg19.exons[j])$mean_coverage,digits=2,format="f")
      chromstart <- as.data.frame(hg19.exons)[j,][2]
      chromend <- as.data.frame(hg19.exons)[j,][3]
			cat(paste(chr,chromstart,chromend,meanCoverage,'\n', sep="\t"))
		}
		
		cat('\n')
	} else {
			cat(paste("Gene symbol",genes[i],"not found."))
			cat('\n')
	}
}
