# Rscript for coverage statistics and number of features in a region of a bam file
# Rscript coverage.R '/export/astrakanfs/mpesj/Agilent/1040PRN0046_GRCh37.gatk.bam' 7 36399490 55889334 8

.libPaths('/home/fro061/R/library')
suppressMessages(require(Rsamtools,quiet=TRUE))
suppressMessages(require(multicore,quiet=TRUE))
suppressMessages(require(GenomicFeatures,quiet=TRUE))
suppressMessages(require(GenomicRanges,quiet=TRUE))
suppressMessages(require(SynergizeR,quiet=TRUE))

VERSION <- '1.1-Mar2012'

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
    # coverage(ranges)
    #coverage(IRanges(x[["pos"]], width = x[["qwidth"]]))
    ### coverage IRLe object
}


args <- commandArgs(trailingOnly = TRUE)
bam <- args[1]
chr <- args[2]
start <- as.integer(args[3])
end <- as.integer(args[4])
minCoverage <- as.integer(args[5])

# ccds <- parseBedFile('ccds.bed')
# bam <- '/export/astrakanfs/mpesj/Agilent/1040PRN0046_GRCh37.gatk.bam'
# chr <- 7
# start<-36399490
# end<-55889334
# minCoverage<-8
createBamIndex(bam)
bamRegion <- getSpecificRegion(chr,start,end,bam)

# txdb <- makeTranscriptDbFromUCSC(genome='hg19',tablename='ccdsGene')
# saveFeatures(txdb,'/export/astrakanfs/mpesj/reference/ccdsGene.hg19.mar2012.sqlite')
# stop('finished')
FEATURES <- '/export/astrakanfs/mpesj/reference/ccdsGene.hg19.mar2012.sqlite'
txdb = loadFeatures(FEATURES)

#Calulate transcript overlaps with the full genomic range
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

transcripts.less_than_minimum.covered <- length(which(elementMetadata(hg19.transcripts)[, "mean_coverage"] < minCoverage))
exons.less_than_minimum <- hg19.exons[which(elementMetadata(hg19.exons)[, "mean_coverage"] < minCoverage)]
exons.less_than_minimum.covered <- length(exons.less_than_minimum)

transcripts.pct = round(transcripts.less_than_minimum.covered/number.transcripts*100)
exons.pct = round(exons.less_than_minimum.covered/number.exons*100)

#exon to transcripts mapping
exons_to_transcripts <- as.matrix(findOverlaps(exons.less_than_minimum,hg19.transcripts))
exons_to_transcripts[,2] <- elementMetadata(hg19.transcripts[exons_to_transcripts[, 2]])[, "tx_name"]
#which(exons_to_transcripts[,1] == '181')
#translate to hgnc_symbols from ccds transcript id's
transcripts_to_genes <- data.frame(transcript=elementMetadata(hg19.transcripts)[, "tx_name"], hgnc=synergizer(authority="ensembl", species="Homo sapiens", domain="ccds", range="hgnc_symbol", ids=sub("\\..*",'',as.vector(elementMetadata(hg19.transcripts)[, "tx_name"])))[,2])
number.genes<-length(unique(transcripts_to_genes$hgnc))

cat('Script version', VERSION,'\n')
cat('Database', FEATURES,'\n\n')
cat(bam,'\n')
cat(paste('Candidate Region: ',chr,':',start,'-',end,sep=''))
cat('\n\n')
cat('Total length of Exons in candidate region overlapping capture enrichment:\n')
cat(paste(size.transcripts,'bp, ',number.genes,' genes, ',number.transcripts,' transcripts, ',number.exons,' exons, mean coverage: ',formatC(coverage.exon.means,digits=2,format='f'),'X\n',sep=''))
cat('List of genes in region:\n')
cat(as.vector(unique(transcripts_to_genes$hgnc[transcripts_to_genes$hgnc != "NA"])))
cat('\n')
# cat('Number and % of transcripts covered at mean < 8X:\n')
# cat(paste('n=',transcripts.less_than_minimum.covered,', ',transcripts.pct,'% of target.\n',sep=''))
cat(paste('Number and % of exons covered at mean < ',minCoverage,'X:\n',sep=''))
cat(paste('n=',exons.less_than_minimum.covered,', ',exons.pct,'% of target.\n',sep=''))

if (exons.less_than_minimum.covered > 0) {
    cat('\nList of poor coverage exons in region:\n')
    cat(paste('chromosome','start','end', 'mean_coverage','hgnc_symbol','transcript_id\n',sep='\t'))

    results <- as.data.frame(exons.less_than_minimum)
    for(i in 1:exons.less_than_minimum.covered) {
        meanCoverage <- formatC(elementMetadata(exons.less_than_minimum[i])[1,"mean_coverage"],digits=2)
        transcripts <- exons_to_transcripts[which(exons_to_transcripts[,1] == i),2]
        genes <- c()
        for (transcript in exons_to_transcripts[which(exons_to_transcripts[,1] == i),2]) {
            genes <- append(genes, as.vector(transcripts_to_genes$hgnc[transcripts_to_genes$transcript==transcript]))
        }
        genes <- unique(na.omit(genes[genes != "NA"]))
        # chromosome <- results[i,][1]
		chromosome <- chr
        chromstart <- results[i,][2]
        chromend <- results[i,][3]

        cat(paste(chromosome,chromstart,chromend,meanCoverage,'', sep="\t"))
        cat(genes)
        cat('\t')
        cat(transcripts)
        cat('\n')
    }
}


