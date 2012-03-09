createBamIndex <- function
### creates index for bam file if it does not exist
(bamFile
### path to bam file
) {
    library(Rsamtools)
    baiFile = paste(bam, '.bai', sep='')
    if (!file.exists(baiFile)) {
        print('Bam index file does not exist, creating.')
        indexBam(bamFile)
    } else{
        print('Bam index file exists, skipping.')
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
    library(Rsamtools)
    library(multicore)
    # param <- ScanBamParam(what = c("pos", "qwidth"),
    #                     which = GRanges(chr, IRanges(chrStart, chrEnd)),
    #                     flag = scanBamFlag(isUnmappedQuery = FALSE)
    #                     )
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

# Rscript coverage.R '/export/astrakanfs/mpesj/36-46/1040PRN0046_GRCh37.gatk.bam' 7 36399490 55889334 8
args <- commandArgs(trailingOnly = TRUE)
bam <- args[1]
chr <- args[2]
start <- as.integer(args[3])
end <- as.integer(args[4])
minCoverage <- as.integer(args[5])

# ccds <- parseBedFile('ccds.bed')
# bam <- '/export/astrakanfs/mpesj/36-46/1040PRN0046_GRCh37.gatk.bam'
createBamIndex(bam)
bamRegion <- getSpecificRegion(chr,start,end,bam)

library(GenomicFeatures)
library(GenomicRanges)

# txdb <- makeTranscriptDbFromUCSC(genome='hg19',tablename='ccdsGene')
# saveFeatures(txdb,'ccdsGene.hg19.sqlite')
# txdb <- makeTranscriptDbFromUCSC(genome='hg19',tablename='knownGene')
# saveFeatures(txdb,'knownGene.hg19.sqlite')
# txdb <- makeTranscriptDbFromUCSC(genome='hg19',tablename='ensGene')
# saveFeatures(txdb,'ensGene.hg19.sqlite')
txdb = loadFeatures('ccdsGene.hg19.sqlite')
# txdb = loadFeatures('knownGene.hg19.sqlite')

hg19.transcripts <- transcriptsByOverlaps(txdb,bamRegion)
hg19.exons <- exonsByOverlaps(txdb,bamRegion)

size.transcripts <- sum(width(hg19.transcripts))
number.transcripts <- length(hg19.transcripts)
number.exons <- length(hg19.exons)
# coverage.region <- mean(coverage(bamRegion)[paste('chr',chr,sep='')])



# coverage <- coverage(bamRegion)
# results <- slice(coverage, upper=8)
# results[paste('chr',chr,sep='')] <- restrict(results[paste('chr',chr,sep='')],start,end)
# subsetByOverlaps(transcripts(txdb), results)
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

cat(paste('Candidate Region: ',chr,':',start,'-',end,sep=''))
cat('\n\n')
cat('CCDS transcripts total length in candidate region:\n')
cat(paste(size.transcripts,'bp, ',number.transcripts,' transcripts, ',number.exons,' exons, mean coverage in ccds region: ',formatC(coverage.exon.means,digits=2,format='f'),'X\n',sep=''))
cat('Number and % of transcripts covered at mean < 8X:\n')
cat(paste('n=',transcripts.less_than_minimum.covered,', ',transcripts.pct,'% of target.\n',sep=''))
cat('Number and % of exons covered at mean < 8X:\n')
cat(paste('n=',exons.less_than_minimum.covered,', ',exons.pct,'% of target.\n',sep=''))

if (exons.less_than_minimum.covered > 0) {
    cat('\nList of poor coverage exons in region:\n')
    cat(paste('chromosome','start','end', 'mean_coverage','transcript_id\n',sep=' '))

    results <- as.data.frame(exons.less_than_minimum)
    for(i in 1:exons.less_than_minimum.covered) {
        meanCoverage <- formatC(elementMetadata(exons.less_than_minimum[i])[1,"mean_coverage"],digits=2)
        transcripts <- exons_to_transcripts[which(exons_to_transcripts[,1] == i),2] 
        chromosome <- results[i,][1]
        chromstart <- results[i,][2]
        chromend <- results[i,][3]
    
        cat(paste(chromosome,chromstart,chromend,meanCoverage,'', sep=" "))
        cat(transcripts)
        cat('\n')
    }
}
