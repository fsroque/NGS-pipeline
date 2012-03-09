#number of cpus to use when multithreading
N_CPUS = 1
#path to binaries
PICARD = '/usr/local/esysbio/src/picard-tools'
GATK = '/usr/local/esysbio/src/GenomeAnalysisTK/dist/GenomeAnalysisTK.jar'
#reference files
REFERENCE = '/export/astrakanfs/fro061/reference/human_g1k_v37.clean.fasta'
DBSNP = '/export/astrakanfs/fro061/reference/dbsnp_132_GRCh37.vcf'
#EXOME = '/export/astrakanfs/fro061/reference/Agilent_SureSelect_v1_GRCh37.bed'
EXOME = '/export/astrakanfs/fro061/reference/Nimblegen_SeqCap_EZ_Exome_v2_37_targetRegOnly_wingspan.bed'
#should be twice the average read depth
#This prevents SNP calling in regions with extreme depth, such as centromeres.
# setting to 1000 to disable it
VARFILTER = 1500

TARGETS := $(wildcard *.gatk.bam)

indels: $(patsubst %.gatk.bam,%.indels,$(TARGETS))
snps: $(patsubst %.gatk.bam,%.snps,$(TARGETS))

.PHONY: indels snps
.PRECIOUS: %.gatk.bam %.bcf

%.bcf: %.gatk.bam
	@samtools mpileup -ugf $(REFERENCE) $< | bcftools \
	view -bvcg - > $@

%.snp.vcf: %.bcf
	@bcftools view -I $< | vcfutils.pl varFilter -D $(VARFILTER) > $@

%.snp.recode.vcf: %.snp.vcf
	@vcftools --vcf $< --out $*.snp --recode --bed $(EXOME) --keep-INFO-all \
	>/dev/null
	rm $*.snp.log

%.indel.vcf: %.bcf
	@bcftools view $< | vcfutils.pl varFilter -D $(VARFILTER) | \
	egrep '#|INDEL' > $@

%.indels: %.indel.vcf
	@vcftools --vcf $< --out $@ --recode --bed $(EXOME) --keep-INFO-all \
	>/dev/null
	mv $@.recode.vcf $@
	rm $@.log 

%.snps: %.snp.recode.vcf %.gatk.bam
	@java -Xmx4g -jar $(GATK) -l INFO -T VariantAnnotator \
	--assume_single_sample_reads sample \
	-R $(REFERENCE) --dbsnp $(DBSNP) -I $*.gatk.bam --variant $< \
	-o $@ -A AlleleBalance -A MappingQualityZero -A LowMQ \
	-A RMSMappingQuality -A HaplotypeScore -A QualByDepth \
	-A DepthOfCoverage -A HomopolymerRun > /dev/null
	rm $@.idx $<.idx
