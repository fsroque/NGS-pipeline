#!/usr/bin/env python
"""

    pipeline_multisample.py
			[--bamdir PATH]
			[--groups INT]
                        [--log_file PATH]
                        [--verbose]
                        [--target_tasks]
                        [--jobs]
                        [--just_print]
                        [--flowchart]
                        [--key_legend_in_graph]
                        [--forced_tasks]

"""
import sys
import os
import glob

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   user definable options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#path to binaries
script_path = os.path.dirname(os.path.realpath(__file__))
picard = os.path.join(script_path,'../src/picard-tools')
gatk = os.path.join(script_path,'../src/GenomeAnalysisTK/GenomeAnalysisTK.jar')
#reference files
reference = os.path.join(script_path,'../reference/human_g1k_v37.clean.fasta')
dbsnp = os.path.join(script_path,'../reference/dbsnp_137.b37.vcf')
exome = os.path.join(script_path,'../reference/Nimblegen_SeqCap_EZ_Exome_v2_37_targetRegOnly_wingspan.bed')
capture = os.path.join(script_path,'../reference/Nimblegen_SeqCap_EZ_Exome_v2_37_targetRegOnly_g1k.bed')
hapmap = os.path.join(script_path,'../reference/hapmap_3.3.b37.sites.vcf')
omni = os.path.join(script_path,'../reference/1000G_omni2.5.b37.sites.vcf')
mills = os.path.join(script_path,'../reference/Mills_and_1000G_gold_standard.indels.b37.vcf')

n_cpus = 2


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    from optparse import OptionParser
    import StringIO

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --bamdir BAM_DIR --groups NUMBER [more_options]")
    parser.add_option("-b", "--bamdir", dest="bam_dir",
                        metavar="FILE",
                        type="string",
                        help="Directory containing all the bams for analysis.")
                        
    parser.add_option("-g", "--groups", dest="groups",
                        type="int",
                        default=1,
                        help="Split dataset into smaller groups for multisample snp caling.")



    #
    #   general options: verbosity / logging
    #
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="count", default=0,
                      help="Print more verbose messages for each additional verbose level.")
    parser.add_option("-L", "--log_file", dest="log_file",
                      metavar="FILE",
                      type="string",
                      help="Name and path of log file")




    #
    #   pipeline
    #
    parser.add_option("-t", "--target_tasks", dest="target_tasks",
                        action="append",
                        default = list(),
                        metavar="JOBNAME",
                        type="string",
                        help="Target task(s) of pipeline.")
    parser.add_option("-j", "--jobs", dest="jobs",
                        default=1,
                        metavar="N",
                        type="int",
                        help="Allow N jobs (commands) to run simultaneously.")
    parser.add_option("-n", "--just_print", dest="just_print",
                        action="store_true", default=False,
                        help="Don't actually run any commands; just print the pipeline.")
    parser.add_option("--flowchart", dest="flowchart",
                        metavar="FILE",
                        type="string",
                        help="Don't actually run any commands; just print the pipeline "
                             "as a flowchart.")

    #
    #   Less common pipeline options
    #
    parser.add_option("--key_legend_in_graph", dest="key_legend_in_graph",
                        action="store_true", default=False,
                        help="Print out legend and key for dependency graph.")
    parser.add_option("--forced_tasks", dest="forced_tasks",
                        action="append",
                        default = list(),
                        metavar="JOBNAME",
                        type="string",
                        help="Pipeline task(s) which will be included even if they are up to date.")
    parser.add_option("--rebuild_mode", dest="rebuild_mode",
                        action="store_false", default=True,
                        help="gnu_make_maximal_rebuild_mode")

    # get help string
    f =StringIO.StringIO()
    parser.print_help(f)
    helpstr = f.getvalue()
    (options, remaining_args) = parser.parse_args()


    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    #
    #   Add names of mandatory options,
    #       strings corresponding to the "dest" parameter
    #       in the options defined above
    #
    mandatory_options = ['bam_dir']

    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    def check_mandatory_options (options, mandatory_options, helpstr):
        """
        Check if specified mandatory options have b een defined
        """
        missing_options = []
        for o in mandatory_options:
            if not getattr(options, o):
                missing_options.append("--" + o)
    
        if not len(missing_options):
            return
    
        raise Exception("Missing mandatory parameter%s: %s.\n\n%s\n\n" %
                        ("s" if len(missing_options) > 1 else "",
                         ", ".join(missing_options),
                         helpstr))
    check_mandatory_options (options, mandatory_options, helpstr)


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   imports


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from ruffus import *
import subprocess
# import drmaa
import resource

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Functions


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
def setlimits():
    """Set maximum meomory to be used by a child process"""
    resource.setrlimit(resource.RLIMIT_AS, (100000000000,100000000000))

def run_cmd(cmd_str):
    """
    Throw exception if run command fails
    """
    process = subprocess.Popen(cmd_str, stdout = subprocess.PIPE,
                                stderr = subprocess.PIPE, shell = True, preexec_fn=setlimits)
    stdout_str, stderr_str = process.communicate()
    if process.returncode != 0:
        raise Exception("Failed to run '%s'\n%s%sNon-zero exit status %s" %
                            (cmd_str, stdout_str, stderr_str, process.returncode))

def rename(file,new_file):
    """rename file"""
    os.rename(file,new_file)

def remove(file):
    """remove file"""
    os.remove(file)
                            
def split_seq(seq, num_pieces):
    """ split a list into pieces passed as param """
    start = 0
    for i in xrange(num_pieces):
        stop = start + len(seq[i::num_pieces])
        yield seq[start:stop]
        start = stop

def index_bam(bam):
    """Use samtools to create an index for the bam file"""
    run_cmd('samtools index %s' % bam)
    
def find_realigns(bam, intervals_file):
    """Finds regions needing realignment"""
    run_cmd("java -Xmx4g -jar %s \
             -T RealignerTargetCreator \
             -I %s \
             -R %s \
             -o %s \
             -nt %d" 
             % (gatk, bam, reference, intervals_file, n_cpus))
             
def run_realigner(bam, intervals_file, realigned_bam):
    """Performs local realignment of reads based on misalignments due to the presence of indels"""
    run_cmd("java -Xmx4g -jar %s \
             -T IndelRealigner \
             -I %s \
             -R %s \
             -targetIntervals %s \
             -o %s" 
             % (gatk, bam, reference, intervals_file, realigned_bam))

# def dup_removal_picard(bam,output):
#     """Use Picard to remove duplicates"""
#     run_cmd('java -Xmx4096m -jar %s/CleanSam.jar INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR' 
#                 % (picard,bam,output))
# 
# def dup_removal_samtools(bam,output):
#     """Use samtools for dupremoval"""
#     run_cmd('samtools rmdup %s %s' % (bam,output))
    
def dup_mark_picard(bam,output):
    """Use Picard to mark duplicates"""
    run_cmd('java -Xmx4096m -jar %s/MarkDuplicates.jar TMP_DIR=/export/astrakanfs/stefanj/tmp REMOVE_DUPLICATES=true INPUT=%s OUTPUT=%s METRICS_FILE=%s.dup_metrics VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR CREATE_INDEX=true'
                % (picard,bam,output,bam))

def base_recalibrator(bam, recal_data):
    """First pass of the recalibration step"""
    run_cmd("java -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx8g -jar %s \
            -T BaseRecalibrator \
            -R %s \
            -knownSites %s \
            -I %s \
            -o %s \
            -nct %d"
            % (gatk, reference, dbsnp, bam, recal_data, n_cpus))

def print_recalibrated(bam, recal_data, output):
    """uses GATK to rewrite quality scores using the recal_data"""
    run_cmd("java -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx4g -jar %s \
            -T PrintReads \
            -R %s \
            -I %s \
            --out %s \
            -BQSR %s \
            -nct %d" 
            % (gatk, reference, bam, output, recal_data, n_cpus) )
            
# def count_covariates(bam, recal_data):
#     """Uses GATK to count covariates"""
#     run_cmd("java -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx4g -jar %s \
#             -T CountCovariates \
#             -l INFO \
#             -R %s \
#             -knownSites %s \
#             --default_platform illumina \
#             -cov ReadGroupCovariate \
#             -cov QualityScoreCovariate \
#             -cov CycleCovariate \
#             -cov DinucCovariate \
#             -I %s \
#             -recalFile %s"
#             % (gatk, reference, dbsnp, bam, recal_data))
            
# def table_recalibration(bam, recal_data, output):
#     """uses GATK to rewrite quality scores using the recal_data"""
#     run_cmd("java -Xmx4g -jar %s \
#             -T TableRecalibration \
#             --default_platform illumina \
#             -R %s \
#             --preserve_qscores_less_than 5 \
#             -l INFO \
#             -I %s \
#             --out %s \
#             -recalFile %s" 
#             % (gatk, reference, bam, output, recal_data) )

def reduce_reads(bam, output):
    """Reduces the BAM file using read based compression that keeps only essential information for variant calling"""
    run_cmd("java -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx15g -jar %s \
            -T ReduceReads \
            -R %s \
            -I %s \
            -o %s"
            % (gatk, reference, bam, output))
            

def bam_quality_score_distribution(bam,qs,pdf):
    """Calculates quality score distribution histograms"""
    run_cmd("java -jar {picard}/QualityScoreDistribution.jar \
             CHART_OUTPUT={chart} \
             OUTPUT={output} \
             INPUT={bam} \
             VALIDATION_STRINGENCY=SILENT".format(
                 picard=picard,
                 chart=pdf,
                 output=qs,
                 bam=bam
             ))

def bam_alignment_metrics(bam,metrics):
    """Collects alignment metrics for a bam file"""
    run_cmd("java -jar {picard}/CollectAlignmentSummaryMetrics.jar \
             REFERENCE_SEQUENCE={reference} \
             OUTPUT={output} \
             INPUT={bam} \
             VALIDATION_STRINGENCY=SILENT".format(
                 picard=picard,
                 reference=reference,
                 output=metrics,
                 bam=bam
             ))
             
def bam_coverage_statistics(bam, statistics):
    run_cmd("java -Xmx4g -jar {gatk} \
            -R {reference} \
            -T DepthOfCoverage \
            -o {output} \
            -I {input} \
            -L {capture} \
            -ct 8 -ct 20 -ct 30 \
            --omitDepthOutputAtEachBase --omitLocusTable \
            ".format(gatk=gatk,
                reference=reference,
                output=statistics,
                input=bam,
                capture=capture
            ))

def filter_by_exome_region(vcf, output):
    """Apply filters to the vcf file to limit calling to exome region"""
    run_cmd("vcftools --vcf %s \
             --out %s \
             --recode \
             --bed %s \
             --keep-INFO-all "
            % (vcf, output, exome))

def multisample_variant_call(files, output):
    """Perform multi-sample variant calling using GATK"""
    cmd = "nice java -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx8g -jar %s \
            -T UnifiedGenotyper \
            -R %s \
            -o %s \
            -glm BOTH \
            -nt %s \
            --dbsnp %s " % (gatk, reference, output, options.jobs, dbsnp)
    for file in files:
        cmd = cmd + '-I {} '.format(file)
    #log the results
    cmd = cmd + '&> {}.log'.format(output)
    run_cmd(cmd)
    
def merge_batch_vcf(vcfs, output):
    """Merges vcf files from the batch run"""
    if len(vcfs) == 1:
        run_cmd('cp {vcf} {output}'.format(vcf = vcfs[0], output = output))
    else:
        merging = ''
        for i in range(len(vcfs)):
            merging = merging + ' -V:batch{number} {file}'.format(number=i,file=vcfs[i])
        run_cmd("java -Xmx4g -jar {gatk} \
                -R {reference} \
                -T CombineVariants \
                -o {output} \
                {files}".format(
                    gatk=gatk,
                    reference=reference,
                    output=output,
                    files=merging
                )) 
            
def merge_indel_and_snp_vcf(snp,indel, output):
    """Merges vcf files from the batch run"""
    run_cmd("java -Xmx4g -jar {gatk} \
            -R {reference} \
            -T CombineVariants \
            -V:SNP {snp} \
            -V:INDEL {indel} \
            -o {output}".format(
                gatk=gatk,
                reference=reference,
                snp=snp,
                indel=indel,
                output=output
            )) 

def split_variants_by_sample(vcf, sample, output):
    run_cmd("java -Xmx2g -jar {} -R {} \
            -T SelectVariants \
            --variant {} -sn {} \
            -o {}".format(gatk,reference,vcf,sample,output))

def split_snps(vcf,snp_file):
    """Select for snps in the vcf file"""
    run_cmd("java -Xmx2g -jar {} -R {} -T SelectVariants \
            --variant {} -selectType SNP\
            -o {}".format(gatk,reference,vcf,snp_file))
            
def split_indels(vcf,indel_file):
    """Select for indels in the vcf file"""
    run_cmd("java -Xmx2g -jar {} -R {} -T SelectVariants \
            --variant {} -selectType INDEL\
            -o {}".format(gatk,reference,vcf,indel_file))

def variant_annotator(bam, vcf, output):
    run_cmd("java -Xmx16g -jar {gatk} \
            -R {reference} \
            -T VariantAnnotator \
            -I {bam} \
            -o {output} \
            -A Coverage \
            --variant {vcf} \
            --dbsnp {dbsnp} \
            ".format(gatk=gatk,
                bam=bam,
                reference=reference,
                output=output,
                vcf=vcf,
                dbsnp=dbsnp
            ))

def recalibrate_snps(vcf,output):
    """Runs VariantRecalibrator on the snps file"""
    cmd = "java -Xmx16g -jar {gatk} \
            -T VariantRecalibrator \
            -R {reference}  \
            -input {input} \
            --maxGaussians 6 \
            -percentBad 0.01 -minNumBad 1000 \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} \
            -resource:omni,known=false,training=true,truth=false,prior=12.0 {omni} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} \
            -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
            -mode SNP \
            -recalFile {output} \
            -tranchesFile {tranches} \
            -rscriptFile {plots}\
            -nt {num_jobs}".format(
                gatk=gatk,
                reference=reference,
                hapmap=hapmap,
                omni=omni,
                dbsnp=dbsnp,
                input=vcf,
                output=output,
                tranches=output+'.tranches',
                plots=output+'.plots.R',
                num_jobs=options.jobs
            )
    if get_num_files() > 10:
        cmd += " -an InbreedingCoeff"
    run_cmd(cmd)

def recalibrate_indels(vcf,output):
    """Runs VariantRecalibrator on the indels file"""
    cmd = "java -Xmx16g -jar {gatk} \
            -T VariantRecalibrator \
            -R {reference}  \
            -input {input} \
            --maxGaussians 4 \
            -percentBad 0.01 \
            -minNumBad 1000 \
            -resource:mills,known=false,training=true,truth=true,prior=12.0 {mills} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp} \
            -an DP -an FS -an ReadPosRankSum -an MQRankSum \
            -mode INDEL \
            -recalFile {output} \
            -tranchesFile {tranches} \
            -rscriptFile {plots}\
            -nt {num_jobs}".format(
                gatk=gatk,
                reference=reference,
                mills=mills,
                dbsnp=dbsnp,
                input=vcf,
                output=output,
                tranches=output+'.tranches',
                plots=output+'.plots.R',
                num_jobs=options.jobs
            )
    if get_num_files() > 10:
        cmd += " -an InbreedingCoeff"
    run_cmd(cmd)

def apply_recalibration_snps(vcf,recal,tranches,output):
    """Apply the recalibration tranch file calculated in recalibrate_snps"""
    run_cmd("java -Xmx16g -jar {gatk} \
            -T ApplyRecalibration \
            -R {reference} \
            -input {vcf} \
            --ts_filter_level 99.9 \
            -tranchesFile {tranches}  \
            -recalFile {recal} \
            -mode SNP \
            -nt {num_jobs} \
            -o {output}".format(
                gatk=gatk,
                reference=reference,
                vcf=vcf,
                tranches=tranches,
                recal=recal,
                output=output,
                num_jobs=options.jobs
            ))

def apply_recalibration_indels(vcf,recal,tranches,output):
    """Apply the recalibration tranch file calculated in recalibrate_snps"""
    run_cmd("java -Xmx16g -jar {gatk} \
            -T ApplyRecalibration \
            -R {reference} \
            -input:name,VCF {vcf} \
            --ts_filter_level 99.9 \
            -tranchesFile {tranches}  \
            -recalFile {recal} \
            -mode INDEL \
            -nt {num_jobs} \
            -o {output}".format(
                gatk=gatk,
                reference=reference,
                vcf=vcf,
                tranches=tranches,
                recal=recal,
                output=output,
                num_jobs=options.jobs
            ))

# def filter_indels(vcf, output):
#     """filter indels vcf"""
#     run_cmd('java -Xmx4g -jar {gatk} \
#             -T VariantFiltration \
#             -o {output} \
#             --variant {input} \
#             --filterExpression "QD < 2.0" \
#             --filterExpression "ReadPosRankSum < -20.0"   \
#             --filterExpression "InbreedingCoeff < -0.8"   \
#             --filterExpression "FS > 200.0"   \
#             --filterName QDFilter   \
#             --filterName ReadPosFilter   \
#             --filterName InbreedingFilter   \
#             --filterName FSFilter \
#             -R {reference}'.format(
#                 gatk=gatk,
#                 output=output,
#                 input=vcf,
#                 reference=reference
#             ))
# 
# def remove_filtered(vcf,output):
#     """Remove filtered variants"""
#     run_cmd("java -Xmx2g -jar {gatk} \
#             -T SelectVariants \
#             -R {reference} \
#             --variant {input} \
#             -o {output} \
#             -env -ef".format(
#                 gatk=gatk,
#                 reference=reference,
#                 input=vcf,
#                 output=output
#             ))

def cleanup_files():
    run_cmd("rm -rf */*.recal_data.csv */*.realign* */*.rmdup* */*.log *.log \
            *.to_filter* multisample.gatk.snp.recal batch* \
            multisample.gatk.recalibratedSNPS.rawIndels.vcf \
            multisample.gatk.indel.model.* multisample.gatk.snp.model.* \
            multisample.gatk.analysisReady.vcf \
            multisample.gatk.analysisReady.vcf.vcfidx \
            multisample.gatk.analysisReady.vcf.idx \
            multisample.gatk.recalibratedSNPS.rawIndels.vcf.idx \
            ")
    

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Logger


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

if __name__ == '__main__':
    import logging
    import logging.handlers

    MESSAGE = 15
    logging.addLevelName(MESSAGE, "MESSAGE")

    def setup_std_logging (logger, log_file, verbose):
        """
        set up logging using programme options
        """
        class debug_filter(logging.Filter):
            """
            Ignore INFO messages
            """
            def filter(self, record):
                return logging.INFO != record.levelno

        class NullHandler(logging.Handler):
            """
            for when there is no logging
            """
            def emit(self, record):
                pass

        # We are interesting in all messages
        logger.setLevel(logging.DEBUG)
        has_handler = False

        # log to file if that is specified
        if log_file:
            handler = logging.FileHandler(log_file, delay=False)
            handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)6s - %(message)s"))
            handler.setLevel(MESSAGE)
            logger.addHandler(handler)
            has_handler = True

        # log to stderr if verbose
        if verbose:
            stderrhandler = logging.StreamHandler(sys.stderr)
            stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
            stderrhandler.setLevel(logging.DEBUG)
            if log_file:
                stderrhandler.addFilter(debug_filter())
            logger.addHandler(stderrhandler)
            has_handler = True

        # no logging
        if not has_handler:
            logger.addHandler(NullHandler())


    #
    #   set up log
    #
    module_name = "exome"
    logger = logging.getLogger(module_name)
    setup_std_logging(logger, options.log_file, options.verbose)

    #
    #   Allow logging across Ruffus pipeline
    #
    def get_logger (logger_name, args):
        return logger

    from ruffus.proxy_logger import *
    (logger_proxy,
     logging_mutex) = make_shared_logger_and_proxy (get_logger,
                                                    module_name,
                                                    {})


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Pipeline


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#       Put pipeline code here
# bam_file = os.path.abspath(options.bam_file)
# prefix to use during the runs will be the basename of the file minus extension
# prefix = os.path.splitext(os.path.basename(bam_file))[0]

def generate_parameters():
    files = glob.glob(options.bam_dir + '/*.bam')
    parameters = []
    for file in files:
        prefix = os.path.splitext(os.path.basename(file))[0]
        parameters.append([None,prefix + '/' + prefix + '.bam', [prefix,os.path.abspath(file)]])
    for job_parameters in parameters:
            yield job_parameters

def get_num_files():
    files = glob.glob(options.bam_dir + '/*.bam')
    return len(files)

@files(generate_parameters)
def link(none, bam, extra):
    """Make working directory and make symlink to bam file"""
    if not os.path.exists(extra[0]):
        os.mkdir(extra[0])
    if not os.path.exists(bam):
        os.symlink(extra[1], bam) 

@follows(link)
@transform(link, suffix(".bam"), '.bam.bai')
def index(input, output):
    """create bam index"""
    index_bam(input)

@follows(index)
@transform(link, suffix(".bam"), '.dedup.bam')
def mark_dups(input, output):
    """Mark dups"""
    dup_mark_picard(input, output)
    # remove(input)

@follows(mark_dups)
@transform(mark_dups, suffix(".dedup.bam"), '.realign.intervals')
def find_realignment_intervals(bam,output):
   """Find regions to be re-aligned due to indels"""
   find_realigns(bam, output)

@follows(find_realignment_intervals)
@transform(find_realignment_intervals, suffix(".realign.intervals"), '.realigned.bam', r'\1.bam')
def indel_realigner(input, output, bam):
   """Re-aligns regions around indels"""
   run_realigner(bam,input,output)
   # remove(input)

@follows(indel_realigner)
@transform(indel_realigner, suffix('.realigned.bam'), '.gatk.bam.recal_data.grp')
def recalibrate_baseq1(input, output):
    """Base quality score recalibration in bam file
        Part 1: count covariates"""
    index_bam(input)
    base_recalibrator(input, output)

@follows(recalibrate_baseq1)
@transform(indel_realigner, suffix('.realigned.bam'), add_inputs(r'\1.gatk.bam.recal_data.grp'),'.gatk.bam')
def recalibrate_baseq2(inputs, output):
    """Base quality score recalibration in bam file
        Part 2: rewrite quality scores into a new bam file"""   
    print_recalibrated(inputs[0], inputs[1], output)
    # remove(inputs[1])

@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.quality_score')
def metric_quality_score_distribution(input,output):
    """docstring for metrics1"""
    bam_quality_score_distribution(input, output, output + '.pdf')

@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.metrics')
def metric_alignment(input,output):
    """docstring for metrics1"""
    bam_alignment_metrics(input, output)

@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.coverage')
def metric_coverage(input,output):
    bam_coverage_statistics(input,output)

@jobs_limit(6)
@follows(recalibrate_baseq2)
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.reduced.bam')
def reduce_bam(input, output):
    reduce_reads(input,output)

@merge(reduce_bam, 'multisample.gatk.vcf')
def call_variants(infiles, output):
    """Splits the files into s number of batches, calls the variants, and merges them back"""
    #splits the calculations into batches
    counter = 1
    batches = []
    for batch in split_seq(infiles, options.groups):
        multisample_variant_call(batch, 'batch{}.vcf'.format(counter))
        batches.append('batch{}.vcf'.format(counter))
        counter += 1
    merge_batch_vcf(batches,output)
    for batch in batches:
        remove(batch)
        remove(batch + '.idx')

# @follows('call_variants')
# @files('multisample.gatk.vcf', ['multisample.gatk.snp.vcf','multisample.gatk.indel.vcf'])
# def split_snps_and_indels(input,output):
#     """Separates the vcf file into only snp calls"""
#     split_snps(input,output[0])
#     split_indels(input,output[1])

@follows('call_variants')
@files('multisample.gatk.vcf', ['multisample.gatk.snp.model','multisample.gatk.snp.model.tranches'])
def find_snp_tranches_for_recalibration(input,output):
    """Run variantRecalibration for trusted sites"""
    recalibrate_snps(input,output[0])

@follows('call_variants')
@files('multisample.gatk.vcf', ['multisample.gatk.indel.model','multisample.gatk.indel.model.tranches'])
def find_indel_tranches_for_recalibration(input,output):
    """Run variantRecalibration for trusted sites"""
    recalibrate_indels(input,output[0])

@follows('find_snp_tranches_for_recalibration')
@files(['multisample.gatk.vcf','multisample.gatk.snp.model','multisample.gatk.snp.model.tranches'],'multisample.gatk.recalibratedSNPS.rawIndels.vcf')
def apply_recalibration_filter_snps(input,output):
    apply_recalibration_snps(input[0],input[1],input[2],output)
    # remove(input[0])
    #     remove(input[1])
    #     remove(input[2])

@follows('apply_recalibration_filter_snps')
@files(['multisample.gatk.recalibratedSNPS.rawIndels.vcf','multisample.gatk.indel.model','multisample.gatk.indel.model.tranches'],'multisample.gatk.analysisReady.vcf')
def apply_recalibration_filter_indels(input,output):
    apply_recalibration_indels(input[0],input[1],input[2],output)

# @follows('split_snps_and_indels')
# @files('multisample.gatk.indel.vcf', 'multisample.gatk.indel.to_filter.vcf')
# def apply_indel_filter(input,output):
#     """docstring for apply_indel_filter"""
#     filter_indels(input,output)
#     # remove(input)

# @follows('apply_recalibration_filter','apply_indel_filter')
# @transform(['multisample.gatk.indel.to_filter.vcf','multisample.gatk.snp.to_filter.vcf'], suffix(".to_filter.vcf"), '.filtered.vcf')
# def filter_variants(infile,outfile):
#     """Filter variants that passed"""
#     remove_filtered(infile,outfile)
#     # remove(infile)
# 
# @follows('filter_variants')
# @merge(['multisample.gatk.indel.filtered.vcf','multisample.gatk.snp.filtered.vcf'],'multisample.gatk.variants.vcf')
# def merge_variants(filtered,variants):
#     """Merge snp and indel files"""
#     merge_indel_and_snp_vcf(filtered[1],filtered[0],variants)
#     # remove(filtered[0])
#     #     remove(filtered[1])


@follows('apply_recalibration_filter_indels')
@transform(apply_recalibration_filter_indels,suffix('.analysisReady.vcf'),'.analysisReady.exome.vcf')
def final_calls(input,output):
    """Produce the final variant calls"""
    output = output[:-10]
    filter_by_exome_region(input, output)
    rename('multisample.gatk.analysisReady.recode.vcf','multisample.gatk.analysisReady.exome.vcf')


def split_snp_parameters():
    files = glob.glob(options.bam_dir + '/*.bam')
    exome_vcf = 'multisample.gatk.analysisReady.exome.vcf'
    for file in files:
        prefix = os.path.splitext(os.path.basename(file))[0]
        yield [exome_vcf, prefix + '/' + prefix + '.exome.vcf', prefix]

@posttask(cleanup_files)
@follows('final_calls')
@files(split_snp_parameters)
def split_snps(input, output, sample):
    split_variants_by_sample(input, sample, output)



#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Main logic


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if __name__ == '__main__':
    if options.just_print:
        pipeline_printout(sys.stdout, options.target_tasks, options.forced_tasks,
                            verbose=options.verbose)

    elif options.flowchart:
        pipeline_printout_graph (   open(options.flowchart, "w"),
                                    # use flowchart file name extension to decide flowchart format
                                    #   e.g. svg, jpg etc.
                                    os.path.splitext(options.flowchart)[1][1:],
                                    options.target_tasks,
                                    options.forced_tasks,
                                    no_key_legend   = not options.key_legend_in_graph)
    else:
        pipeline_run(options.target_tasks, options.forced_tasks,
                            multiprocess    = options.jobs,
                            logger          = stderr_logger,
                            verbose         = options.verbose,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode)

