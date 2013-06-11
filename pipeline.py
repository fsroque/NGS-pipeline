#!/usr/bin/env python
"""

    pipeline.py
			[--bam GLOB]
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
#   Just uncomment the desired option, and comment the ones that are not needed


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
n_cpus = 2
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

bed = os.path.basename(exome)
sys.stderr.write("WARNING: using %s as the bed file\n" % bed)

#should be twice the average read depth
#This prevents SNP calling in regions with extreme depth, such as centromeres.
# setting to 1000 to disable it
varfilter = 1500


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    from optparse import OptionParser
    import StringIO

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --bam ALIGNMENT_BAM [more_options]")
    parser.add_option("-b", "--bam", dest="bam_file",
                        metavar="FILE",
                        type="string",
                        help="Path of alignment file in BAM format. Can be specified as a glob (*.bam) to run with multiple files.")
    



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
    mandatory_options = ['bam_file']

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

def call_variants(input, output):
    """Perform multi-sample variant calling using GATK"""
    cmd = "nice java -Djava.io.tmpdir=/export/astrakanfs/stefanj/tmp -Xmx8g -jar %s \
                 -T HaplotypeCaller \
                 -R %s \
                 -I %s \
                 -stand_emit_conf 10.0 \
                 -L %s \
                 -o %s \
                 --dbsnp %s " % (gatk, reference, input, capture, output, dbsnp)
    #log the results
    cmd = cmd + '&> {}.log'.format(output)
    run_cmd(cmd)


def indel_filter(vcf, output):
    """Apply filters to the vcf file to generate the final call of indels"""
    run_cmd("vcftools --vcf %s \
             --out %s \
             --recode \
             --bed %s \
	     --keep-INFO-all "
            % (vcf, output, exome))
    run_cmd("mv %s.recode.vcf %s" %(output, output))

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

def filter_snps(input, output):
    """docstring for filter_snps"""
    run_cmd('java -Xmx2g -jar {gatk} \
       -R {reference} \
       -T VariantFiltration \
       -o {output} \
       --variant {input} \
       --filterExpression "QD < 2.0" \
       --filterExpression "MQ < 40.0" \
       --filterExpression "FS > 60.0" \
       --filterExpression "HaplotypeScore > 13.0" \
       --filterExpression "MQRankSum < -12.5" \
       --filterExpression "ReadPosRankSum < -8.0" \
       --filterName "QDFilter" \
       --filterName "MQFilter" \
       --filterName "FSFilter" \
       --filterName "HaplotypeScoreFilter" \
       --filterName "MQRankSumFilter" \
       --filterName "ReadPosRankSumFilter"'.format(
           gatk=gatk,
           reference=reference,
           output=output,
           input=input
           )
    )
           
def filter_indels(vcf, output):
    """filter indels vcf"""
    run_cmd('java -Xmx4g -jar {gatk} \
            -T VariantFiltration \
            -o {output} \
            --variant {input} \
            --filterExpression "QD < 2.0" \
            --filterExpression "ReadPosRankSum < -20.0"   \
            --filterExpression "FS > 200.0"   \
            --filterName QDFilter   \
            --filterName ReadPosFilter   \
            --filterName FSFilter \
            -R {reference}'.format(
                gatk=gatk,
                output=output,
                input=vcf,
                reference=reference
            ))

def remove_filtered(vcf,output):
    """Remove filtered variants"""
    run_cmd("java -Xmx2g -jar {gatk} \
            -T SelectVariants \
            -R {reference} \
            --variant {input} \
            -o {output} \
            -env -ef".format(
                gatk=gatk,
                reference=reference,
                input=vcf,
                output=output
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
            
def cleanup_files():
    run_cmd("rm -rf */*.recal_data.csv */*.realign* */*.dedup* */*.log *.log \
            ")

def shared_snps_btw_samples(inpufile_list,pathdir):
    '''
    Creates intersections and complements of two or more VCF files. 
    Given multiple VCF files, it can output the list of positions which are shared by at least N files, at most N files, exactly N files, etc. 
    It receives a list of input files
    (Not used)
    '''

    outfile_name=pathdir+'.vcf.isec.gz'

    args = ['vcf-isec', '-n', '=', len(inpufile_list)]+ inpufile_list + ['|','bgzip','-c']
    f = open(outfile_name, "w")
    retcode = subprocess.call(args, stdout=f)
    f.close()

    return retcode

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
    files = glob.glob(options.bam_file)
    parameters = []
    for file in files:
        prefix = os.path.splitext(os.path.basename(file))[0]
        parameters.append([None,prefix + '/' + prefix + '.bam', [prefix,os.path.abspath(file)]])
    for job_parameters in parameters:
            yield job_parameters

@files(generate_parameters)
def link(none, bam, extra):
    """Make working directory, chdir and make symlink to bam file"""
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
    
    
@follows('reduce_bam')
@transform(reduce_bam, suffix('.reduced.bam'), '.vcf')
def call_variants_hc(bam, vcf):
    """Call snps on the data, using samtools"""
    call_variants(bam, vcf)    

@follows('call_variants_hc')
@transform(call_variants_hc, suffix('.vcf'), '.indel.vcf')
def vcf_to_indels(vcf, indel):
    """Converts bcf files to vcf, only indels"""
    split_indels(vcf, indel)
    
@follows('vcf_to_indels')
@transform(vcf_to_indels, suffix('.indel.vcf'), '.indel.filtered.vcf')
def filter_indel_file(vcf, output):
    """Converts bcf files to vcf, only indels"""
    filter_indels(vcf, output)

@follows('call_variants_hc')
@transform(call_variants_hc, suffix('.vcf'), '.snp.vcf')
def vcf_to_snps(vcf, snp):
    """Converts bcf files to vcf, only indels"""
    split_snps(vcf, snp)

@follows('vcf_to_snps')
@transform(vcf_to_snps, suffix('.snp.vcf'), '.snp.filtered.vcf')
def filter_snp_file(vcf, output):
    """Converts bcf files to vcf, only indels"""
    filter_snps(vcf, output)

@posttask(cleanup_files)
@merge([filter_snp_file, filter_indel_file],'gatk.variants.vcf')
def merge_variants(filtered,variants):
    """Merge snp and indel files"""
    merge_indel_and_snp_vcf(filtered[0],filtered[1],variants)


# @follows('filter_snps')
# @transform(filter_snps, suffix('.snp.recode.vcf'), add_inputs(r'\1.gatk.bam'), '.snps')
# def annotate_vcf(inputs, annotated_vcf):
#     """Annotates snp vcf file"""
#     annotate_vcf_gatk(inputs[0], inputs[1], annotated_vcf)

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
                            verbose         = options.verbose)

