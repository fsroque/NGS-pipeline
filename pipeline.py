#!/usr/bin/env python
"""

    ruffus_template.py
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
#number of cpus to use when multithreading
cpus = 1 
#path to binaries
picard = '/usr/local/esysbio/src/picard-tools'
gatk = '/usr/local/esysbio/src/GenomeAnalysisTK/dist/GenomeAnalysisTK.jar'
#reference files
reference = '/export/astrakanfs/mpesj/reference/human_g1k_v37.clean.fasta'
dbsnp = '/export/astrakanfs/mpesj/reference/dbsnp_132_GRCh37.vcf'
#exome = '/export/astrakanfs/mpesj/reference/Agilent_SureSelect_v1_GRCh37.bed'
exome = '/export/astrakanfs/mpesj/reference/Nimblegen_SeqCap_EZ_Exome_v2_37_targetRegOnly_wingspan.bed'
#exome = '/export/astrakanfs/fro061/reference/Nimblegen3.0_new_samples_intersect.bed'

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
                        help="Path of alignment file in BAM format. ")
    



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

        raise Exception("Missing mandatory parameter%s: %s.\n\n%s\n\n" % ("s" if len(missing_options) > 1 else "",", ".join(missing_options),helpstr))
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
             -o %s" 
             % (gatk, bam, reference, intervals_file))
             
def run_realigner(bam, intervals_file, realigned_bam):
    """Performs local realignment of reads based on misalignments due to the presence of indels"""
    run_cmd("java -Xmx4g -jar %s \
             -T IndelRealigner \
             -I %s \
             -R %s \
             -targetIntervals %s \
             -o %s" 
             % (gatk, bam, reference, intervals_file, realigned_bam))

def dup_removal_picard(bam,output):
    """Use Picard to remove duplicates"""
    run_cmd('java -Xmx4096m -jar %s/CleanSam.jar INPUT=%s OUTPUT=%s VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR' 
                % (picard,bam,output))

def dup_removal_samtools(bam,output):
    """Use samtools for dupremoval"""
    run_cmd('samtools rmdup %s %s' % (bam,output))
                
def count_covariates(bam, recal_data):
    """Uses GATK to count covariates"""
    run_cmd("java -Xmx4g -jar %s \
            -T CountCovariates \
            -nt %s \
            -l INFO \
            -R %s \
            -knownSites %s \
            --default_platform illumina \
            -cov ReadGroupCovariate \
            -cov QualityScoreCovariate \
            -cov CycleCovariate \
            -cov DinucCovariate \
            -I %s \
            -recalFile %s"
            % (gatk, cpus, reference, dbsnp, bam, recal_data))
            
def table_recalibration(bam, recal_data, output):
    """uses GATK to rewrite quality scores using the recal_data"""
    run_cmd("java -Xmx4g -jar %s \
            -T TableRecalibration \
            --default_platform illumina \
            -R %s \
            --preserve_qscores_less_than 5 \
            -l INFO \
            -I %s \
            --out %s \
            -recalFile %s" 
            % (gatk, reference, bam, output, recal_data) )

def snp_call_samtools(bam, pileup):
    """Use samtools for snp calling (deprecated)"""
    prior="0.001"
    run_cmd("samtools pileup -vcs -r %s -l %s -f %s %s > %s" 
            % (prior, exome, reference, bam, pileup) )
            
def snp_call_mpileup(bam, bcf):
    """Uses samtools mpileup command for snp calling
        -C 50 reduces the effect of reads with excessive mismatches. 
        This aims to fix overestimated mapping quality and appears to be preferred for BWA-short
        Removed the -l option"""
    # get the calls with everything
    run_cmd("samtools mpileup -ugf %s -C 50 %s | bcftools view -bvcg - > %s"
            % (reference, bam, bcf))


def bcf_to_snp(bcf, snps):
    """Converts bcf to vcf file, only snp calls"""
    run_cmd("bcftools view -I %s | vcfutils.pl varFilter -D %d > %s" % (bcf, varfilter, snps))
    
def bcf_to_indels(bcf, indels):
    """Converts bcf to vcf, only indels"""
    run_cmd("bcftools view %s | vcfutils.pl varFilter -D %d | egrep '#|INDEL' > %s" % (bcf, varfilter, indels))

def snp_filter(vcf, output):
    """Apply filters to the vcf file to generate the final call of snps"""
    run_cmd("vcftools --vcf %s \
             --out %s \
             --recode \
             --bed %s \
             --keep-INFO-all "
            % (vcf, output, exome))

def indel_filter(vcf, temp, final_calls):
    """Apply filters to the vcf file to generate the final call of indels"""
    run_cmd("vcftools --vcf %s \
             --out %s \
             --recode \
             --bed %s \
	     --keep-INFO-all "
            % (vcf, temp, exome))
    run_cmd("mv %s.recode.vcf %s" %(temp, final_calls))

def pileup_to_vcf(pileup, filtered, prefix):
    """Convert the pileup file to vcf"""
    run_cmd("sam2vcf.pl --refseq %s < %s > %s.vcf.tmp"
             % (reference, pileup, pileup))
    # filter lines that contain D or N in the 5th column, substitute data for the sample name on the header
    run_cmd("awk {'if ($5 !~ /[(,D)N]/) print'} %s.vcf.tmp | awk '{gsub(\"data\",\"%s\",$0);print}' > %s.vcf.tmp2" % (pileup,prefix,pileup))
    #convert to vcd v4 (gatk does not accept v3 any more)
    run_cmd("cat %s.vcf.tmp2 | vcf-convert -r %s > %s" % (pileup,reference, filtered)) 

def annotate_vcf_gatk(vcf, bam, annotated_vcf):
    """Annotate the vcf file using gatk"""
    run_cmd("java -Xmx4g -jar %s \
                -l INFO \
                -T VariantAnnotator \
                -R %s \
                --dbsnp %s \
                -I %s \
                --variant %s \
                -o %s \
                -A AlleleBalance \
                -A MappingQualityZero \
                -A LowMQ \
                -A RMSMappingQuality \
                -A HaplotypeScore \
                -A QualByDepth \
                -A DepthOfCoverage \
                -A HomopolymerRun"
            % (gatk, reference, dbsnp, bam, vcf, annotated_vcf) )
            
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
@transform(link, suffix(".bam"), '.realign.intervals')
def find_realignment_intervals(bam,output):
   """Find regions to be re-aligned due to indels"""
   find_realigns(bam, output)

@follows(find_realignment_intervals)
@transform(find_realignment_intervals, suffix(".realign.intervals"), '.realigned.bam', r'\1.bam')
def indel_realigner(input, output, bam):
   """Re-aligns regions around indels"""
   run_realigner(bam,input,output)

@follows(indel_realigner)
@transform(indel_realigner, suffix(".realigned.bam"), '.rmdup.bam')
def remove_dups(input, output):
    """Remove dups"""
    dup_removal_samtools(input, output)

@follows('remove_dups')
@transform(remove_dups, suffix('.rmdup.bam'), '.gatk.bam.recal_data.csv')
def recalibrate_baseq1(input, output):
    """Base quality score recalibration in bam file
        Part 1: count covariates"""
    index_bam(input)
    count_covariates(input, output)

@follows(recalibrate_baseq1)
@transform(remove_dups, suffix('.rmdup.bam'), add_inputs(r'\1.gatk.bam.recal_data.csv'),'.gatk.bam')
def recalibrate_baseq2(inputs, output):
    """Base quality score recalibration in bam file
        Part 2: rewrite quality scores into a new bam file"""   
    table_recalibration(inputs[0], inputs[1], output)

@follows('recalibrate_baseq2')
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.bcf')
def call_snps(bam, bcf):
    """Call snps on the data, using samtools"""
    snp_call_mpileup(bam, bcf)

@follows('call_snps')
@transform(call_snps, suffix('.bcf'), '.snp.vcf')
def snps_to_vcf(bcf, vcf):
    """Converts bcf files to snp vcf"""
    bcf_to_snp(bcf,vcf)

@follows('call_snps')
@transform(call_snps, suffix('.bcf'), '.indel.vcf')
def indels_to_vcf(bcf, indel):
    """Converts bcf files to vcf, only indels"""
    bcf_to_indels(bcf, indel)

@follows('snps_to_vcf')
@transform(snps_to_vcf, suffix('.snp.vcf'), '.snp.recode.vcf',r'\1.snp')
def filter_snps(vcf, snps,extra):
    """Use vcftools for filtering the vcf output"""
    snp_filter(vcf, extra)

@follows('indels_to_vcf')
@transform(indels_to_vcf, suffix('.indel.vcf'), '.snp.recode.vcf',r'\1.indel')
def filter_indels(vcf, indel_calls,extra):
    """Use vcftools for filtering the vcf output, indels"""
    indel_filter(vcf, extra, indel_calls)

@follows('recalibrate_baseq2')
@transform(recalibrate_baseq2, suffix('.gatk.bam'), '.gatk.bam.bai')
def index2(bam, index):
    index_bam(bam)

@follows('filter_snps')
@transform(filter_snps, suffix('.snp.recode.vcf'), add_inputs(r'\1.gatk.bam'), '.snps')
def annotate_vcf(inputs, annotated_vcf):
    """Annotates snp vcf file"""
    annotate_vcf_gatk(inputs[0], inputs[1], annotated_vcf)

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

