#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
Program: rnaSeqPipelineGLBRC.py

Purpose: Implement the currently used Gasch lab RNA-Seq pipeline. 

Input : text file with RNA-Seq fastq files to be processed
        to generate:  /bin/ls *.fastq > input.txt
        optional parameters: -htseq reverse   ( for HTSeq )
   
Output: Each step has its own output see below.

Steps:
    Trimmomatic         -- http://www.usadellab.org/cms/?page=trimmomatic
    Fastqc              -- http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    bowtie2             -- http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    bwa mem             -- http://bio-bwa.sourceforge.net/
    bam2wig             -- http://search.cpan.org/~tjparnell/Bio-ToolBox-1.24001/lib/Bio/ToolBox.pm
    picard/CleanSam.jar -- http://broadinstitute.github.io/picard/
    samtools            -- http://www.htslib.org/
    HTSeq               -- http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
    Normalization(RPKM) 

*******************************************************************************
Trimmomatic:
        HEADCROP       = 5  
        LEADING        = 3 
        TRAILING       = 3
        SLIDINGWINDOW  = 3:30   (window = 3, min avg quality for window = 30)
        MINLEN         = 36        
        -phred33 
        -threads 8
        SE is single-end
        PE is paired-end
        -trimlog <logfile>

    Single-end
    java -jar <path to trimmomatic jar> SE [-threads <threads>] -phred33 
    [-trimlog <logFile>] <input> <output> 

    on GLBRC scarcity:
    java -jar  /opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar

    OUTPUT: trimmed fastq files
    
*******************************************************************************

FastQC:
    To run non-interactively simply a list of files to process on the commandline

    fastqc somefile.txt someotherfile.txt

    You can specify as many files to process in a single run as you like

    If you want to save your reports in a folder other than the folder 
    
    --outdir=/some/other/dir/
    -quiet   = only report errors
    
    on GLBRC scarcity:
    /opt/bifxapps/FastQC/fastqc --help
    
    OUTPUT: QC report after trimming with Trimmomatic

*******************************************************************************
MAPPING:  default is Bowtie2, but user may choose bwa mem
*******************************************************************************
Bowtie2:
    -p 8           # specified number of parallel search threads
    --phred33 
    -N 1           # Sets the number of mismatches to allowed in a seed alignment
                     during multiseed alignment. Can be set to 0 or 1.
    -x $REF        # index reference

    -U $READS      # read
    -S $OUT.sam    # output sam file

    on GLBRC scarcity:
    /opt/bifxapps/bin/bowtie2
    
    OUTPUT: sam file

Bwa mem:
    -t 8          # number of threads
    -M $REFERENCE # reference genome file
  example:
    bwa mem -t 8 -M $REFERENCE $file 1> $out.sam 2>>$OUTPUT_DIR/Bwa_run.log
    
    OUTPUT: sam file

*******************************************************************************
picard

 Clean the SAM file
    This soft-clips an alignment that hangs off the end of its reference sequence.
    This will print out all the errors that it ignores (MAPQ errors)

java -Xmx15g -jar /opt/bifxapps/picard/CleanSam.jar I=$OUT.sam O=$OUT.cleaned.sam
rm $OUT.sam

 Add the RG header and sort the SAM file
    This will print out all the errors that it ignores (MAPQ errors)

java -Xmx15g -jar /opt/bifxapps/picard/AddOrReplaceReadGroups.jar I=$OUT.cleaned.sam 
O=$OUT.final.sam SO=coordinate LB=$REF.fasta PL=ILLUMINA PU=unknown SM=$OUT 
VALIDATION_STRINGENCY=LENIENT

rm $OUT.cleaned.sam
    on GLBRC scarcity:
    java -Xmx15g -jar /opt/bifxapps/picard/AddOrReplaceReadGroups.jar
    
    OUTPUT: sam file
*******************************************************************************
samtools

Make the BAM file, sort and index it
    samtools view -uS -t $REF.fasta.fai $OUT.final.sam | samtools sort - $OUT.sorted

samtools index $OUT.sorted.bam 

    on GLBRC scarcity:
    /opt/bifxapps/bin/samtools
    
    OUTPUT: sorted bam file and index for bam file
*******************************************************************************
bam2wig.pl

This script will convert alignments from a Bam file into enumerated point data
in a wig format.

bam2wig.pl --in bamFile --pos mid --strand --rpm --out

This works on  scarcity-1, scarcity-5, scarcity-6 but may fail elsewhere

    OUTPUT:  gzipped wig file

*******************************************************************************
HTSeq

Given a file with aligned sequencing reads and a list of genomic features, 
count how many reads map to each feature.

htseq-count -t CDS -i Parent samFile  gff 

    OUTPUT: htseq text file

*******************************************************************************
RPKM

RPKM normalization - Normalize.jar from Nikolay

java -Xmx8g -jar Normalize.jar cwd  gff RPKM.results --gene

    OUTPUT: RPKM text file , one column per sample
    
*******************************************************************************

REQUIRED:
    Picard tools
    bowtie2
    fastqc
    trimmomatic
    python (snakemake)
    bowtie2-build <reference.in> <basename for index files>
    samtools faidx <referenc.in>

@author: mplace
"""
import os
import re
import sys
import subprocess      
import argparse            

# reference dictionary to stores locations of reference fasta, gff, dict files
# key = reference name, value is a tuple where the order is defined as:
#   [0] = bowtie2 reference
#   [1] = bwa mem reference
#   [2] = gff file
#   [3] = samtools index
#   [4] = picard used for nameing
# default reference = R64 (SGD R64-1-1 = UCSC sacCer3)
# Y22 = reference of S. cerevisiae Y22-3 GLBRC sequenced strain

ref = { 'R64' : ( "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/s.cerevisiae-R64-1-1",
                  "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta",
                  "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff",
                  "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa.fai",
                  "/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta"),
        
        'Y22' : ( "/home/GLBRCORG/mplace/data/reference/Y22-3/Y22-3-bowtie",
                  "/home/GLBRCORG/mplace/data/reference/Y22-3/Y22-3.fasta",
                  "/home/GLBRCORG/mplace/data/reference/Y22-3/Y22-3_Final_GFF.gff",
                  "/home/GLBRCORG/mplace/data/reference/Y22-3/Y22-3.fasta.fai",
                  "/home/GLBRCORG/mplace/data/reference/Y22-3/Y22-3.fasta" )
        }

def runTrimmomatic( fastq ):
    """
    Run trimmomatic on fastq file
    java -Xmx6g -jar ~/bin/trimmomatic SE -phred33 
    -trimlog trimlog.out run333.YPS1009.10kreads.fastq.gz
    trimmed.out.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:3:30 MINLEN:36
    
    """
    program = '/opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar'
    outFile = re.sub(r"fastq","trim.fastq", fastq)
    logfile = fastq.rstrip('fastq')
    logfile += "trim.log"    
    
    cmd = ['java', '-Xmx6g', '-jar', program , 'SE', '-phred33', 
            fastq, outFile, 'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:3:30', 'MINLEN:36' ] 
                      
    output = subprocess.Popen( cmd, stderr=subprocess.PIPE).communicate()
    result = output[1].decode('utf-8')              # must explicitly convert bytes to unicode
    with open('Trimmomatic.log', 'a') as log:
        log.write(result)
        log.write("\n")
    return outFile                                  # name of trimmed fastq

def runFastqc( fastq ):
    """
    run Fastqc on the trimmomatic results
    /home/mplace/bin/fastqc
    """   
    program = '/opt/bifxapps/FastQC/fastqc'
    cmd =  [ program , fastq ]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    result1 = output[0].decode('utf-8')
    result2 = output[1].decode('utf-8')
    with open('Fastqc.log', 'a') as log:
        log.write(result1)
        log.write(result2)
        log.write("")
        log.write("")
                           
def runBowtie2( fastq, refer ):
    """
    Run bowtie2 for alignment
    """    
    program = 'bowtie2'
    outFile = re.sub(r"fastq", "sam", fastq )
    cmd = [ program , '-p',  '8', '--phred33',  '-N',  '1', '-x', ref[refer][0],
           '-U', fastq, '-S', outFile ]
    output = subprocess.Popen( cmd, stderr=subprocess.PIPE ).communicate()
    result = output[1].decode( 'utf-8' )
    with open( 'Bowtie2.log', 'a' ) as log:
        log.write("%s\n" %(outFile) )
        log.write(result)
        log.write("\n\n")
                       
    return outFile             

def runBwaMem( fastq, refer ):
    """
    Run default bwa mem for alignment
    
    bwa mem -t 8 -M $REFERENCE $file 1> $out.sam
    reference: /home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta
    """
    program = 'bwa'
    outFile = re.sub(r"fastq", "sam", fastq )    
    cmd     = [ program, 'mem', '-t', '8', '-M', ref[refer][1], fastq ]
#    with open(outFile, "w") as out:
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    result = output[0].decode( 'utf-8' )
    log    = output[1].decode( 'utf-8' )    
    # write sam file
    with open( outFile, 'w' ) as out:
        out.write(result)
    out.close()
    # write bwa mem log results
    with open( 'bwamem.log', 'a' ) as logout:
        logout.write( "%s\n" %(outFile) )
        logout.write( log )
        logout.write( "\n\n" )
    logout.close()
    
    return outFile        

def bam2Wig( bamFile ):
    """
    Run bam2wig.pl, convert alignments from a Bam file into enumerated
    point data in a wig format. 
    bam2wig.pl --in run333.YPS163.10kreads.sort.bam --pos mid --strand --rpm --out YPS163.wig
    """
    program = '/opt/bifxapps/biotoolbox/scripts/bam2wig.pl'
    outFile = re.sub(r"bam", "wig", bamFile)
    cmd = [ program , '--in', bamFile, '--pos', 'mid', '--strand', '--rpm', '--out',  outFile]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
    result1 = output[0].decode( 'utf-8')
    result2 = output[1].decode( 'utf-8' )
    
    with open( 'bam2wig.log', 'a' ) as log:
        log.write("%s\n" %(outFile) )
        log.write(result1)
        log.write(result2)
        log.write("\n\n")    
    

def runPicard( samFile, refer ):
    """
    Run picard tools on a sam file created by Bowtie2.
    java -Xmx8g -jar /home/mplace/bin/picard-tools-1.119/CleanSam.jar
    #'/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta',
    """
    program = '/opt/bifxapps/picard-tools/CleanSam.jar'
    sample  = samFile.split('.')
    outFile = re.sub( r"sam", "clean.sam", samFile )
    cmd = [ 'java', '-Xmx8g', '-jar', program ,'I=', samFile, 'O=', outFile]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
    result = output[1].decode( 'utf-8' )
    with open( 'Picard.log', 'a' ) as log:
        log.write("%s\n" %(outFile) )
        log.write(result)
        log.write("\n\n")

    program = '/opt/bifxapps/picard-tools/AddOrReplaceReadGroups.jar'
    inFile  = outFile
    outFile = re.sub(r"clean.sam", "final.sam", inFile )
    cmd = [  'java', '-Xmx8g', '-jar', program, 'I=', inFile, 'O=', outFile, 'SO=', 'coordinate', 'LB=',
             ref[refer][4], 'PL=', 'ILLUMINA', 'PU=unk', 'RGSM=', sample[1] ]
    
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
    result = output[1].decode( 'utf-8' ) 
    with open( 'Picard.log', 'a' ) as log:
        log.write("%s\n" %(outFile) )
        log.write(result)
        log.write("\n\n")
        
     #remove initial samFile here      
    if os.path.exists( inFile ):
        os.remove( inFile )
 
    return outFile

def runSamtools( samFile, refer ):
    """
    Run samtools to sort, index and produce a bam file
    samtools view -uS -t $REF.fasta.fai $OUT.final.sam | samtools sort - $OUT.sorted
    samtools index $OUT.sorted.bam
    #'/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa.fai',
    """
    sortSam   = re.sub(r"trim.final", "sort",     samFile)
    bam       = re.sub(r"trim.final.sam", "bam",  samFile)
    sortBam   = re.sub(r"trim.final.sam", "sort", samFile)
    deleteSam = re.sub(r"final.", "", samFile)

    # convert sam to bam
    cmd = [ 'samtools', 'view', '-bS', '-t', ref[refer][3], '-o', bam,  samFile ]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
    result = output[1].decode( 'utf-8' ) 
    with open( 'samtools.log', 'a' ) as log:
        log.write("convert to bam %s\n" %(samFile) )
        log.write(result)
        log.write("\n\n")
        
    os.unlink( samFile )
    # sort bam file
    cmd = [ 'samtools', 'sort', bam, sortBam ]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
    result = output[1].decode( 'utf-8' ) 
    with open( 'samtools.log', 'a' ) as log:
        log.write("sort bam %s\n" %(bam) )
        log.write(result)
        log.write("\n\n")
    
    # create wigfile from sorted bam file
    wigBam = sortBam + ".bam"
    bam2Wig( wigBam )
    
    # index bam file
    cmd = ['samtools', 'index', wigBam]
    output = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
    result = output[1].decode( 'utf-8' )
    with open( 'samtools.log', 'a' ) as log:
        log.write("samtools index %s\n" %(wigBam) )
        log.write(result)
        log.write("\n\n")

    # delete first bam file
    os.unlink(bam)
    os.unlink(deleteSam)    

def runHTSeq( cwd, refer, rvse ):
    """
    Call HTSeq-count
    /home/GLBRCORG/mplace/test/forkevin 
    add -s reverse as parameter
    #"/home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208_noFasta.gff"
    """
    gff    = ref[refer][2] 
    program  = "/opt/bifxapps/python/bin/htseq-count"
    
    bamFiles = [ fn for fn in os.listdir(cwd) if fn.endswith(".bam") ]         # Get a list of all bam files in current directory
    
    for i in bamFiles:
        samName  = re.sub(r"bam", "sam", i ) 
        htseqOut    = i + "_HTseqOutput.txt"         
        samcmd  = [ 'samtools', 'view', '-o', samName, i ]
        output  = subprocess.Popen( samcmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
        
        if rvse == 1:
            cmd     = [ program, '-t', 'CDS', '-i', 'Parent', '-s', 'reverse' , samName, gff ]
        else:
            cmd     = [ program, '-t', 'CDS', '-i', 'Parent' , samName, gff ]
        output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
        result1 = output[0].decode( 'utf-8' )
        result2 = output[1].decode( 'utf-8' )
        
        with open( htseqOut, 'w' ) as out:
            out.write( result1 )
        
        with open( 'HTSeq.log', 'a' ) as log:
            log.write( "HTSeq for %s\n" %( samName ) )
            log.write(" ".join(cmd) )
            log.write( result2 )
            log.write("\n")
        
        os.unlink(samName)
        
def runRPKM( cwd, refer):
    """
    Run RPKM normalization on all HTSeq files in current working directory.
    java -Xmx6g -jar ~/bin/Normalize.jar . /home/GLBRCORG/mplace/data/reference/S288C_reference_genome_R64-1-1_20110203/saccharomyces_cerevisiae_R64-1-1_20110208.gff test.out --gene
    """
    program = '/home/GLBRCORG/mplace/bin/Normalize.jar'
    gff     =  ref[refer][2] 
    cmd     = [ 'java', '-Xmx8g', '-jar', program, cwd, gff, 'RPKM.results', '--gene' ]
    output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE ).communicate()
    result  = output[0].decode( 'utf-8' )    
    with open( 'RPKM.log', 'w' ) as log:
        log.write(result)
     
def cleanUp( cwd ):
    """
    Clean up and move output files.
    os.mkdir()
    os.rename( currentPath/, newPath )
    """
    cwd = cwd + "/"
    # move bam files to alignment directories
    os.mkdir( "alignment" )
    bamDir = cwd + "/alignment/"
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".bam") ]
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".bai") ] 
    # make wig directory
    os.mkdir( "wig" )
    widDir = cwd + "/wig/"
    [ os.rename( (cwd + fn), (widDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".wig.gz") ]
    # make HTseq directory
    os.mkdir( "htseq" )
    htsDir = cwd + "/htseq/"
    [ os.rename( (cwd + fn), (htsDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("_HTseqOutput.txt") ]
    # make log file directory
    os.mkdir( "log" )
    logDir = cwd + "/log/"
    [ os.rename( (cwd + fn), (logDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".log") ]
    # make fastqc dir
    os.mkdir( "fastqc" )
    qcDir = cwd + "/fastqc/" 
    [ os.rename( (cwd + fn), (qcDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("_fastqc.zip") ]
    [ os.rename( (cwd + fn), (qcDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("_fastqc") ] 
    # make a directory for sequence reads
    os.mkdir( "fastq" )
    seqDir = cwd + "/fastq/"
    [ os.rename( (cwd + fn), (seqDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".fastq") ]

def replaceWig( cwd ):
    """
    Run Kevin's script to replace chromosome names in wig file, for compatiblity with Mochi view.
    """
    wigDir  = cwd + "/wig/"
    os.chdir(wigDir)
    program = '/home/GLBRCORG/mplace/scripts/findreplace_WIG.pl'
    cmd     = [ program, wigDir ]
    output  = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    result  = output[0].decode( 'utf-8' )
    with open( cwd + "/log/" + 'findreplace_WIG.log', 'w' ) as log:
        log.write(wigDir)
        log.write("\n")
        log.write("\n".join(cmd) )
        log.write("\n")
        log.write(result)

def main():
    
    cmdparser = argparse.ArgumentParser(description="RNA-Seq alignment, HTSeq & RPKM pipeline.",
                                        usage='%(prog)s -f <fastq file list.txt> [optional args: -r -d ]' ,prog='rnaSeqPipelineGLBRC.py'  )
    cmdparser.add_argument('-f', '--file',    action='store', dest='FILE',    help='Text file, one fastq file name per line.')
    cmdparser.add_argument('-r', '--reverse', action='store_true', dest='REVERSE', help='HTSeq -s reverse, for Biotech GEC data, optional.')
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest='DETAIL',  help='Print a more detailed description of program.')
    cmdparser.add_argument('-a', '--aligner',  action='store', dest='ALIGNER', help='Default aligner is Bowtie2, to use Bwa mem: -a bwamem')
    cmdparser.add_argument('-ref', '--reference', action='store', dest='REFERENCE', help='To Change default (SGD R64-1-1) reference: -ref Y22-3 ')
    cmdparser.add_argument('-trim', '--trimAs', action='store_true', dest='TRIM', help='Trim consecutive 3 prime A\'s from reads ')
    cmdResults = vars(cmdparser.parse_args())
        
    fastq = []                      # list of fastq files to process
    cwd   = os.getcwd() 
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
        
    if cmdResults['DETAIL']:
        print("\nrnaSeqPipelineGLBRC.py")
        print("\nPurpose: Implementation of the Gasch lab RNA-Seq pipeline.")
        print("\nInput : A text file with RNA-Seq fastq files to be processed.")
        print("\nPlease use a dedicated directory for running pipeline.")
        print("Create your directory and copy your fastq files into that directory.")
        print("Generate input file by moving into working directory and running /bin/ls *.fastq > input.txt\n")
        print("Required Parameters: -f input.txt")
        print("To run default enter:  /home/GLBRCORG/mplace/scripts/rnaSeqPipelineGLBRC.py -f input.txt\n")
        print("Optional Parameters:") 
        print("\t-r  this will use \"-s reverse\" parameter for HTSeq.\n")
        print("\t-a  bwamem changes aligner from default bowtie2 to bwa mem\n")
        print("\t-ref  change default reference, usage:  -ref y22-3")
        print("\t    Current Reference List:")
        print("\t\tR64-1-1 -- default equivalent to UCSC sacCer3" )
        print("\t\tY22-3   -- GLBRC assembly \n")       
        print("Outline of steps & commands used in pipeline:\n")
        print("  1) Trimmomatic ")
        print("\t/opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar SE -phred33 input_fastq ")
        print("\toutFile LEADING:3 TRAILING:3 SLIDINGWINDOW:3:30 MINLEN:36\n")
        print("  2) Fastqc ")
        print("\t/opt/bifxapps/FastQC/fastqc trimmed_fastq\n")
        print("  3) Alignment ")
        print("\tBowtie2 -p 8 --phred33 -N 1 -x  referenceFile -U  trimmed_fastq -S outFile" )
        print("\tor")
        print("\tbwa mem -t 8 -M referenceFile trimmed_fastq\n")
        print("  4) Picard-tools ")
        print("\t/opt/bifxapps/picard-tools/CleanSam.jar I= samFile O=outFile")
        print("\t/opt/bifxapps/picard-tools/AddOrReplaceReadGroups.jar 'I=', inFile ")
        print("\tO=outFile SO=coordinate LB=S288C_reference_sequence_R64-1-1_20110203.fasta\n\tPL=ILLUMINA PU=unk RGSM=SampleName\n")
        print("  5) samtools ")
        print("\tsam to bam: samtools view -bS -t reference.fsa.fai -o bam samFile")
        print("\t  sort bam: samtools sort input_bam sorted_Bam ")
        print("\t index bam: samtools index sorted_Bam\n")
        print("  6) HTSeq ")
        print("\t/opt/bifxapps/python/bin/htseq-count -t CDS -i Parent inputFile Ref_gff")
        print("\t -- if '-r' specfied -s reverse will be used as well\n")
        print("  7) RPKM ")
        print("\t/home/GLBRCORG/mplace/bin/Normalize.jar dir Ref_gff RPKM.results --gene\n")
        print("  8) bam2wig.pl ")
        print("\t/opt/bifxapps/biotoolbox/scripts/bam2wig.pl --in bamFile --pos mid --strand --rpm --out outFile\n")
        print("  9) findreplace_WIG.pl")
        print("\t/home/GLBRCORG/mplace/scripts/findreplace_WIG.pl\n")
        print("The results are organized into the following subdirectories:")
        print("\t alignments/ fastq/ fastqc/ htseq/ log/ wig/ ")
        print("\t RPKM results are written to a file called: RPKM.results\n")
        print("\n")
        print("This script is designed to run on the GLBRC scarcity servers 6-10.")
        print("See Mike Place for problems with this script.")
        sys.exit(1)

    if cmdResults['REFERENCE'] == 'Y22-3':
        reference = 'Y22'
    else:
        reference = 'R64'

    if cmdResults['ALIGNER'] is not None:
        aligner = cmdResults['ALIGNER']
    else:
        aligner='bowtie2'

    if cmdResults['FILE'] is not None:
        inFile = cmdResults['FILE']        
        with open( inFile, 'r') as f:
            for line in f:
                line = line.rstrip()
                fastq.append(line)
        # process each fastq file in list
        for item in fastq:
            trimFastq    = runTrimmomatic( item )
            runFastqc( trimFastq  )
            if re.match(aligner, 'bwamem'):
                samFile      = runBwaMem( trimFastq, reference)
            else:
                samFile      = runBowtie2( trimFastq, reference )  # add parameter
            finalSamFile = runPicard( samFile, reference )
            runSamtools( finalSamFile, reference )
        
        if cmdResults['REVERSE']:
            runHTSeq( cwd, reference, rvse=1 )
        else:
            runHTSeq  ( cwd, reference, rvse=0)
        runRPKM   ( cwd, reference )
        cleanUp   ( cwd )  
        replaceWig( cwd )

if __name__ == "__main__":
    main()

