# RNA-Seq-Pipeline

rnaSeqPipelineGLBRC.py

Purpose: Implementation of the Gasch lab RNA-Seq pipeline.

Input : A text file with RNA-Seq fastq files to be processed.

Please use a dedicated directory for running pipeline.
Create your directory and copy your fastq files into that directory.
Generate input file by moving into working directory and running /bin/ls *.fastq > input.txt

Required Parameters: -f input.txt
To run default enter:  /home/GLBRCORG/mplace/scripts/rnaSeqPipelineGLBRC.py -f input.txt

Optional Parameters:
	-r  this will use "-s reverse" parameter for HTSeq.

	-a  bwamem changes aligner from default bowtie2 to bwa mem

	-ref  change default reference, usage:  -ref y22-3
	    Current Reference List:
		R64-1-1 -- default equivalent to UCSC sacCer3
		Y22-3   -- GLBRC assembly 

Outline of steps & commands used in pipeline:

  1) Trimmomatic 
	/opt/bifxapps/Trimmomatic-0.30/trimmomatic-0.30.jar SE -phred33 input_fastq 
	outFile LEADING:3 TRAILING:3 SLIDINGWINDOW:3:30 MINLEN:36

  2) Fastqc 
	/opt/bifxapps/FastQC/fastqc trimmed_fastq

  3) Alignment 
	Bowtie2 -p 8 --phred33 -N 1 -x  referenceFile -U  trimmed_fastq -S outFile
	or
	bwa mem -t 8 -M referenceFile trimmed_fastq

  4) Picard-tools 
	/opt/bifxapps/picard-tools/CleanSam.jar I= samFile O=outFile
	/opt/bifxapps/picard-tools/AddOrReplaceReadGroups.jar 'I=', inFile 
	O=outFile SO=coordinate LB=S288C_reference_sequence_R64-1-1_20110203.fasta
	PL=ILLUMINA PU=unk RGSM=SampleName

  5) samtools 
	sam to bam: samtools view -bS -t reference.fsa.fai -o bam samFile
	  sort bam: samtools sort input_bam sorted_Bam 
	 index bam: samtools index sorted_Bam

  6) HTSeq 
	/opt/bifxapps/python/bin/htseq-count -t CDS -i Parent inputFile Ref_gff
	 -- if '-r' specfied -s reverse will be used as well

  7) RPKM 
	/home/GLBRCORG/mplace/bin/Normalize.jar dir Ref_gff RPKM.results --gene

  8) bam2wig.pl 
	/opt/bifxapps/biotoolbox/scripts/bam2wig.pl --in bamFile --pos mid --strand --rpm --out outFile

  9) findreplace_WIG.pl
	/home/GLBRCORG/mplace/scripts/findreplace_WIG.pl

The results are organized into the following subdirectories:
	 alignments/ fastq/ fastqc/ htseq/ log/ wig/ 
	 RPKM results are written to a file called: RPKM.results



This script is designed to run on the GLBRC scarcity servers 6-10.
See Mike Place for problems with this script.
