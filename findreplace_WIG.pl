#!/usr/bin/perl
use strict;
use warnings;

#This script will collect all .WIG files and convert the chromosome names for use in MochiView
#This can also be modified to find and replace any characters in all files of a file type in a directory
#Written by Kevin Myers (kmyers2@wisc.edu) 12-18-2014

#Modified 3/17/2015 by Mike Place for use with the RNA-Seq pipeline
#The wig files output by the rnaSeqPipelineGLBRC.py are gzipped, so we 
#need to decompress them prior to replacing the chromosome names.

#Input arguments and variables and collect WIG files in directory.

my $INPUT_DIR = shift;
die "Undefied input directory!\n" unless defined $INPUT_DIR;

opendir my ($IN), $INPUT_DIR or die "Cannot open input directory '$INPUT_DIR'\n";

my $cmd = 'gzip -d *.wig.gz';
system($cmd);

# Find all WIG files from the current directory and store in @INfiles

my @INfiles = grep { /\.wig$/} readdir $IN;
die "No .wig files found in this directory!\n" unless @INfiles > 0;

#Open and define input and output files from the array containing WIG files in the specified directory
foreach my $file (@INfiles) {
	if ($file =~/(.*)\.wig$/) {
	open(IN,"<$file") or die "$file not found!\n";
	my $file_out = "$file"."_mochiview.wig";
	open(OUT,">$file_out") or die " $file_out not found!\n";

#Search for the GenBank chromosome names and replace them with MochiView chromosome names
	while (<IN>) {
   		s/chrom\=ref\|NC_001133\|/chrom\=chr1/g; # do the first replacement
   		s/chrom\=ref\|NC_001134\|/chrom\=chr2/g; # do the second replacement and so on
   		s/chrom\=ref\|NC_001135\|/chrom\=chr3/g;
   		s/chrom\=ref\|NC_001136\|/chrom\=chr4/g;
   		s/chrom\=ref\|NC_001137\|/chrom\=chr5/g;
   		s/chrom\=ref\|NC_001138\|/chrom\=chr6/g;
   		s/chrom\=ref\|NC_001139\|/chrom\=chr7/g;
   		s/chrom\=ref\|NC_001140\|/chrom\=chr8/g;
   		s/chrom\=ref\|NC_001141\|/chrom\=chr9/g;
   		s/chrom\=ref\|NC_001142\|/chrom\=chr10/g;
   		s/chrom\=ref\|NC_001143\|/chrom\=chr11/g;
   		s/chrom\=ref\|NC_001144\|/chrom\=chr12/g;
   		s/chrom\=ref\|NC_001145\|/chrom\=chr13/g;
   		s/chrom\=ref\|NC_001146\|/chrom\=chr14/g;
   		s/chrom\=ref\|NC_001147\|/chrom\=chr15/g;
   		s/chrom\=ref\|NC_001148\|/chrom\=chr16/g;
   		s/chrom\=ref\|NC_001224\|/chrom\=chrM/g;

   print OUT; # print to the modified output file
		}
	# overwrite the original file
	rename $file_out, $file;
	}
}
close(IN);
close(OUT);
