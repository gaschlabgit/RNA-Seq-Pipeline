#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
@Program: remove3primeAFastq.py

@Purpose: Remove all consecutive 3 prime A's from fastq file, also trims quality scores.
          Remove all consecutive 3 prime A's from fastq file where pattern is something like:
          
          nnnnnnnnnnnnnnnAAAAAAAGATCGGAAGAGC
          
@Input:  Input file listing one fastq file per line, unzipped and the number of A's to search for

@Dependencies: Python 3

@Output:  Fastq file with 3 prime A's removed
@author: Mike Place
@Date:   2/5/2015

Example:

~/scripts/remove3primeAFastq.py <listfile> <num A> 

ORIGINAL:

@HWI-D00256:264:HBC0UADXX:2:1101:11620:2235 1:N:0:CGTGAT
GGGGCAGCTTTCAGATTGTACTAAGCATAATTACGTGTTTTCATAGTTTAACGCTTTCAGAACTACTTATTTAATTTTGTAAGAAGTAATTTGAGTCACAT
+
?@@?D@@;CFBFDD<CBC<<<AA<GGHEGBHIIDG1???FGBBDGGIGGGED?@GDFFHCGDHGIGGCDEEEEC?CHHEHHB>?DE@@ADEE@C@CCCCAC

@HWI-D00256:264:HBC0UADXX:2:1101:14412:2175 1:N:0:CGTGAT
GGGATGTTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
+
CCCFFFFFHHHHHJJJJJJJJJJHFDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

TRAILING A's removed:

@HWI-D00256:264:HBC0UADXX:2:1101:11620:2235 1:N:0:CGTGAT
GGGGCAGCTTTCAGATTGTACTAAGCATAATTACGTGTTTTCATAGTTTAACGCTTTCAGAACTACTTATTTAATTTTGTAAGAAGTAATTTGAGTCACAT
+
?@@?D@@;CFBFDD<CBC<<<AA<GGHEGBHIIDG1???FGBBDGGIGGGED?@GDFFHCGDHGIGGCDEEEEC?CHHEHHB>?DE@@ADEE@C@CCCCAC


@HWI-D00256:264:HBC0UADXX:2:1101:14412:2175 1:N:0:CGTGAT
GGGATGTTT  -- -- SEE THE A's have been removed
+
CCCFFFFFH  -- Quality score also adjusted




"""
import itertools
import re
import sys

def stripA( dna, qual ):
    """
    strip A's from the end of string
    """
    m = re.search("A+$", dna)
    if m:
        start = m.start()
        trimA = dna[:start]
        trimQ = qual[:start]
        return trimA, trimQ
    else:
        return dna, qual
        
def stripAplusEnd( dna, qual, numA ):
    """
    strip   a pattern like:   NNNNNN.....AAAAAAAAGATCGGAAGA from 3' end of seq
    """
    #******  CHANGE THE DESIRED NUMBER OF A's HERE 
    exp = r"[A]{%d,}[GATC]{5,}$" %(numA)    
    m = re.search(exp, dna)
    if m:
        start = m.start()
        trimA = dna[:start]
        trimQ = qual[:start]
        return trimA, trimQ
    else:
        return dna, qual
    
def main():
    """
    main() 
    """
  
    # if no args print help
    if len(sys.argv) <= 2:
        print("Input file and Number of A's to remove are required.")
        sys.exit(1)
    else:
        inFile = sys.argv[1]
        numA   = int(sys.argv[2])
    
    fastqList = []
    
    with open( inFile, 'r') as f:
            for line in f:
                line = line.rstrip()
                fastqList.append(line)   
    
    for item in fastqList:
        removeAFastq = item.replace("fastq", "rmA.fastq")
    
        with open(removeAFastq, "w") as out:
            with open(item) as f:
                for line1, line2, line3, line4 in itertools.zip_longest(*[f]*4):
                    header = line1.rstrip()
                    dna    = line2.rstrip()
                    qual   = line4.rstrip()
            
                    dna_A, qual_A = stripA(dna, qual)                  # strip A's first
                    
                    dna_End, qual_End = stripAplusEnd(dna_A, qual_A, numA)   # strip AAAGATC... ends
            
                    out.write("%s\n" %(header) )
                    out.write("%s\n" %(dna_End) )
                    out.write("+\n")
                    out.write("%s\n" %(qual_End))  
        out.close()
    

if __name__ == "__main__":
    main()

