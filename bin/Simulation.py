#!/usr/bin/python
import sys
import os
import re
#import classes.utility
import classes.simulateRNAseq
#import classes.filter_reads
if len(sys.argv)<8:
	print "python", sys.argv[0], "outDir genome1Name genome2Name genome1Fasta genome2Fasta readlength insert "
	sys.exit(0);
#reuse the simulation fuction in Allim
simulated_reads_dir=sys.argv[1]
genome1=sys.argv[2]
genome2=sys.argv[3]
transcript_fasta_Parent1=sys.argv[4]
transcript_fasta_Parent2=sys.argv[5]
readlength=sys.argv[6]
readqual="e"
error=0.0
readtype="pe"
insert=sys.argv[7]
classes.simulateRNAseq.simulateRNAseq.simulate_reads(simulated_reads_dir,genome1,transcript_fasta_Parent1,genome1,readlength,readqual,error,readtype,insert)
classes.simulateRNAseq.simulateRNAseq.simulate_reads(simulated_reads_dir,genome2,transcript_fasta_Parent2,genome2,readlength,readqual,error,readtype,insert)
