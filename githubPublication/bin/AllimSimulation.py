#!/usr/bin/python
import sys
import os
import re
#import classes.utility
import classes.simulateRNAseq
#import classes.filter_reads
#if len(sys.argv)<3:
#	print "python", sys.argv[0], "outDir genome1Name genome2Name genome1Fasta genome2Fasta readlength readqual error readtype insert "
#	os.exit(0);
simulated_reads_dir="/home/uqhshao/30days/FruitFly/pipeline/github/bin/"
genome1="P1"
transcript_fasta_Parent1="/home/uqhshao/30days/FruitFly/pipeline/github/TestData/line190_GATK_noINDEL.fasta"
readlength=151
readqual="e"
error=0.0
readtype="pe"
insert=350
#reuse the simulation fuction in Allim
classes.simulateRNAseq.simulateRNAseq.simulate_reads(simulated_reads_dir,genome1,transcript_fasta_Parent1,genome1,readlength,readqual,error,readtype,insert)
