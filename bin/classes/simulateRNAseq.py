#!/usr/bin/python
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import sys
import os
import re # to support perl like regular expression
import numpy
import collections

class simulateRNAseq:
    
    @classmethod
    def read_gtf(cls,gtffile):
        
        if not os.path.isfile(gtffile):
            sys.stderr.write("\nERROR --> "+str(gtffile)+" file does not exists\n\n")
            sys.exit(1)
                        
        fh = open(gtffile,"r")
        transcript = collections.defaultdict(lambda:"")
    
        for l in fh:
            l = l.strip()
            #print l
            a = l.split("\t")
            last_col = a[-1]
            start = a[3]
            end = a[4]
            strand = a[6]
            chromosome = a[0]
            
            if (a[2]=="exon"):
                if strand=="+" or strand=="-":
                    (geneid,transid) = simulateRNAseq.parse_anno(last_col)
                    key=geneid+"\t"+transid
                    transcript[key] += str(start)+".."+str(end)+".."+str(strand)+".."+str(chromosome)+"#"
        
        longest_isoform = simulateRNAseq.get_longest_isoform(transcript)
    
        #for k in longest_isoform.keys():
        #    print k,longest_isoform[k]
        return longest_isoform
    
    @classmethod
    def get_longest_isoform(cls,transcript):
        
        genes = collections.defaultdict(dict)
        genes_coord = collections.defaultdict(dict)
        
        for k in transcript.keys():
            
            geneid,tranid = k.split("\t")
            exons = transcript[k].split("#")
    
            start=[]
            end=[]
            strand=""
            chromosome=""
            for exon in exons:
                if exon:
                    st,nd,strand,chromosome = exon.split("..");
                    start.append(int(st))
                    end.append(int(nd))
            
            start =  sorted(start,reverse=False) 
            end = sorted(end,reverse=False)
    
            mrna_coord = []
            index=0
            for s in start:
                mrna_coord.append(str(s)+".."+str(end[index]))
                index +=1
            
            mrna_coord.append(str(strand))
            mrna_coord.append(str(chromosome))
            mrna_start = int(start[0])
            mrna_end = int(end[-1])
            mrna_length = (mrna_end-mrna_start)+1
            
            genes[geneid][tranid]=int(mrna_length)
            genes_coord[geneid][tranid]=mrna_coord
            
            mrna_coord=[]
            
        ### getting the longest isoform and strand
        longest_isoform = collections.defaultdict(dict)
        for k in genes.keys():
            longest_transcript =  max(genes[k], key=genes[k].get)
            longest_isoform[k][longest_transcript] = genes_coord[k][longest_transcript]
        
        ## clear the intermediate dicts
        transcript.clear()
        genes.clear()
        genes_coord.clear()
    
        return  longest_isoform

    @classmethod
    def parse_anno(cls,string):
        geneid = ""
        transid = ""
        a = string.split(";")
        for f in a :
            f = f.strip()
            genematch = re.match(r'gene_id',f)
            transmatch = re.match(r'transcript_id',f)
                
            if genematch:
                geneid = f
                geneid = geneid.replace("gene_id","")
                geneid = geneid.replace("\"","")
                geneid = geneid.strip()
            if transmatch:
                transid = f
                transid = transid.replace("transcript_id","")
                transid = transid.replace("\"","")
                transid = transid.strip()
        
        return (geneid,transid)
    
            
    @classmethod
    def modify_genome(cls,refseq,snpfile,Parent1_outout_file_genome,Parent2_outout_file_genome,genome1,genome2):
        
        genome_dict = SeqIO.to_dict(SeqIO.parse(refseq, "fasta"))
        genome_dict_mut = {}
        genome_Parent1 = {}
        genome_Parent2 = {}
        for chro, r in genome_dict.items():
            genome_dict_mut[chro] = list(r.seq)
            genome_Parent1[chro] = list(r.seq)
            genome_Parent2[chro] = list(r.seq)
        
        fh = open(snpfile,"r")
        for l in fh:
            l = l.strip()
            headerline = re.match(r'chr',l)
            if not headerline:
                #2	231	A	G	2	2
                chro,pos,base1,base2,cov1,cov2 = l.split("\t")
                #print chro,pos,base1,base2,cov1,cov2
                base1 = base1.strip()
                base2 = base2.strip()
                chro = chro.strip()
                pos = int(pos.strip())-1
                
                #print "before",genome_Parent1[chro][pos],base1,l
                genome_Parent1[chro][pos] = base1
                genome_Parent2[chro][pos] = base2
                #print "after",genome_Parent1[chro][pos],base1
        
        
        ## Write replaced sequences
        outgenome_Parent1 = [SeqRecord(Seq("".join(l)),id=str(chro)+"_"+str(genome1), name="",description="") for chro, l in genome_Parent1.items()]
        SeqIO.write(outgenome_Parent1,Parent1_outout_file_genome,"fasta")
        
        outgenome_Parent2 = [SeqRecord(Seq("".join(l)),id=chro+"_"+str(genome2), name="",description="") for chro, l in genome_Parent2.items()]
        SeqIO.write(outgenome_Parent2,Parent2_outout_file_genome,"fasta")
        
        fh.close()
    @classmethod
    def get_transcripts(cls,longest_isoform,genome_fasta,transcript_fasta,genome_strain):
        
        genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
        genome = {}
        for chro, r in genome_dict.items():
	    if re.search("_"+str(genome_strain),chro):
		chro = re.sub("_"+str(genome_strain),"",chro)
		
            genome[chro] = r.seq
       
        
        transcript = {}
        chromosome = ""
        strand = ""
	
	#if "FBgn0250158" in longest_isoform:
	#    print longest_isoform["longest_isoform"]
	#sys.exit()
        # Gene loop
        for gene in longest_isoform.keys():
            #print gene,longest_isoform[gene]
            transcript_seq1 = ""
            transcript_len1 = int(0)
            # Isoform loop
            for t in longest_isoform[gene].keys():
    
                #['4947174..4947275', '4949887..4950285', '4950551..4950611', '4950673..4951055', '4951126..4951239', '-', '4_group4']
                exons = longest_isoform[gene][t]
		#print gene,t,exons
                chromosome = exons.pop(-1)
                strand = exons.pop(-1)
		if chromosome in genome:
		    for exon in exons:
			estart,eend = exon.split("..")
			exon_len = (int(eend)-int(estart))+1
			transcript_len1 += exon_len
			exon_len=0
			exon_seq1 = genome[chromosome][int(estart)-1:int(eend)]
			transcript_seq1 += exon_seq1
                    
            # reverse complement the transcript sequence if it is on minus strand
            
            
	    if transcript_seq1!="":
		if str(strand)=="-":
		    transcript_seq1 = transcript_seq1.reverse_complement()
		    
		transcript[str(t)+"-"+str(chromosome)]=transcript_seq1
            transcript_seq1=""
            chromosome = ""
            strand = ""
            
        ## Write replaced sequences
        # Parent1
        outtranscript=[]
        for transcriptid in transcript.keys():
            
            record=SeqRecord(transcript[transcriptid],id=transcriptid, name="",description="")
            outtranscript.append(record)
        SeqIO.write(outtranscript,transcript_fasta,"fasta")
        


    @classmethod
    def simulate_reads(cls,outputdir,outputfile,inputfile,parentid,readlength,readqual,error,readtype,insert):
        
        if (readtype=="pe" or readtype=="PE" or readtype=="Pe"):
            
            #------------------
            # create read1 and read 2 fastq files
            #read1_fastq = outputdir+"/"+outputfile+"_"+str(parentid)+"_err"+str(error)+"_"+str(readtype)+"_1.fq"
            #read2_fastq = outputdir+"/"+outputfile+"_"+str(parentid)+"_err"+str(error)+"_"+str(readtype)+"_2.fq"
            
	    read1_fastq = str(outputdir)+"/"+outputfile+"_1.fq"
            read2_fastq = str(outputdir)+"/"+outputfile+"_2.fq"
	    
            ofh_A1 = open(read1_fastq,"w")
            ofh_A2 = open(read2_fastq,"w")
            
            transcript_dict = SeqIO.to_dict(SeqIO.parse(inputfile, "fasta"))
            for transid, r in transcript_dict.items():
                simulateRNAseq.simulate_reads_with_SNPs_PE(transid,r.seq,ofh_A1,ofh_A2,parentid,readlength,readqual,error,insert)
            
           # print "\n***** Simulated read fastq files for "+str(parentid)+" could be found in following loaction *****\n"
	   # print "\t"+str(read1_fastq)
	   # print "\t"+str(read2_fastq)
            
        
        elif (readtype=="se" or readtype=="SE" or readtype=="Se"):
        
            read1_fastq = str(outputdir)+"/"+str(outputfile)+".fq"
        
            ofh_A1 = open(read1_fastq,"w")
            
            transcript_dict = SeqIO.to_dict(SeqIO.parse(inputfile, "fasta"))
            for transid, r in transcript_dict.items():
                simulateRNAseq.simulate_reads_with_SNPs_SE(transid,r.seq,ofh_A1,parentid,readlength,readqual,error)
            
            #print "\n***** Simulated read fastq files for "+str(parentid)+" could be found in following loaction *****\n"
	    #print "\t"+str(read1_fastq)
            
            
        else:
            sys.stderr.write("\nERROR --> Read type is not valid; valid options are 'pe OR se'.\n\n")
            sys.exit(1)
            
    @classmethod
    def include_errors_in_reads(cls,seq, p):
            a="CGT"
            c="AGT"
            g="ACT"
            t="ACG"
            seq_new=""
#            if p==0:
#                    return seq # reduce the time for no error sequence
            for i in seq:
                    r = numpy.random.random()
                    if r < p:
                            if i == "A":
                                    seq_new+=a[numpy.random.randint(0,3)]
                            elif i == "C":
                                    seq_new+=c[numpy.random.randint(0,3)]
                            elif i == "G":
                                    seq_new+=g[numpy.random.randint(0,3)]
                            elif i == "T":
                                    seq_new+=t[numpy.random.randint(0,3)]
                            else:
                                    seq_new+=i
                    else:
                            seq_new+=i
                    
            return seq_new

    @classmethod
    def simulate_reads_with_SNPs_PE(cls,trans_id, trans_seq,ofh_A1,ofh_A2,parentid,readlength,readqual,error,insert):
    
        readlength = int(readlength)
        readqual = str(readqual)
        error = float(error)
        insert = int(insert)
        
        # simulate reads from modified transcripts
        for i,b in enumerate(trans_seq): # i starts with 0
            if (i+readlength+readlength+insert) <= len(trans_seq):
                s1=i # start pos 1st pe
                s2=i+readlength+insert #start pos 2nd pe
                
                ofh_A1.write( '@'+trans_id+"_"+str(i+1)+"_"+parentid+"#0/1\n" )
                ofh_A1.write( simulateRNAseq.include_errors_in_reads(trans_seq[s1:s1+readlength],error) +"\n" )
                ofh_A1.write( '+'+trans_id+"_"+str(i+1)+"_"+parentid+"#0/1\n" )
                ofh_A1.write( readqual * readlength +"\n" )
                
                
                ### read 2 should be reversed complemented
                read2_seq = Seq(simulateRNAseq.include_errors_in_reads(trans_seq[s2:s2+readlength],error))
                read2_seq = read2_seq.reverse_complement()

                ofh_A2.write( '@'+trans_id+"_"+str(i+1)+"_"+parentid+"#0/2\n" )
                ofh_A2.write( str(read2_seq) +"\n" )
                ofh_A2.write( '+'+trans_id+"_"+str(i+1)+"_"+parentid+"#0/2\n" )
                ofh_A2.write( readqual * readlength +"\n" )

    @classmethod
    def simulate_reads_with_SNPs_SE(cls,trans_id,trans_seq,ofh_A1,parentid,readlength,readqual,error):
        readlength = int(readlength)
        readqual = str(readqual)
        error = float(error)
        # simulate reads from modified transcripts
        for i,b in enumerate(trans_seq): # i starts with 0
            if (i+readlength) <= len(trans_seq):
                s1=i # start pos
                ofh_A1.write( '@'+trans_id+"_"+str(i+1)+"_"+parentid+"#0/1\n" )
                ofh_A1.write( simulateRNAseq.include_errors_in_reads(trans_seq[s1:s1+readlength],error) +"\n" )
                ofh_A1.write( '+'+trans_id+"_"+str(i+1)+"_"+parentid+"#0/1\n" )
                ofh_A1.write( readqual * readlength +"\n" )




