# A pipeline for Detecting Allele Specific Expression
Data:   RNA-seq data
 
Source: [Allim](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3739924/) and [Feat et.al 2016](https://www.genetics.org/content/203/3/1177).
 
Goal:	This pipeline is designed to detect allele specific expression difference in population RNA-seq data.

### Installation
 
Download the files.

### Requirement
bamtools
 
samtools
 
bedtools
 
hisat2
 
biopython

### Pipeline workflow

1 Extract all the homologous single nucleotide mismatches from two parents

2 Replace all these sites into 'N' and generate a masked reference

3 Estimate mapping bias by simulation.

	3.1 Generate reference sequence for two parents by replacing mismatches
	
	3.2 Simulate non-error short reads from two parents
	
	3.3 Align two parents simulated reads to the masked reference
	
	3.4 Count the number of read for each single nucleotide mismatch

4 Estimate the expression for each exon mismatch

	3.1 Align all the raw data to the masked reference
	
	3.2 Count the number of read for each single nucleotide mismatch

### Prepare data

reference genome

vcf for parent1 and parent2

### Usage
```
Usage: ./pipeline.pl --vcf1 --vcf2 --ref --fastq --outDir --bin --exon
Options:
--vcf1 STR          Parent1 genome vcf file
--vcf2 STR          Parent2 genome vcf file
--ref STR           Genome reference file
--fastq STR         A file contains all the fastqs.
                    Header format: fastq1	fastq2	Info1	Info2	Info3...
                    Body format:   fq1	fq2	Info1Value	Info2Value	Info3Value...
                    Note1: Info[1-n] and their values will appear in the final table header and body, respectively.
                    Note2: The combination of 'Info1-Info2-...-Info[n]' should be unique.
--bin STR           Bin directory. 'InstallationPath/bin'
--exon STR          A non-overlap file contains all the exon and gene infomation.
                    Format: chromosome	start	end	exonName	geneName
--mode (run|script) run: run the pipeline directly [run]
                    script: generate script instead of running (step1.sh step2[id1,id2,...].sh step3.sh)
--vcf1Name STR      Name for vcf1. [190]
--vcf2Name STR      Name for vcf2. [226]
--outDir STR        Output directory
--cpu INT           Number of cpu. [1]
--readlength INT    Maximum read length for simulation. [151]
--insertsize INT    Read insert size for simulation.    [350]
--hisat2 STR        Path for hisat2 [hisat2]
--hisat2-build STR  Path for hisat2-build [histat2-build]
--samtools STR      Path for samtools [samtools]
--bedtools STR      Path for bedtools [bedtools]
```

### Contact
Email: uqhshao at uq.edu.au or haojingshao at gmail.com

### Citation
Manuscript in preparation.
