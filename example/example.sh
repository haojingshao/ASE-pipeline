#!/bin/sh
module load biopython
module load samtools
module load bedtools
module load hisat2
dir=`pwd`
echo "fastq1 fastq2 Run line id" > list
echo "$dir/1.1.fq $dir/1.2.fq Run_1 33 1" >> list
echo "$dir/2.1.fq $dir/2.2.fq Run_1 33 2" >> list
../pipeline.pl --vcf1=$dir/P1.vcf --vcf2=$dir/P2.vcf --ref=$dir/ref.fa --fastq=$dir/list --vcf1Name=P1 --vcf2Name=P2 --outDir=$dir/TestOut --cpu=8 --path=$dir/../ --exon=$dir/exon.bed
