#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use Cwd;
my ($vcf1,$vcf2,$ref,$fastq,$vcf1Name,$vcf2Name);
my $usage="Usage: $0 --vcf1 --vcf2 --ref --fastq --outDir --path --exon
Options:
--vcf1 STR          Parent1 genome vcf file
--vcf2 STR          Parent2 genome vcf file
--ref STR           Genome reference file
--fastq STR         A file contains all the fastqs.
                    Header format: fastq1\tfastq2\tInfo1\tInfo2\tInfo3...
                    Body format:   fq1\tfq2\tInfo1Value\tInfo2Value\tInfo3Value...
                    Note1: Info[1-n] and their values will appear in the final table header and body, respectively.
                    Note2: The combination of 'Info1-Info2-...-Info[n]' should be unique.
--path STR          Installation directory. 'InstallationPath'
--exon STR          A non-overlap file contains all the exon and gene infomation.
                    Format: chromosome\tstart\tend\texonName\tgeneName
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
Example: 
$0 --vcf1=TestData/P1.vcf --vcf2=TestData/P2.vcf --ref=TestData/ref.fa --fastq=TestData/list --vcf1Name=190 --vcf2Name=226 --outDir=TestOut --cpu=8 --path=/home/uqhshao/30days/FruitFly/pipeline/github/bin/ --exon=TestData/exon.bed --hisat2=/sw/QFAB/installs/hisat2/hisat2-2.1.0/hisat2 --hisat2-build=/sw/QFAB/installs/hisat2/hisat2-2.1.0/hisat2-build --samtools=/opt/biotools/samtools/1.3/bin/samtools --bedtools=/opt/biotools/bedtools/bin/bedtools\n";

our %options;
GetOptions( \%options, 'path=s','vcf1=s' ,'vcf2=s','ref=s','fastq=s','vcf1Name:s','vcf2Name:s','outDir=s','cpu:i','exon=s','hisat2:s','hisat2-build:s','samtools:s','bedtools:s','mode:s','readlength:i','insertsize:i');
if($options{"vcf1Name"} eq ""){$options{"vcf1Name"}="vcf1";};
if($options{"vcf2Name"} eq ""){$options{"vcf2Name"}="vcf2";};
if($options{"hisat2"} eq ""){$options{"hisat2"}="hisat2";};
if($options{"hisat2-build"} eq ""){$options{"hisat2-build"}="hisat2-build";};
if($options{"samtools"} eq ""){$options{"samtools"}="samtools";};
if($options{"bedtools"} eq ""){$options{"bedtools"}="bedtools";};
if($options{"cpu"} ==0){$options{"cpu"}=1;};
if($options{"readlength"} ==0){$options{"readlength"}=151;};
if($options{"insertsize"} ==0){$options{"insertsize"}=350;};
if($options{"vcf1"} eq "" ){die "Missing vcf1\n$usage"};
if($options{"vcf2"} eq "" ){die "Missing vcf2\n$usage"};
if($options{"ref"} eq "" ){die "Missing ref\n$usage"};
if($options{"ref"}!~/\.fasta$/ and $options{'ref'}!~/\.fa$/){die "Ref must end with \".fasta\" or \".fa\" \n$usage"};
if($options{"fastq"} eq "" ){die "Missing fastq\n$usage"};
if($options{'outDir'} eq "" ){die "Missing outDir\n$usage"};
if($options{"exon"} eq "" ){die "Missing exon bed\n$usage"};
our $mode="run";
if($options{"mode"} eq "script" ){$mode="script"};

#my  = getcwd; # path for the data
#="/";
#my $out=FindFullPath($options{'outDir'},);
#print $out,"\n";
#die;
my $info=ReadList($options{"fastq"});
my $cmd="mkdir -p $options{'outDir'}";
Cmd($cmd,1);
$cmd="perl $options{'path'}/bin/CreatMaskedReference.pl $options{'ref'} $options{'vcf1'} $options{'vcf2'} $options{'outDir'}/$options{'vcf1Name'}_$options{'vcf2Name'}";
DeleteOldScripts();
Cmd($cmd,1);
if($mode eq "run"){print "#Making masked referece\n";}
my $ref="$options{'outDir'}/$options{'vcf1Name'}_$options{'vcf2Name'}.masked.fa";
$cmd="$options{'hisat2-build'} $ref $ref";
if($mode eq "run"){print "#Indexing masked referece\n";}
Cmd($cmd,1);
$cmd="";
open my $simulateh,">$options{'outDir'}/qsub-simulate.sh" or die $!;
print $simulateh "#!/bin/sh\n";
print $simulateh "mkdir -p $options{'outDir'}/simulation ;cd $options{'outDir'}/simulation\n";
print $simulateh "$options{'samtools'} faidx $options{'ref'}\n";
my $dict=$options{'ref'}; 
$dict=~s/\.fasta/.dict/g;
$dict=~s/\.fa/.dict/g;
if(not -e $dict){
	print $simulateh "java -jar $options{'path'}/bin/picard.jar CreateSequenceDictionary R=$options{'ref'} O=$dict\n";
}
print $simulateh "java -jar $options{'path'}/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R $options{'ref'} -V $options{'vcf1'} -o $options{'vcf1Name'}.fasta ; perl -i -p -e 's/>\\d+ />/;s/:1//' $options{'vcf1Name'}.fasta\n";
print $simulateh "java -jar $options{'path'}/bin/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R $options{'ref'} -V $options{'vcf2'} -o $options{'vcf2Name'}.fasta ; perl -i -p -e 's/>\\d+ />/;s/:1//' $options{'vcf2Name'}.fasta\n";
print $simulateh "python $options{'path'}/bin/Simulation.py ./ $options{'vcf1Name'} $options{'vcf2Name'} $options{'vcf1Name'}.fasta $options{'vcf2Name'}.fasta $options{'readlength'} $options{'insertsize'}\n";
print $simulateh "cd ../..\n".GenerateAlignmentShell("$options{'outDir'}/simulation/$options{'vcf1Name'}_1.fq","$options{'outDir'}/simulation/$options{'vcf1Name'}_2.fq",$options{'vcf1Name'},"$ref");
print $simulateh "cd ../..\n".GenerateAlignmentShell("$options{'outDir'}/simulation/$options{'vcf2Name'}_1.fq","$options{'outDir'}/simulation/$options{'vcf2Name'}_2.fq",$options{'vcf2Name'},"$ref");
print $simulateh "cd ../\necho \"chromosome site Parent1Allele1 Parent1Allele2 Parent2Allele1 Parent2Allele2 exonID geneID\" > simulation.table\n";
print $simulateh "paste $options{'vcf1Name'}/Fear2016.snp $options{'vcf1Name'}/Fear2016.snp |awk \'{print \$1,\$2,\$3,\$4,\$10,\$11,\$6,\$7}\' >> simulation.table\n";
close $simulateh;
$cmd.="sh $options{'outDir'}/qsub-simulate.sh\n";
if($mode eq "run"){print "#Generating shell for each alignment. You may qsub each shell.\n";}
open my $snph, 	">snp.list"  or die $!;
open my $depthh, ">depth.list" or die $!;
for (my $i=1;$i<scalar @$info;$i++){
	my $id="";
	for (my $j=2;$j< scalar @{$$info[$i]};$j++){
		$id.="$$info[$i][$j]-";
	}
	chop $id;
	open my $outh,">$options{'outDir'}/qsub-$id.sh" or die $!;
	$cmd.="sh $options{'outDir'}/qsub-$id.sh\n";
	print $outh GenerateAlignmentShell($$info[$i][0],$$info[$i][1],$id,"$ref");
	close $outh;
	print $snph  "$options{'outDir'}/$id/Fear2016.snp\n";
	print $depthh "$options{'outDir'}/$id/exon.table\n";
}
Cmd($cmd,2);
close $snph;
close $depthh;
$cmd="perl $options{'path'}/bin/MakeSNPTable.pl snp.list depth.list $options{'outDir'}/ASE.table $options{'fastq'}\n";
Cmd($cmd,3);
my $len=(scalar @$info)-1;
$cmd="perl $options{'path'}/bin/FilterExpressExon.pl $options{'outDir'}/ASE.table $len > $options{'outDir'}/ASE.table.expressedExon 2> $options{'outDir'}/ASE.table.stat \n";
Cmd($cmd,3);
if($mode eq "run"){print "End of the program.\n";}
else{print "Scripts (step*sh) are generated.\n"}
exit 1;

#input 	fastq list
#output fastq array
#header	= T
sub GenerateAlignmentShell{
	my $fq1=shift;
	my $fq2=shift;
	my $id=shift;
	my $ref=shift;
	my $prefix=$ref;
	$prefix=~s/\.masked\.fa$//g;
	my $out=
"#!/bin/bash
dir=$options{'outDir'}/$id
mkdir -p \$dir
cd \$dir
fq1=$fq1
fq2=$fq2
ref=$ref
$options{'hisat2'} -q -x \$ref -1 \$fq1 -2 \$fq2  -S RNA_hisat.sam -p $options{'cpu'}
$options{'samtools'} view -S -bh -T \$ref  RNA_hisat.sam  > RNA_hisat.bam
rm -rf RNA_hisat.sam
$options{'samtools'} sort -m 10G -@ 3 -o RNA_hisat.sorted.bam RNA_hisat.bam
rm -rf RNA_hisat.bam ref.fa*
$options{'samtools'} view -bh -F 256 -f 2 -q 30 RNA_hisat.sorted.bam > final.bam
$options{'samtools'} index final.bam
$options{'bedtools'} intersect -c -a $options{'exon'} -b final.bam > exon.table
#exon.tab is overall expressed depth for each exon including the non-SNP expressed exon.
$options{'samtools'} view final.bam |perl $options{'path'}/bin/GetRegionCount.pl - $options{'exon'} $prefix  Fear2016
#Fear2016.snp is homozygous SNP read counts.
\n";
	return $out;
}

sub ReadList{
	my $list=shift;
	my @all;
	my $i=0;
	open my $h, $list or die $!;
	while(<$h>){
		chomp;
		@{$all[$i]}=split;
		$i++;
	}
	close $list;
	return \@all;
}

sub Cmd{
	my $cmd=shift;
	my $step=shift;
	if($mode eq "run"){
		print "#Running cmd: \n$cmd\n";
		system($cmd);
	}else{
		my $h;
		open $h, ">>step$step.sh" or die $!;
		print $h $cmd."\n";
		close $h;
	}
}
sub DeleteOldScripts{
	if($mode eq "script" and (-e "step1.sh" or -e "step2.sh" or -e "step3.sh")){
		system("rm step1.sh step2.sh step3.sh");
		print "Delete old scripts.\n";
	}
	return;
}
