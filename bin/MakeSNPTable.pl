#!/usr/bin/perl
use strict;
if(@ARGV<3 or @ARGV>4){
	die "perl $0 snpList nonList(depth) outPrefix append\nperl $0 snp.list non.list table append\n";
}
#my ($snpList,$nonList,$out)=@ARGV;
my $snpList=shift;
my $nonList=shift;
my $out=shift;
my $append=shift;
my $appendArray=ReadAppend($append);
my $s=OpenList($snpList);
my $n=OpenList($nonList);
my $outHandleSNP;
open $outHandleSNP, ">$out.SNP.Fear2016.hisat2.table" or die $!;

#print $outHandleSVExon "feature_id\tParent1_hybrid\tParent2_hybrid\tSample\tExon_expression\tRun\tLine_sex\tSex\tLine\n";
print $outHandleSNP "chromosome\tsite\tParent1_hybrid\tParent2_hybrid\tOtherAllele\texon_id\tgene_id\tExon_expression\tGene_expression\t$$appendArray[0]\n";
for (my $counter=0;$counter<@$s;$counter++){
#	CheckID($$s[$counter],$$n[$counter]);#don't check
	my $outSNP=ReadFile($$s[$counter],$$n[$counter],$$appendArray[$counter+1]);
	print $outHandleSNP $outSNP;
}
sub ReadFile{
	my $s=shift; # actually I don't plan to read the gene file it was contained by the exon file
	my $n=shift; # file of depth
	my $append=shift;
	my ($exon,$gene,$region)=ReadNonFile($n); # 
	my $outSeq="";
	#my ($type,$id1,$id2)=GetInfo2($s);
	my %printedExon;
	open my $h,$s or die $!;
	while(<$h>){
		chomp;
		my @e=split;
		#$outSeq.="$_\t$$exon{$e[5]}\t$$gene{$e[6]}\t$id1\t$id2\n";
		$outSeq.="$_\t$$exon{$e[5]}\t$$gene{$e[6]}\t$append\n";
		$printedExon{$e[5]}=1;
	}
	close $h;
	for my $exonID (sort {$a<=>$b} keys %$exon){
		if (not exists $printedExon{$exonID}){
			my @e=split(/\t/,$$region{$exonID});
			$outSeq.="$e[0]\t$e[1]\t0\t0\t0\t$exonID\t$e[2]\t$$exon{$exonID}\t$$gene{$e[2]}\t$append\n";
		}
	}
	return $outSeq;
}
#############
## read a non file
## return 3 hash
## %exon $exon{$exonID}=$depth
## %gene $gene{$geneID}=$geneDepth
## %region $region{$exonID}="$chr\t$pos\t$geneID"; # for this exon
sub ReadNonFile{
	my $file=shift;
	my %exon;
	my %gene;
	my %region;
	my $oldID="";
	my $gExpression1=0;
	my $gExpression2=0;
	my $gDepth=0;
	open my $h, $file or die $!;
	while(<$h>){
		chomp;
		my @r=split(/\t/,$_);
		$region{$r[3]}="$r[0]\t$r[1]\t$r[4]";
		$exon{$r[3]}=$r[5];
		if ($oldID ne $r[4]){
			if ($oldID ne ""){
				$gene{$oldID}=$gDepth;
			}
			$gDepth=$r[5];
		}else{
			$gDepth+=$r[5];
		}
		$oldID=$r[4];
	}
	$gene{$oldID}=$gDepth;
	close $h;
	return \%exon,\%gene,\%region;
}
sub ReadFiles{
	my ($g,$e,$r)=@_; # actually I don't plan to read the gene file it was contained by the exon file
	my $outSVExon;
	my $outSVGene;
	my $outMAExon;
	my $outMAGene;
	my ($type,$id1,$id2)=GetInfo2($g);
	my ($seq1,$seq2)=ReadTwoFile($e,$r,$id1,$id2);
	if($type eq "Standing_Variation"){
		$outSVExon.=$seq1;
		$outSVGene.=$seq2;
	}elsif($type eq "Mutation_Accumulation"){
		$outMAExon.=$seq1;
		$outMAGene.=$seq2;
	}else{
		die "Error in type:$type\n";
	}
	return $outSVExon,$outSVGene,$outMAExon,$outMAGene;
}
sub OpenList{
	my $file=shift;
	open my $h, $file or die $!;
	my @array=<$h>;
	close $h;
	@array=sort @array;
	return \@array;
}
sub ReadAppend{
	my $file=shift;
	my $i=0;
	my @out;
	open my $h, $file or die $!;
	while(<$h>){
		my @d=split;
		for (my $j=2;$j<scalar @d;$j++){
			$out[$i].=$d[$j]."\t";
		}
		chop $out[$i];
		$i++;
	}
	close $h;
	return \@out;
}
