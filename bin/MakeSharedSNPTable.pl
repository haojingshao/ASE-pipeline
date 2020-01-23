#!/usr/bin/perl
use strict;
if(@ARGV!=4){
	die "perl $0 snpList nonList(depth) SharedSNP outPrefix\nperl $0 snp.list non.list SV.overlap.homo table\n";
}
my ($snpList,$nonList,$sharedSNP,$out)=@ARGV;
my $s=OpenList($snpList);
my $n=OpenList($nonList);
my $snpH=OpenHash($sharedSNP);
my $outHandleSNP;
open $outHandleSNP, ">$out.SNP.Fear2016.hisat2.table" or die $!;

#print $outHandleSVExon "feature_id\tParent1_hybrid\tParent2_hybrid\tSample\tExon_expression\tRun\tLine_sex\tSex\tLine\n";
print $outHandleSNP "chromosome\tsite\tParent1_hybrid\tParent2_hybrid\tOtherAllele\texon_id\tgene_id\tExon_expression\tGene_expression\tSample\tRun\tLine_sex\tSex\tLine\n";
for (my $counter=0;$counter<@$s;$counter++){
	CheckID($$s[$counter],$$n[$counter]);#die if wrong
	my $outSNP=ReadFile($$s[$counter],$$n[$counter],$snpH);
	print $outHandleSNP $outSNP;
}
sub OpenHash{
	my $file=shift;
	my %h;
	my $handle;
	open $handle, $file or die $!;
	while(<$handle>){
		chomp;
		my @c=split;
		$h{$c[0]}{$c[1]}="$c[5]\t$c[6]";
	}
	close $handle;
	return \%h;
}
sub ReadFile{
	my $s=shift; # actually I don't plan to read the gene file it was contained by the exon file
	my $n=shift; # file of depth
	my $snpH=shift;
	my ($exon,$gene,$region)=ReadNonFile($n); # 
	my $outSeq="";
	my ($type,$id1,$id2)=GetInfo2($s);
	my %printedExon;
	my %printedSNP;
	open my $h,$s or die $!;
	while(<$h>){
		chomp;
		my @e=split;
		if(exists $$snpH{$e[0]}{$e[1]}){
			$outSeq.="$_\t$$exon{$e[5]}\t$$gene{$e[6]}\t$id1\t$id2\n";
			#$outSeq.="$_\t$id1\t$id2\n";
			$printedExon{$e[5]}=1;
			$printedSNP{$e[0]}{$e[1]}=1;
		}
	}
	close $h;
	for my $chr (sort {$a<=>$b} keys %$snpH){
		for my $pos (sort {$a<=>$b} keys %{$$snpH{$chr}}){
			if(not exists $printedSNP{$chr}{$pos}){
				my @info=split(/\t/,$$snpH{$chr}{$pos});
				my $exonID=$info[0];
				my $geneID=$info[1];
				$outSeq.="$chr\t$pos\t0\t0\t0\t$$snpH{$chr}{$pos}\t$$exon{$exonID}\t$$gene{$geneID}\t$id1\t$id2\n";
				$printedExon{$exonID}=1;
			}
		}
	}
#	for my $exonID (sort {$a<=>$b} keys %$exon){
#		if (not exists $printedExon{$exonID}){
#			my @e=split(/\t/,$$region{$exonID});
#			$outSeq.="$e[0]\t$e[1]\t0\t0\t0\t$exonID\t$e[2]\t$$exon{$exonID}\t$$gene{$e[2]}\t$id1\t$id2\n";
#		}
#	}
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
####### the ID should be kept in the file name
sub CheckID{
	my ($g,$e)=@_;
	my $gInfo=GetInfo($g);
	my $eInfo=GetInfo($e);
	if($gInfo eq $eInfo){
		return 1;
	}else{
		die "g=$gInfo\te=$eInfo\n";
		return 0;
	}
}
sub GetInfo2{
	my $seq=shift;
	if ($seq=~/Fear2016\/(.*)\/Run_(\d+)\/line(\d+)\/([^\/]*)\//){
		my $type=$1;
		my $run=$2;
		my $line=$3;
		my $id=$4;
		$id=~/^([^\_]+)_.*$/;
		$id=$1;
		$id=~/\d+(M|F)(M|F)/;
		my $g1=$1;
		my $g2=$2;
		return $type,$id,"$run\t$g1\t$g2\t$line";
	}else{
		die "Error in getting the Run info at :$seq\n";
		return;
	}

}
sub GetInfo{
	my $seq=shift;
	if ($seq=~/Fear2016\/(.*)\/Run_(\d+)\/line(\d+)\/([^\/]*)\//){
		return "$1\t$2\t$3\t$4";
	}else{
		die "Error in getting the Run info at :$seq\n";
	}
}
sub OpenList{
	my $file=shift;
	open my $h, $file or die $!;
	my @array=<$h>;
	close $h;
	@array=sort @array;
	return \@array;
}
