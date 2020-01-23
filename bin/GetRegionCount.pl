#!/usr/bin/perl
use strict;
#####################################
# global setting
######################################
my $seqError=0.001;

######################################
## read ARGV
#####################################
if(@ARGV!=4){
	die "perl $0 inSam inRegion inHomo outPrefix\nExample:\nmodule load samtools;samtools view /home/uqhshao/30days/FruitFly/GROUP/Fear2016/Standing_Variation/Run_1/line100/100FM2_S212/final.bam |perl $0 -  non-redundant.bed /home/uqhshao/30days/FruitFly/GROUP/Fear2016/Standing_Variation/Run_1/line100/100FM2_S212/190_100 out\n";
}
my ($sam,$bed,$mask,$out)=@ARGV;
####################################################
## Read Bed =>
## BedArray keep the exon and gene ID
## BedHash keep the chr and position
####################################################
my %bed;
my @bed;
my $count=0;
open BED,$bed or die $!;
while(<BED>){
	chomp;
	if($_=~/^#/){
		next;
	}
	my @c=split;
	push @bed,"$c[3]\t$c[4]";
	for (my $i=$c[1];$i<=$c[2];$i++){
		$bed{$c[0]}{$i}=$count;
	}
	$count++;
}
close BED;
print "Reading BED is done\n";
#####################################
# Read Homo/Hetero SNP in the sample
# Asign the exon/gene name to each SNP
# Store in %hash{$chr}{$pos}="P1Allele\tP2Allele\texonID\tgeneID";
#####################################
my %hash;
open HOMO, "$mask.homo" or die $!;
while(<HOMO>){
	chomp;
	my @c=split;
	my ($chr,$pos,$pos2,$p1,$p2)=@c;
	my $region=GetGeneExonName($chr,$pos,\%bed,\@bed);
	$hash{$chr}{$pos}="$p1\t$p2\t$region\n";
}
close HOMO;
#open HETERO, "$mask.hetero" or die $!;
#while(<HETERO>){
#	chomp;
#	my @c=split;
#	my ($chr,$pos,$pos2,$p1,$p2)=@c;
#	my $region=GetGeneExonName($chr,$pos,\%bed,\@bed);
#	$hash{$chr}{$pos}="$p1\t$p2\t$region\n";
#}
#close HETERO;
print "Reading SNP is done\n";
#########################################################
# readSam
#
#
#########################################################
my $right=0;
my $miss=0;
#my %exon;#final result for exon hash-array   $exon{"exonID"}[0]=P1   [1]=P2
#my %gene;#final result for gene
my %SNP;
$count=0;
open IN, "$sam" or die $!;
while(<IN>){
	$count++;if($count%1000000==0){print STDERR "Reading sam: $count reads are read.\n"};
	chomp;
	my @c=split;
	my $chr  =$c[2];
	my $start=$c[3];
	my $cigar=$c[5];
	my $seq  =$c[9];
	my $md   ="";
	for(my $i=9;$i<@c;$i++){
		if($c[$i]=~/MD:Z:(.*)/){ # MD:Z only care about the site for reference and ignore insertion in reads!!
			$md=$1;
		}
	}
	if($md=~/N/){
		my $n=getNSite($start,$cigar,$md); #get All N in Array $n
		if($cigar=~/^(\d+)S/){ #chop the S sequence;
			my $s=$1;
			$cigar=~s/^\d+S//;
			$seq=substr($seq,$s);
		}
		for (my $j=0;$j<@$n;$j++){ # checking All the Ns in a read
			my @d=split(/\t/,$$n[$j]);
			my ($refPos,$readPos)=@d; #check wether the readPos is all correct? Answer: look OK
			if(exists $hash{$chr}{$d[0]}){
				$right++;
				#print "Right in N: $chr\t$d[0]\t$d[1]\n$_\n";
				chomp $hash{$chr}{$d[0]};
				my @e=split(/\t/,$hash{$chr}{$d[0]});
				my ($P1,$P2,$exon,$gene)=@e;
				my $siteSeq=substr($seq,$readPos,1);
				#print "$siteSeq\t$P1\t$P2\t$exon\t$gene";
				if($exon eq "None"){
					next;
				}else{
				#	Modify here
				#	$tmpExon{$exon}.="$siteSeq\t$P1\t$P2\t";
				#	$tmpGene{$gene}.="$siteSeq\t$P1\t$P2\t";
				#	Add value to $P1 or $P2
					if($siteSeq eq $P1){
						$SNP{$chr}{$d[0]}[0]++;
					}elsif($siteSeq eq $P2){
						$SNP{$chr}{$d[0]}[1]++;
					}else{
						$SNP{$chr}{$d[0]}[2]++;
					}
				}
			}else{
				$miss++;
#				print STDERR "Error in N: $chr\t$d[0]\t$d[1]\n$_\n";
			}
		}
	}
}
close IN;
print "Reading Sam is done\n";
#####################################################
# WriteOutput to table
####################################################
#WriteTable($out.".exon",\%exon,\@bed,0);
#WriteTable($out.".gene",\%gene,\@bed,1);
WriteTableSNP($out.".snp",\%SNP,\%hash);#Write the snp table

#print STDERR printError($right,$miss);
################################################
## input startSite CigarString MD:Z-String
## output \@out out[0]="outRefPos0\toutReadPos0"
################################################
sub getNSite{
	my ($start,$cigar,$md)=@_;
	#print "$start\t$cigar\t$md\n";
	my @out;
	my $refPos=$start-1;#0 base but the MD and cigar are 1 base The start will add 1 by itself
	my $readPos=-1;#0 base
	#my $cigarOffset=GetCigarOffset($cigar);#Cigar contain N that increase the reference pos
	my $numberOfOffset=0;##test whether offset is added!!!
	while($md=~/(\d+)([ATGCN]+|\^[ATGCN]+)/g){
		if(substr($2,0,1) eq "^"){
			if($2 eq "^N"){#in rare rare situation the "masked SNP" is deleted
				$readPos+=$1;
				$refPos+=$1+1;
				#skip
			}elsif(substr($2,-1,1) eq "N"){#in rare situation the deletion come with N
				$readPos+=$1+1;
				$refPos+=$1+length($2)-1-1+1;
				($readPos,$refPos,$numberOfOffset)=GetCigarOffset($cigar,$readPos,$refPos,$numberOfOffset);
				#push @out, "$refPos\t$readPos";
			}else{
				$readPos+=$1;
				$refPos+=$1+length($2)-1;
				($readPos,$refPos,$numberOfOffset)=GetCigarOffset($cigar,$readPos,$refPos,$numberOfOffset);
			}
		}elsif($2 eq "N"){
			$readPos+=$1+1;
			$refPos+=$1+1;
			($readPos,$refPos,$numberOfOffset)=GetCigarOffset($cigar,$readPos,$refPos,$numberOfOffset);
			push @out, "$refPos\t$readPos";
			#print "$refPos\t$readPos\t$numberOfOffset\n";
		}elsif($2 eq "A" or $2 eq "G" or $2 eq "C" or $2 eq "T"){
			$readPos+=$1+1;
			$refPos+=$1+1;
			($readPos,$refPos,$numberOfOffset)=GetCigarOffset($cigar,$readPos,$refPos,$numberOfOffset);
		}else{
			print STDERR "NewMD:"."$start\t$cigar\t$md\n";
		}
	}
	#print "END\n";
	return \@out;
}
#############################################
## GetTheOffset due to Cigar such as "10S20M50N30M" reduce the readlength by 10 and add the reference pos 50
##############################################
sub GetCigarOffset{
	my $cigar=shift;
	my $readPos=shift;
	my $refPos=shift;
	my $numberOfOffset=shift;
	my $currentOffset=0;
	if($cigar=~/N/){
		my $currentReadPos=-1;
		my $refOffset=0;
		my $readOffset=0;
		while($cigar=~/(\d+)([SHMIDN])/g){
		#	print "cigar=$cigar\t1=$1\t2=$2\trRefOffset=$refOffset\treadOffset=$readOffset\treadPos=$readPos\tcurrentReadPos=$currentReadPos\n";
			if($2 eq "H"){
				#skip
			#}elsif($2 eq "S"){ There is no S
			#	if($numberOfOffset == 0){ #add for the first S
			#		$readOffset+=$1; #??
			#	}
			#}
			}elsif($2 eq "M"){
				$currentReadPos+=$1;
			}elsif($2 eq "I"){
				$currentReadPos+=$1;
				$currentOffset++;
				if($currentOffset>$numberOfOffset){
					$readPos+=$1;#have to add read pos as well
				}
			}elsif($2 eq "N"){ #if there is two splicing!?
				$currentOffset++;
				if($currentOffset>$numberOfOffset){
					$refOffset+=$1;
				}
			}
			if($currentReadPos>=$readPos){
				return $readPos+$readOffset,$refPos+$refOffset,$currentOffset;
			}
		}

	}
	return $readPos,$refPos,0;
}
###################
##	check the chr pos is contained by exon or not
## if true return "exonID\tgeneID"
###################
sub GetGeneExonName{
	my ($chr,$pos,$bedH,$bedA)=@_;
	if(not exists $$bedH{$chr}{$pos}){
		return "None\tNone";
	}else{
		return "$$bedA[$$bedH{$chr}{$pos}]";
	}
	#print join "\t",@{$$bed{$chr}},"\n";
	#for (my $i=0;$i<@{$$bed{$chr}};$i++){
	#	print $i,"\t";
	#	my @c=split(/\t/,$$bed{$chr}[$i]);
	#	my ($start,$end,$exon,$gene)=@c;
	#	if( $start<=$pos and $pos<=$end){
	#		return "$exon\t$gene";
	#	}
	#}
	#return "None\tNone";
}
#####################
##	 Print the number of missing in the sam input
#####################
sub printError{
	my $right=shift;
	my $miss=shift;
	my $tmp=$right/($right+$miss);
	if($right+$miss >0){
		return "Right=$right\tMiss=$miss\tRate=$tmp\n";
	}else{
		return "The sam input is empty!\n";
	}
}
###########
## input info "seq0\tP10\tP20\tseq1\tP11\tP22\t..."
## output string P1|P2|None mean P1|P2|None is the best 
###########
sub MaximumLikelihood{
	my $info=shift;
	my $seqError=shift;
	chomp $info;
	my $logP1=0;
	my $logP2=0;
	my @c=split(/\t/,$info);
	for (my $i=0;$i<@c;$i+=3){
		if($c[$i] eq $c[$i+1]){ #seq eq P1
			$logP1+=log(1-$seqError);
			$logP2+=log($seqError);
		}elsif($c[$i] eq $c[$i+2]){#seq eq P2
			$logP1+=log($seqError);
			$logP2+=log(1-$seqError);
		}else{ #seq ne P1 ne P2
			$logP1+=log($seqError);
			$logP2+=log($seqError);
		}
	}
	my $return="";
	if ($logP1-$logP2>log(100)){
		$return="P1";
	}elsif($logP2-$logP1>log(100)){
		$return="P2";
	}else{
		$return="None";
	}
	if(0){#print for checking the ML method
		my ($out0,$out1,$out2);
		for (my $i=0;$i<@c;$i+=3){
			$out0.=$c[$i];
			$out1.=$c[$i+1];
			$out2.=$c[$i+2];
		}
		print "$out0\t$out1\t$out2\t$logP1\t$logP2\t$return\n";
	}
	return $return;
}
########################
### write hash to file
########################
sub WriteTable{
	my ($file,$hash,$array,$pos)=@_;
	my $h;
#	my @keys=keys @{$$hash};
#	print join "\t",@keys,"\n";
	my $oldName="";
	open $h,">$file" or die $!;
	for (my $i=0;$i<@$array;$i++){
		my @c=split(/\t/,$$array[$i]);
		my $name=$c[$pos];
		if($oldName eq $name){
			next;
		}
		if(exists $$hash{$name}){
			if(not exists $$hash{$name}[0] ){$$hash{$name}[0]=0;}
			if(not exists $$hash{$name}[1] ){$$hash{$name}[1]=0;}
			print $h "$name\t$$hash{$name}[0]\t$$hash{$name}[1]\n";
		}else{
			print $h "$name\t0\t0\n";
		}
		$oldName=$name;
	}
	close $h;
}
#WriteTableSNP($out.".snp",\%SNP,\%hash);#Write the snp table
sub WriteTableSNP{
	my ($out,$SNP,$hash)=@_;
	my $h;
	open $h, ">$out" or die $!;
	for my $chr (sort {$a<=>$b} keys %$SNP){
		for my $pos(sort {$a<=>$b} keys %{$$SNP{$chr}}){
			if(not exists $$SNP{$chr}{$pos}[0]){$$SNP{$chr}{$pos}[0]=0};
			if(not exists $$SNP{$chr}{$pos}[1]){$$SNP{$chr}{$pos}[1]=0};
			if(not exists $$SNP{$chr}{$pos}[2]){$$SNP{$chr}{$pos}[2]=0};
			chomp $$hash{$chr}{$pos};
			my @e=split(/\t/,$$hash{$chr}{$pos});
			my ($P1,$P2,$exon,$gene)=@e;
			print $h "$chr\t$pos\t$$SNP{$chr}{$pos}[0]\t$$SNP{$chr}{$pos}[1]\t$$SNP{$chr}{$pos}[2]\t$exon\t$gene\n";
		}
	}
	close $h;
}
