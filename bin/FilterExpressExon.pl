#!/usr/bin/perl
use strict;
if(@ARGV!=2){
	die "perl $0 inTable NumberOfSharedSample(320) > outTable\nfilter the exon that don't have SNP in all sample\n";
}
my $h;
my %exon;
my ($file,$numberOfSample)=@ARGV;
open $h, $ARGV[0] or die $!;
while(<$h>){
	chomp;
	my @c=split;
	if($c[2]>0 or $c[3]>0){
		my $id="";
		for (my $i=10;$i<@c;$i++){
			$id.="$c[$i]-";
		}
		chop $id;
		$exon{$c[5]}{$id}=1;
	}
}
close $h;
print STDERR "ExonID\tNumberOfExpressedSample\n";
for my $exonID (keys %exon){
	print STDERR $exonID,"\t",scalar(keys %{$exon{$exonID}}),"\n";
#	if ((scalar (keys %{$exon{$exonID}})) == $numberOfSample){
#		print STDERR $exonID,"\n";
#	}else{
#		print STDERR $exonID,"\t",scalar(keys %{$exon{$exonID}}),"\n";
#	}
}
#close;

open $h, $ARGV[0] or die $!;
$_=<$h>;
print $_;
while(<$h>){
	chomp;
	my @c=split;
	if(exists $exon{$c[5]} and scalar (keys %{$exon{$c[5]}}) == $numberOfSample){
		print $_,"\n";
	}else{
#		$c[2]=0;
#		$c[3]=0;
#		$c[4]=0;
#		print join "\t",@c;
#		print "\n";
	}
}
close $h;



