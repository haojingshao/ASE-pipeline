#!/usr/bin/perl
use strict;
if(@ARGV!=4){
    die "perl $0 inref inVCF1 inVCF2 outPrefix> STDOUT\noutput: prefix.masked.fa --masked reference for both homo and hetero sites\noutput: prefix.homo --homo info in \"bed\" format(chr start end VCF1 VCF2)\noutput: prefix.hetero --hetero info in \"bed\" format(chr start end VCF1 VCF2)\n";
}
my $refH=readRef($ARGV[0]);
print STDERR "Read Ref done\n";
my %empty;
my $h=ReadVCF($ARGV[1],\%empty,1,$refH);
print STDERR "Read VCF1 done\n";
$h=ReadVCF($ARGV[2],$h,2,$refH);
print STDERR "Read VCF2 done\n";
MakeReference($h,$refH,$ARGV[3].".masked.fa",$ARGV[3].".homo",$ARGV[3].".hetero");
sub MakeReference{
	my ($h,$refH,$outRef,$outHomo,$outHetero)=@_;
	open OUT, ">$outRef" or die $!;
	open HOMO, ">$outHomo" or die $!;
	open HETE, ">$outHetero" or die $!;
	for my $chr (sort {$a<=>$b} keys %$h){
		for my $pos (sort {$a<=>$b} keys %{$$h{$chr}}){
			chomp $$h{$chr}{$pos};
			my @c=split (/\ /,$$h{$chr}{$pos});
	#		print scalar @c,"\t$$h{$chr}{$pos}\n";
			if(scalar @c eq 1){ #one vcf is reference
				if(length $c[0] eq 2){ # one HOMO SNP
					my $refSeq=substr($$refH{$chr},$pos,1);
					substr($$refH{$chr},$pos,1)="N";
					$c[0]=~/^(\d)(.)$/;
					if($1 eq "1"){
						print HOMO "$chr\t$pos\t$pos\t$2\t$refSeq\n";
					}else{
						print HOMO "$chr\t$pos\t$pos\t$refSeq\t$2\n";
					}
				}else{ #one HETERO SNP
					my $refSeq=substr($$refH{$chr},$pos,1);
					substr($$refH{$chr},$pos,1)="N";
					$c[0]=~/^(\d)(.+)$/;
					if($1 eq "1"){
						print HETE "$chr\t$pos\t$pos\t$2\t$refSeq\n";
					}else{
						print HETE "$chr\t$pos\t$pos\t$refSeq\t$2\n";
					}
				}
			}else{ #both vcf is not reference
				if((length $c[0] eq 2) and (length $c[1] eq 2)){ # two HOMO SNP
					my $refSeq=substr($$refH{$chr},$pos,1);
					$c[0]=~/^(\d)(.+)$/;
					my $snp1=$2;
					$c[1]=~/^(\d)(.+)$/;
					my $snp2=$2;
					if($snp1 eq $snp2){
						substr($$refH{$chr},$pos,1)=$snp1;
					}else{
						print HOMO "$chr\t$pos\t$pos\t$snp1\t$snp2\n";
						substr($$refH{$chr},$pos,1)="N";
					}
				}else{ #two HETERO SNP
					my $refSeq=substr($$refH{$chr},$pos,1);
					substr($$refH{$chr},$pos,1)="N";
					$c[0]=~/^(\d)(.+)$/;
					my $snp1=$2;
					$c[1]=~/^(\d)(.+)$/;
					my $snp2=$2;
					print HETE "$chr\t$pos\t$pos\t$snp1\t$snp2\n";
				}

			}
		}
	}
	#print new reference
	for my $chr (sort {$a<=>$b} keys %$refH){
		print OUT ">$chr\n";
		for (my $pos=1;$pos<length $$refH{$chr};$pos+=60){
			print OUT substr($$refH{$chr},$pos,60),"\n";
		}
	}
	close OUT;
	close HOMO;
	close HETE;
}

sub ReadVCF{
	my $vcf=shift;
	my $hash=shift;
	my $id=shift;
	my $refH=shift;
	my $count=0;
	open IN, $vcf or die $!;
	while(<IN>){
		$count++;if($count%10000000==0){print STDERR "VCF $id ";print STDERR $count," lines are done.\n"};
		chomp;
		if($_=~/^#/){
			next;
		}
		my @c=split;
	#		print $c[3],"\t",substr($$refH{$c[0]},$c[1],1),"\n";
		if($c[-1]=~/0\/1/){#heterozygous
			$$hash{$c[0]}{$c[1]}.="$id$c[3]/$c[4] ";
		#	print "hetero:$c[0]\t$c[1]\t$id$c[3]$c[4]\n";
		}elsif($c[-1]=~/1\/1/){
		#	print $c[3],"\t",substr($$refH{$c[0]},$c[1],1),"\n";
			$$hash{$c[0]}{$c[1]}.="$id$c[4] ";
		}
	}
	close IN;
	return $hash;
}

sub readRef{
    my $file=shift;
    my %h;
    my $refChr;
    open FILE,"$file" or die $!;
        while(<FILE>){
        chomp;
                if($_=~/^>([^ ]*)/){
                        $h{$1}="N";#array0 =N
                        $refChr=$1;
                }else{
                        $h{$refChr}.=uc ($_);
                }
        }
        close FILE;
    return \%h;
}
