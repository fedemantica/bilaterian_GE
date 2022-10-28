#!/usr/bin/perl
#use warnings;
#use strict;

$gde=$ARGV[0]; #alignment from mafft
#$output=$ARGV[1]; ##OUTPUT FILE

#`perl /users/mirimia/ymarquez/SOFTWARE/AlignIntronPos.pl $fasta_file`;
#@tmp=split(/\./,$fasta_file);
$id=0;
#my $gde=$tmp[0].".gde";

open (ALN, "$gde");
$rs1="";
$rs2="";
while (<ALN>){
	chomp($_);
	if ($_=~/\%/){
		$_=~s/\%//;
		if (!$n1){									
			$n1=$_;	$id=1;
		} else {	
			$n2=$_; $id=2; 	
		}
	}
	elsif ($id==1){
		$rs1.=$_;
	}
	elsif ($id==2){
		$rs2.=$_;
	}
}

$sim_score=score_proteins($n1,$n2,$rs1,$rs2);			
@scores=split(/\,/,$sim_score);
$l1=length($rs1);
$l2=length($rs2);
#open (OUT, ">$output");
#print OUT "Query\tSubject\tSim_score\n";
print "$n1\t$n2\t$scores[0]\n";
print "$n2\t$n1\t$scores[1]\n";


sub score_proteins{
	my %sim; ##AA similarity

	$sim{"FY"}=1;
	$sim{"YF"}=1;
	$sim{"FW"}=1;
	$sim{"WF"}=1;
	$sim{"YW"}=1;
	$sim{"WY"}=1;
	$sim{"VI"}=1;
	$sim{"IV"}=1;
	$sim{"VL"}=1;
	$sim{"LV"}=1;
	$sim{"LI"}=1;
	$sim{"IL"}=1;
	$sim{"RK"}=1;
	$sim{"RH"}=1;
	$sim{"KR"}=1;
	$sim{"HR"}=1;
	$sim{"KH"}=1;
	$sim{"HK"}=1;
	$sim{"DE"}=1;
	$sim{"ED"}=1;
	$sim{"ST"}=1;
	$sim{"TS"}=1;
	$sim{"NQ"}=1;
	$sim{"QN"}=1;

	my ($n1,$n2,$seq1,$seq2)=@_;
	#print "$n1\t$n2\t$seq1\t$seq2\n";
	my (@s1,@s2);
	@s1=split //,$seq1;
	@s2=split //,$seq2;	
	my $ng=-1;
	my $sim_score=0;
	my $g=0;
	my ($l1,$l2);
	my $t1=$seq1;
	$t1=~s/\-//g;
	my $t2=$seq2;
	$t2=~s/\-//g;
	$l1=length($t1); 
	$l2=length($t2); 
	my ($sim1)=(0,0,0,0);
	my $gap=0;
	my $i;

	##scoring proteins

	for ($i=0; $i<scalar(@s1); $i++){

		if ($sim{$s1[$i].$s2[$i]} || $s1[$i] eq $s2[$i] ){
			$sim_score++;	
			$ng=0;
		}			

		elsif ($s2[$i] eq "-" || $s1[$i] eq "-") { 
			if ($ng==0) { $g++; $ng=1; } 
			else { $e++; $ng=1; } 		
		}

	}
	#print "##$sim_score\t$l1\n";
	$sim1=0; 
	if ($l1>0 && $l2>0) {
		#Only considering the length of the protein of species 1
		$sim1=sprintf("%.6f",(($sim_score/$l1)));
		$sim2=sprintf("%.6f",(($sim_score/$l2)));
		#RETURNING RESULTS OF SIMILARITY P1 VS P2
		$tsim=$sim1.",".$sim2;
		return ($tsim);
		
	}

}
