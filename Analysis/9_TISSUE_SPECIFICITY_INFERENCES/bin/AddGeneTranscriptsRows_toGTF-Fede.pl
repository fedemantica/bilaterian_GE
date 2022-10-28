#!/usr/bin/perl

$GTF=$ARGV[0];
$OUTPUT=$ARGV[1];

open (I, $GTF) || die "It cannot open the GTF $GTF file\n";
while (<I>){
    chomp;
    @t=split(/\t/);

    if ($t[2] eq "gene") {
	die "The GTF $GTF already has gene rows\n";
    }
    elsif ($t[2] eq "exon"){
	($g)=/gene_id \"(.+?)\"/;
	($tr)=/transcript_id \"(.+?)\"/;
	$tr_g{$tr}=$g;
	$i=$t[3];
	$f=$t[4];
	
	$minTR{$tr}=$i if $i <= $minTR{$tr} || !defined $minTR{$tr};
	$maxTR{$tr}=$f if $f >= $maxTR{$tr} || !defined $maxTR{$tr};
	$minG{$g}=$i if $i <= $minG{$g} || !defined $minG{$g};
	$maxG{$g}=$f if $f >= $maxG{$g} || !defined $maxG{$g};
    }
}
close I;


open (I, $GTF) || die "It cannot open the GTF $GTF file\n";
($root)=$GTF=~/(.+?)\.gtf/;
$root=~s/.+\///;
open (O, ">$OUTPUT") || die "It cannot open the output\n";
while (<I>){
    chomp;
    @t=split(/\t/);

    ($g)=/gene_id \"(.+?)\"/;
    ($tr)=/transcript_id \"(.+?)\"/;
    ($name)=/gene_name \"(.+?)\"/;
    $name="NA" if !$name;

    if (!$doneG{$g}){ #chr0protein_codingexon18158161816226.-.
#	print O "$t[0]\t$t[1]\tgene\t$minG{$g}\t$maxG{$g}\t$t[5]\t$t[6]\t$t[7]\tgene_id \"$g\"\; gene_name \"$name\"\;\n";
	$all_data{$g}="$t[0]\t$t[1]\tgene\t$minG{$g}\t$maxG{$g}\t$t[5]\t$t[6]\t$t[7]\tgene_id \"$g\"\; gene_name \"$name\"\;\n";
	$doneG{$g}=1;
    }
    if (!$doneTR{$tr}){
#	print O "$t[0]\t$t[1]\ttranscript\t$minTR{$tr}\t$maxTR{$tr}\t$t[5]\t$t[6]\t$t[7]\tgene_id \"$g\"\; transcript_id \"$tr\"\; gene_name \"$name\"\;\n";
	$all_data{$g}{$tr}.="$t[0]\t$t[1]\ttranscript\t$minTR{$tr}\t$maxTR{$tr}\t$t[5]\t$t[6]\t$t[7]\tgene_id \"$g\"\; transcript_id \"$tr\"\; gene_name \"$name\"\;\n";
	$doneTR{$tr}=1;
    }
    
#    print O "$_\n";
    $all_data{$g}{$tr}.="$_\n";
}
close I;

foreach $g (sort keys %all_data){
    print O "$all_data{$g}";
    foreach $tr (sort keys %{$all_data{$g}}){
	print O "$all_data{$g}{$tr}";
    }
}
close O;
