#!/usr/local/bin/perl -w
#Tang_find_typical_RE: identify if a typical RE insertion is human-speicific by compare it to the chimp genome

use strict;
use warnings;
use Getopt::Std;

my %POS;
sub sub1{
    my $line1 =shift;
    $line1 = reverse($line1);
    $line1 =~ tr/ATCGatcg/TAGCtagc/;
    return($line1);
};


my %opt;
getopts("I:J:h", \%opt);
my $usage= qq(
This script is used to identify whether a typical RE insertion is human-speicific by comparing it to the chimp genome
Usage :  $0  [Option] bam_infiles
Option:  -I         input file name
         -J         sequences file folder name
         -h         Help Information

Author : Daniel W. Tang
Date   : Mar. 1st, 2011
\n);

die "$usage" if (!$opt{I}||!$opt{J}||$opt{h});
my @input;
my %chr;
open(IN,"$opt{I}") or die "can't open $opt{I}:$!\n";
while (<IN>){
    my $line = $_;
    chomp($line);
    push(@input,$line);
#    print $line;
#    exit;
}
close(IN);#read input files and store the data;
foreach my $x (@input){
    my @f = split(/\t/,$x);
    my $c = $f[1];
    if (!$chr{$c}){
	open(REF,"$opt{J}/$c.fa") or die "can't open $opt{J}/$c.fa:$!\n$x\n$f[1]\n";
	while (<REF>){
	    my $line = $_;
	    chomp($line);
	    next if ($line =~ />/);
	    $chr{$c} .= $line;
	}
    }
    my ($id) = $f[0] =~ /(^\d+)/;
    $f[0] = $id;
#    print ">$id","*i\n";
#    print substr($chr{$f[1]},$f[2] - 101,100),"\n";
#    print substr($chr{$f[1]},$f[3],100),"\n";
#    print ">$id","*5\n";
#    print substr($chr{$f[1]},$f[2] - 101,100),"\n";
#    print substr($chr{$f[1]},$f[2] - 1,100),"\n";
#    print ">$id","*3\n";
#    print substr($chr{$f[1]},$f[3] - 100,100),"\n";
#    print substr($chr{$f[1]},$f[3],100),"\n";
    $f[10] =~ tr/\/_/,/; $f[10] =~ s/,+/,/g; $f[11]=~ tr/\/_/,/; $f[11] =~ s/,+/,/g; 
    my @s=split(/,/,$f[10]);
    my @e=split(/,/,$f[11]);
    my @ss = sort { $a <=> $b } @s;
    my @se = sort { $a <=> $b } @e;
    if (!(@s ~~ @ss && @e ~~ @se)){
#	print "$x\n";
	next;
    }
    my $outseq = "";
    for(my $i=0;$i<scalar(@s);$i++){
	if ($i <scalar(@s)-1 && $e[$i] >= $s[$i+1]){
	    $e[$i] = $s[$i+1] - 1;
	}
	$outseq .= substr($chr{$c},$s[$i]-1,$e[$i]-$s[$i]+1);
#	print "$s[$i]-$e[$i]\n";
    }
    my $tmplen = length($outseq);
    my $ID;
    if ($f[6] eq "LTR"){$ID = $f[0]."|".$f[7]."|".$f[6]."|".$f[9];}else{$ID = $f[0]."|".$f[5]."|".$f[6]."|".$f[9]}   
    print ">$ID\n$outseq\n";
}
