#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

my %opt;
GetOptions("-i=s%" => \%opt);

my @genomes = ("hg38","panTro5","gorGor4","ponAbe2","nomLeu3","chlSab2","macFas5","rheMac8","papAnu2","calJac3");

my $usage= qq(
This script is used to cluster all shared active DNA transposons
Usage :  $0  [Option] bam_infiles
Option:  --i hg=(file location..)
         (input,pt,pa,gg,rm,nl)
         -h         Help Information

Author : Daniel W. Tang
Date   : Nov. 22nd, 2011
\n);

my %input;
my %revinput;
my %data;
my %output; 
my %checked;

if ($opt{input}){
    open(IN,$opt{input}) or die "can't open $opt{input}:$!\n";
    while(<IN>){
	my $line = $_;
	chomp($line);
	$line =~ s/\t+/\t/g;
	my @f = split(/\t/,$line);
	$input{$f[9]}{$f[13]} = 1;
	$input{$f[13]}{$f[9]} = 1; 
	$data{$f[9]} = 1;
	$data{$f[13]} = 1;
    }
    close(IN); #Read input data;
}else{
    die "$usage";
}

foreach my $x (keys %input){
    if (!$checked{$x}){
	my %tmp = %{$input{$x}};
	my $condition;
	do {
	    $condition = 0;
	    my $t1 = scalar(keys %tmp);
	    foreach my $y (keys %tmp){
		%tmp = map { $_ => 1 } keys %tmp, keys %{$input{$y}};
	    }
	    my $t2 = scalar(keys %tmp);
	    if ($t2 > $t1){$condition = 1}
	} until($condition == 0);    
	foreach my $y (keys %tmp){
	    $checked{$y} = 1;
	    %{$input{$y}} = %tmp;
	}
    }
}

for (my $v = 0;$v < 10; $v++){
    foreach my $x (keys %input){
	my @t = split(/_/,$x);
	my $g = shift @t;
	my $id = join("_", @t);
	if ($g eq $genomes[$v]){
	    if (!$output{$x}){
		print "$x\t";
		$output{$x} = 1;
		foreach my $y (keys %{$input{$x}}){
		    print "$y\t";
		    $output{$y} = 1;
		}
		print "\n";
	    }
	}
    }
}
