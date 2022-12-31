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
	my %nr;
	my %marker;
	my $ID;
	my $line = $_;
	my @outline;
	my @corrected;
	my @extra;
	chomp($line);
	$line =~ s/\t+/\t/g;
	my @f = split(/\t/,$line);
	foreach my $x (@f){
	    $nr{$x} = 1;
	}
	for (my $i = 0; $i < 10; $i++){
	    foreach my $x (keys %nr){
		my @g = split(/_/,$x);
		if ($g[0] eq $genomes[$i]){
		    if (!$ID){$ID = $x}
		    if (!$marker{$g[0]}){
			$marker{$g[0]} = 1;
			push @outline, $g[0];
		    }else{
			push @extra, $x;
		    }
		}
	    }
	}
	my $min = 10; 
	my $max = -1;
        for (my $i = 0; $i < 10; $i++){
	    if ($marker{$genomes[$i]}){
		if ($i < $min){$min = $i;}
		if ($i > $max){$max = $i;}
	    }
	}
	for (my $i = $min; $i <= $max; $i++){
	    push @corrected,$genomes[$i];
	}
	print "$ID\t",join("_",@outline),"\t",join("_",@corrected),"\n";
#	foreach my $x (@extra){
#	    print "$x\t",join("_",@outline),"\t",join("_",@corrected),"\n";
#	}
    }
    close(IN); #Read input data;
}else{
    die "$usage\n";
}
