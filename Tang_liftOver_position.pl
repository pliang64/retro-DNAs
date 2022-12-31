#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

my @genomes = ("hg38","panTro5","gorGor4","ponAbe2","nomLeu3","chlSab2","macFas5","rheMac8","papAnu2","calJac3");

my %opt;
GetOptions("-i=s%" => \%opt);

my $usage= qq(
This script is used to compare REs position in human genome with gap information 
Usage :  $0  [Option] bam_infiles
Option:  --i hg38=(liftOver file location..) 
             rawlist=(rawlist file location);
             output=(output file location);

Author : Daniel W. Tang
Date   : May. 4th, 2014
\n);

my %POS;
my %input;
my %data5;
my %data3;
my $count = 0;

my @filename = ("rawlist","output");

if (!%opt){die "$usage"};
open(IN,"$opt{$filename[0]}") or die "can't open $opt{$filename[0]}:$!\n";
while(<IN>){
    my $line = $_;
    chomp($line);
    $line =~ s/\t+/\t/g;
    my @f = split(/\t/,$line);
    my ($id) = $f[0];
    $input{$id} = $line;
}
close(IN); #Read input data;

for (my $v = 0;$v < 10; $v++){
    if ($opt{$genomes[$v]}){
	my %t5;
	my %t3;
#	print "$genomes[$v]:$opt{$genomes[$v]}\n";
	open(LO,"$opt{$genomes[$v]}") or die "can't open $opt{$genomes[$v]}}:$!\n";
	while(<LO>){
	    my $line = $_;
	    chomp($line);
	    $line =~ s/\t+/\t/g;
	    my @f = split(/\t/,$line);
	    my @g = split(/\*/,$f[-1]);
	    if ($g[-1] eq "5"){
		${$t5{$g[0]}}{$g[1]} = $line;
	    }elsif($g[-1] eq "3"){
		${$t3{$g[0]}}{$g[1]} = $line;
	    }
	}
	close(LO);
	$data5{$genomes[$v]} = { %t5 };
	$data3{$genomes[$v]} = { %t3 };
    } #Read all new data;
}

open(OUT,">$opt{$filename[1]}") or die "can't open $opt{$filename[1]}:$!\n";

foreach my $k (keys %input){
    my @f = split(/\t/,$input{$k});
    my ($tmp5,$tmp3);
    my @marker;
    my @tmarker;
    my @dmarker;
    my @gmarker;
    my $outline;    
    for (my $v = 0;$v < 10; $v++){	
	if ($opt{$genomes[$v]}){
	    my $r = "NA";
	    my $t = "NA";
	    my @position;
	    my $chr;
	    my $distance = "NA";
	    my (%m5,%m3);
	    %m5 = %{$data5{$genomes[$v]}};
	    %m3 = %{$data3{$genomes[$v]}};
	    if ($m5{$f[0]} && $m3{$f[0]}){
		foreach my $x (sort {$a<=>$b} keys %{$m5{$f[0]}}){
		    my @tmp = split(/\t/,${$m5{$f[0]}}{$x});
#		    print "${$m5{$f[0]}}{$x}\n";
#		    exit;
		    my $count = 0;
		    if ($tmp[0] =~ /chr/ && $tmp[0] !~ /un/){
			$chr = $tmp[0];
			foreach my $x (sort {$a<=>$b} keys %{$m3{$f[0]}}){
			    my @tmp2 = split(/\t/,${$m3{$f[0]}}{$x});
			    $count ++;
			    if ($tmp2[0] eq $chr && $count < 100){
				push (@position,$tmp2[1]);
				push (@position,$tmp2[2]);
				push (@position,$tmp[1]);
				push (@position,$tmp[2]);
				my @tmparray5 = split(/\*/,$tmp[3]);
				my @tmparray3 = split(/\*/,$tmp2[3]);
				$distance = $tmparray5[1] * 100 + $tmparray3[1] * 100;
				$t = "$tmp[0]:$tmp[1]-$tmp[2].$tmp[3]&$tmp2[0]:$tmp2[1]-$tmp2[2].$tmp2[3]";
				last;
			    }
			}
		    }
		    if (@position){last;}
		}
		if (@position){
		    @position = sort {$a<=>$b} @position;
		    $r = "$chr:$position[0]-$position[-1]";
		    $distance = $distance - ($position[-1] - $position[0]);
		    my $tt = $position[-1] - $position[0];
		}
	    }
	    push (@marker,$r);
	    push (@tmarker,$t);
	    push (@dmarker,$distance);
	    push (@gmarker,$genomes[$v]);
	}
    }
    $outline = join ("\t",@f,join("|",@marker),join("|",@tmarker),join("|",@dmarker),join("|",@gmarker));
    print OUT "$outline\n";    
}
close(OUT);
