#!/usr/local/bin/perl -w
#Tang_find_typical_RE: identify if a typical RE insertion is human-speicific by compare it to the chimp genome

use strict;
use warnings;
use Getopt::Std;

my %POS;
sub sub1{
    my $line1 =shift;
    my $value = 0;
    my @f = split(/|/,$line1);
    if ($f[13] eq "T"){$value++};
    if ($f[14] eq "T"){$value++};
    if ($f[15] eq "A"){$value++};
    if ($f[16] eq "A"){$value++};
    if ($f[17] eq "A"){$value++};
    if ($f[18] eq "A"){$value++};
#    if ($f[19] eq "A"){$value++};

    return($value);
};


my %opt;
getopts("F:B:h", \%opt);
my $usage= qq(
This script is used to choose sites 
Usage :  $0  [Option] bam_infiles
Option:  -F         final list file name
         -h         Help Information

Author : Daniel W. Tang
Date   : Mar. 1st, 2011
\n);


die "$usage" if (!$opt{F}||$opt{h});


open(IN,"$opt{F}") or die "can't open $opt{F}:$!\n";
my $i = 0;
while (<IN>){
    my $line = $_;
    chomp($line);
    my @f = split(/\t/,$line);
    if (sub1(uc($f[1])) > sub1(uc($f[2]))){
	print "$f[0]\t$f[1]\n";
    }else{
	print "$f[0]\t$f[2]\n";
    }
}
close(IN);#read final list files and to store the data for Class I entries;
