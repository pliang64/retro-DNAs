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
getopts("F:B:A:h", \%opt);
my $usage= qq(
This script is used to get pre-integration sequence for Class I entries
Usage :  $0  [Option] bam_infiles
Option:  -F         final list file name
         -A         Assembly
         -h         Help Information

Author : Daniel W. Tang
Date   : Mar. 1st, 2011
\n);


my %data;
my @input;

my @genomes = ("hg38","panTro5","gorGor4","ponAbe2","nomLeu3","chlSab2","macFas5","rheMac8","papAnu2","calJac3");
my %app; 

my %tmp;

my %all;

die "$usage" if (!$opt{F}||!$opt{A}||$opt{h});


open(IN,"$opt{F}") or die "can't open $opt{F}:$!\n";
my $i = 0;
while (<IN>){
    $i += 1;
    my $line = $_;
    chomp($line);
    my @f = split(/\t/,$line);
    my ($id) = $f[0] =~ /(^\d+)/;
    $data{$id} = $line;
}
close(IN);#read final list files and to store the data for Class I entries;

foreach my $x (@genomes){
    opendir my $dir,"/work/wt09rh/seq/$x" or die "Cannot open directory: $!"; 
    my @files = readdir $dir;
    close $dir;
    foreach my $y (@files){
	if ($y =~ ".fa" && !($y =~ /random/ || $y =~ /Un/)){
	    open(REFSEQ, "/work/wt09rh/seq/$x/$y") or die("Can't open reference genome at /work/wt09rh/seq/$x/$y");
	    my @tmpchr = split(/\./,$y);
	    my $j = 0;
	    $tmp{$tmpchr[0]} = "";
	    my $line2;
	    while (<REFSEQ>){
		$j += 1;
		if ($j !=1) {$line2 = $_; chomp($line2); $tmp{$tmpchr[0]} .= $line2;}
	    }
	    close(REFSEQ);
	}
    }
    $all{$x} = {%tmp};
}



foreach my $x (sort {$a<=>$b} keys %data){
    my $testgenome;
    my $selected; 
    my $orthregion;

    my @f = split(/\t/,$data{$x});
    if (!($f[10] eq "DNA" || $f[10] eq "SINE" || $f[10] eq "LINE" || $f[10] eq "Retroposon" || $f[10] eq "LTR" )){next;}#only proceed with entries that have TSD mechanism 
#    if ($f[24] eq "NA|NA|NA|NA|NA|NA|NA|NA|NA"){next;}#next if the entry can't find any thing;
    my @g = split(/\|/,$f[24]);
    my @h = split(/\|/,$f[23]);
    my @k = split(/\|/,$f[22]);
    my @l = split(/\|/,$f[25]);
    my @m = split(/\|/,$f[26]);

#    print "array g = $f[24]\n array h = $f[23]\n array k = $f[22]\n array l = $f[25]\n array m = $f[26]\n"; exit;
    my $j = -1;

    my ($OCHR,$OGS,$OGE) = $f[4] =~ /(chr.+?)\:(\d+)\-(\d+)/;
    my ($PCHR,$PGS,$PGE);
    my $PSEQ;

#    for (my $y = 0; $y <= 9; $y++){
##	print "$l[$y]\n";
#	$app{$l[$y]} = 1;
#    }
    $app{$l[0]} = 1;
#    print STDERR "processing $f[0]\n";
#    for (my $y = 0; $y <= 8; $y++){
#	print STDERR "$g[$y]\n$k[$y]\n$m[$y]\n";
#	if ($g[$y] ne "NA" && abs($g[$y]) <= 10000 && !($k[$y] =~ /random/ || $k[$y] =~ /Un/) && $m[$y] eq "+" ){
    if ($g[0] ne "NA" && abs($g[0]) <= 10000 && !($k[0] =~ /random/ || $k[0] =~ /Un/) && $m[0] eq "+" ){
#	    $j = $y;
	$j = 0;
##	    print STDERR "$j\n";
#	    ($PCHR,$PGS,$PGE) = $k[$y] =~ /(chr.+?)\:(\d+)\-(\d+)/;
	($PCHR,$PGS,$PGE) = $k[0] =~ /(chr.+?)\:(\d+)\-(\d+)/;
##	    %tmp = %{$all{$l[$y]}};
##   	    $PSEQ = substr($tmp{$PCHR},$PGS-1,$PGE-$PGS+1);
#	    $PSEQ = substr(${$all{$l[$y]}}{$PCHR},$PGS-1,$PGE-$PGS+1);
#	$PSEQ = substr(${$all{$l[$y]}}{$PCHR},$PGS-1,$PGE-$PGS+1);
	$PSEQ = substr(${$all{$l[0]}}{$PCHR},$PGS-1,$PGE-$PGS+1);
#	$orthregion = "$l[$y]|$PCHR:$PGS-$PGE";
	$orthregion = "$l[0]|$PCHR:$PGS-$PGE";
#	my ($tmp1,$tmp2) = split(/&/,$h[$y]);
	my ($tmp1,$tmp2) = split(/&/,$h[0]);
##	    print "$tmp1\n$tmp2\n";
	my ($t1) = $tmp1 =~ /\*(\d+)\*/;
	my ($t2) = $tmp2 =~ /\*(\d+)\*/;
##	    print "$t1\n$t2\n";
	$OGS = $OGS - ($t1 * 100);
	$OGE = $OGE + ($t2 * 100);
	$selected = $t1 * 100 . "|" . $t2 * 100;
#	last;
    }
#    }
    
    if ($j < 0){next}#next if can't find a suitable fit for out-group genomes;
#    if (!$testgenome){
#	foreach my $y (@genomes){
#	    if (!$app{$y}){
#		$testgenome = $y;
##	    print "$testgenome\n"; exit;
#	    }
#	}
#    }
#    %tmp = %{$all{$testgenome}};   
    
#    my $OSEQ = substr($tmp{$OCHR},$OGS-1,$OGE-$OGS+1);

    $testgenome = $opt{A};

    my $OSEQ = substr(${$all{$testgenome}}{$OCHR},$OGS-1,$OGE-$OGS+1);
    splice @f, 20, 2;
    splice @f, 5, 3;
    splice @f, 2, 2;
    print join("\t",@f),"\t$selected\t$orthregion\n";
    print ">$f[0] insertion allele\n";
    print "$OSEQ\n";
    print ">$f[0] pre-integration site allele\n";
    print "$PSEQ\n";
    print "//\n";
}  
