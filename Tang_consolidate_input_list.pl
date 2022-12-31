#!/usr/bin/perl -w
# Tang_sort_input_new.pl is used to consolidate repeat masker input file to get new version of input file (Containing Alu, L1, SVA, LTR*)
#
#* Adjacent LTRs(within 50bps range) are combined together and treated as one unique entry in later study

use strict;
use DBI;
use Getopt::Std;

my %opt;
getopts("i:p:h:",\%opt);
if ($opt{h} || !$opt{i}){
    print "Usage: $0 [-i inputfile] [-p Percentage cutoff type default is 80] \n"; exit;
}
my $REcutoff;
if (!$opt{p}){
    $REcutoff = 0.8;
}else{
    ($REcutoff) = $opt{p} =~ /(\d+)/;
    if (!$REcutoff || $REcutoff < 0 || $REcutoff > 100){
	$REcutoff = 0.8;
    }else{
	$REcutoff = $REcutoff / 100;
    }
}

my $i = 0;
my $family;
my @g;
my @data;
my @type;
my $totaln = 0;
my $sortn = 0;
my %out;
my %mark;
open(IN,"$opt{i}") or die "Error open $opt{i}: $!\n";
while (<IN>){
    my $line = $_;
    chomp($line);
    my @f = split(/\t/,$line);
    $i += 1;
    push(@g,$line);
#    printf STDERR "reading %9d entries\n", $i;
#    print STDERR "\e[A";
}
close(IN);
foreach my $k (@g){
    my @f = split(/\t/,$k);
    if (!($f[5] =~ /_/i)){ #REs in chr*_random are not selected 
	if ($f[10] =~ /Alu/i || $f[12]=~ /L1/i || $k =~ /LTR/i || $k =~ /SVA/i){
	    push(@data,$k);
	    $sortn ++;
	}
    }
}#attract all REs;
my $j;
my $count = 0;
my $z = 0;
for ($j = 0; $j < scalar(@data); $j++){
    $z ++;
#    printf STDERR "processing %9d entries\n", $z;
#    print STDERR "\e[A";
    if (!$mark{$j}){
	my @f = split(/\t/,$data[$j]);
	my @test;
	my $REratio;
	my $ID;
	if ($f[9] eq "+"){
	    $REratio = ($f[14] - $f[13] + 1)/($f[14] - $f[15]);
	}else{
	    $REratio = ($f[14] - $f[15] + 1)/($f[14] - $f[13]);
	}
	if ($REratio > $REcutoff){
	    $count += 1;
	    $ID = sprintf("%08d", $count);
	    push(@test,$j);
	    foreach my $m (@test){$mark{$m} = 1};
	    $out{$ID} = output($ID,@test);
	}else{
	    $count += 1;
	    $ID = sprintf("%08d", $count);
	    @test = check($j);
	    foreach my $m (@test){$mark{$m} = 1};
	    $out{$count} = output($ID,@test);    
	}# check if entry is complete. if so, treat as one entry;
    }
}

foreach my $k (sort {$a <=> $b} keys %out){
    print "$out{$k}\n";
}
print STDERR "Processed $opt{i}: Total $sortn entries were selected out of $i and combined to $count unique entries\n";


sub check{
    my @tmpout;
    my $x = shift;
    push(@tmpout,$x);
    $mark{$x} = 1;
    my $i;
    for ($i = 1; $i < 21; $i++){#search 20 REs downstream;
	my $y = $x + $i;
	my @f = split(/\t/,$data[$x]);
	if ($y < scalar(@data)){
	    my @g = split(/\t/,$data[$y]);
	    if ($f[5] eq $g[5] && $f[9] eq $g[9] && $f[10] eq $g[10] && ($g[6] - $f[7]) < 50000 && !$mark{$y}){
		if ($f[9] eq "+"){
		    if ($f[14] <= $g[13] || ($f[14] >= $g[13] && $f[13] < $g[13] && (($f[14] - $g[13])/($g[14] - $g[15])) < 0.1 && $g[14] > $f[14])){
			foreach my $k (check($y)){
			    push(@tmpout, $k);
			    $mark{$k} = 1;
			}
			last;
		    }
		}
		else{
#		    if ($g[15] - $g[13] == 0){die "$data[$y]\n";}
		    if ($g[14] <= $f[15] || ($g[14] >= $f[15] && $g[15] < $f[15] && (($g[14] - $f[15])/($g[14] - $g[13])) < 0.1 && $f[14] > $g[14])){
			foreach my $k2 (check($y)){
			    push(@tmpout, $k2);
			    $mark{$k2} = 1;
			}
			last;
		    }
		}
	    }
	}
    }
    return(@tmpout);
}

sub output{
    my ($count,@p) = @_;
    my $line;
    if (!((scalar(@p) > 1))){
	my @t = split(/\t/,$data[$p[0]]);
	if ($t[9] eq "+"){
	    $line = join("\t",($count,@t[5..7,9..12],1,@t[6..7,13..14,15]));
	}else{ 
	    $line = join("\t",($count,@t[5..7,9..12],1,@t[6..7,14..15,13]));
	}
    }else{
	@p = sort {$a <=> $b} @p;
	my $blocknumber = scalar(@p);
	my $genostart;
	my $genoend;
	my $restart;
	my $reend;
	foreach my $k (@p){
	    my @f = split(/\t/,$data[$k]);
	    $genostart .= $f[6].",";
	    $genoend .= $f[7].",";
	    if ($f[9] eq "+"){
		$restart .= $f[13].",";
		$reend .= $f[14].","; 
	    }else{
		$restart .= $f[14].",";
		$reend .= $f[15].",";
	    }
	}
	my $renotmatched;
	my @g = split(/\t/,$data[$p[0]]);
	my @h = split(/\t/,$data[$p[scalar(@p)-1]]);
	if ($h[9] eq "+"){
	    $renotmatched = $h[15];
	}
	else{
	    $renotmatched = $g[13];
	}
	my $tmpn = scalar(@p);
	$line = join("\t",($count,@g[5..6],$h[7],@g[9..12],scalar(@p),$genostart,$genoend,$restart,$reend,$renotmatched));
#	$line = "$count\t$g[5]\t$g[6]\t$h[7]\t$g[9]\t$g[10]\t$g[11]\t$g[12]\t$tmpn\t$genostart\t$genoend\t$restart\t$reend\t$renotmatched\t";
    }
    return($line);
}
