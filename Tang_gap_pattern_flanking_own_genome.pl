#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

my %opt;
GetOptions("-i=s%" => \%opt);

my @genomes = ("hg38","panTro5","gorGor4","ponAbe2","nomLeu3","chlSab2","macFas5","rheMac8","papAnu2","calJac3");

my $usage= qq(
This script is used to compare REs position in human genome with gap information 
Usage :  $0  [Option] bam_infiles
Option:  --i hg=(file location..)
         (input,pt,pa,gg,rm,nl)
         -h         Help Information

Author : Daniel W. Tang
Date   : Nov. 22nd, 2011
\n);

my %POS;
my %input;
my %data;
my $count = 0;

if ($opt{input}){
    open(IN,$opt{input}) or die "can't open $opt{input}:$!\n";
    while(<IN>){
	my $line = $_;
	chomp($line);
	$line =~ s/\t+/\t/g;
	my @f = split(/\t/,$line);
	my ($id) = $f[0];
	$input{$id} = $line;
	$count ++;
#    printf STDERR "Reading %9d entries\n", $count;
#    print STDERR "\e[A";
    }
    close(IN); #Read input data;
}else{
    die "$usage";
}

for (my $v = 0;$v < 10; $v++){
    if ($opt{$genomes[$v]}){
	my %th;
        open(LO,"$opt{$genomes[$v]}") or die "can't open $opt{$genomes[$v]}:$!\n";
	while(<LO>){
	    $count ++;
#	printf STDERR "Reading %9d entries\n", $count;
#	print STDERR "\e[A";
	    my $line = $_;
	    chomp($line);
	    $line =~ s/\t+/\t/g;
	    my @f = split(/\t/,$line);
	    unshift @{$th{$f[0]}},{S=>$f[1],E=>$f[2]};
	}
	close(LO);
	foreach my $chr (keys %th){
	    @{$th{$chr}}=sort { $a->{S} <=> $b->{S} || $a->{E} <=> $b->{E} } @{$th{$chr}};	
	}
	$data{$genomes[$v]} = { %th };
    } #Read all gap data;
}

foreach my $k (keys %input){
    my @f = split(/\t/,$input{$k});
#    my @orthologousregion = split(/\|/,$f[-3]);
#    my @order = split(/\|/,$f[-1]);
    my @pattern;
    my ($chr,$s,$e);
    $chr=$f[5];
    for (my $v = 0;$v < 10; $v++){
	if ($opt{$genomes[$v]}){
	    %POS = %{$data{$genomes[$v]}};
	}
    }

#    print STDERR $f[6]-100,"\t",$f[6]-1,"\n";
#    print STDERR $f[7]+1,"\t",$f[7]+100,"\n";
    push(@pattern,binSearch({S=>$f[6]-100,E=>$f[6]-1},$chr));
    push(@pattern,binSearch({S=>$f[7]+1,E=>$f[7]+100},$chr));
    

    my $new =join("\t",(@f,join("|",@pattern)));
    print "$new\n";
}

#check overlap between current each entry and gap regions
sub binSearch {
    my ($p,$chr)=@_;
    if (!$POS{$chr} || !@{$POS{$chr}}){return "+";} #no gap data
    my @a = @{$POS{$chr}};
    my ($fr,$lr)=(0,$#a);
#    my $i;
#    while ($fr <=$lr){
#        $i= int(($fr+$lr)/2);
#        if ($a[$i]->{S} > $p->{E}){$lr =$i-1}
#        elsif($a[$i]->{E} < $p->{S}){$fr =$i+1}
#        elsif(($a[$i]->{S} < $p->{S}||$a[$i]->{S} == $p->{S}) && ($a[$i]->{E} > $p->{E}||$a[$i]->{E} == $p->{E})){
#	    return "-";
#	}
#	else{
#	    return "-";
#	}
#    }
    foreach my $i (@a){
	if (!($i->{E}<$p->{S}||$p->{E}<$i->{S})){
	    return "-";
	}
    }
#    print STDERR "return is from here, no overlap detected\n"; 
    return "+";
    
}
