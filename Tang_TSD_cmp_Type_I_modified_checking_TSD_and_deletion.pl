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
getopts("S:A:O:h", \%opt);
my $usage= qq(
This script is used to get pre-integration sequence for Class I entries
Usage :  $0  [Option] bam_infiles
Option:  -S         TSD sequence file name
         -A         Genome assebmly version
         -O         output mode, default is 0, without sequences in fasta format
         -h         Help Information

Author : Daniel W. Tang
Date   : Mar. 1st, 2011
\n);
my %chr;
my %data;
my @input;
my %pt3;
my %pa2;
my %gg3;
die "$usage" if (!$opt{S}||!$opt{A}||$opt{h});
if (!$opt{O}){$opt{O} = 0};


open(IN,"$opt{S}") or die "can't open $opt{S}:$!\n";
my @filepre = split(/\./,$opt{S});
my $i = 0;
my $assembly = $opt{A};
my %refseq;
my $errorcount = 0;
my $previousfalsenegative = 0;
while (<IN>){
    my $tfcondition = 0;
    my $RIMD = "NA";
    my $RIMDlen = 0;
    my $seq = "";
    my $line = $_;
    my $conditionRIMD = 0;
    while (!($line =~ /\/\//)){
        if (!($line eq "\n")){$seq .= $line;}
        $line = <IN>;
    }
    if (!($seq =~ /^>/)){$seq = ">".$seq};
    my @f = split(/>/,$seq);
    my @go = split(/\t/,$f[1]);
    my ($id) = $go[0] =~ /(\d+)/;
    $id = "$assembly"."$id";
#    my ($OGS,$OGE) = split(/\|/,$go[-1]);
    my ($OGS,$OGE) = split(/\|/,$go[-2]);
    splice @go, 15, 7;

    $go[15] = "NA";
    $go[16] = 0;
    $go[17] = "NA";
    $go[18] = "NA";
    $go[19] = 0;
    $go[20] = "NA";
    $go[21] = 0;
    $go[22] = "NA";
    $go[23] = 0;#initializing;
    $i++;
#    printf STDERR "Processing %9d entries\n", $i;
#    print STDERR "\e[A";
    my $tmpfile = $id.'insertion.fas';
    my $tmpreferrence = $id.'pre_integration.fas';
    my $bl2seqoutfile = $id.'bl2seqout.psl';
    my @tmpinsertionseq = split(/\n/,$f[2]);
    my $tmpx;
    my ($tmp,$tmp2);
    for ($tmpx = 1;$tmpx < scalar(@tmpinsertionseq);$tmpx++){
	$tmp2 .= $tmpinsertionseq[$tmpx];
    }
    my $ilen = length($tmp2);
#    print "$ilen\n";
    open(TMP, ">$tmpfile") or die "Error open $tmpfile: $!\n";
    print TMP (">$id|insertion\n$tmp2\n");
    close(TMP);
    ($tmp,$tmp2) = split(/\n/,$f[3]);
    open(TMPSEQ2, ">$tmpreferrence");
    print TMPSEQ2 (">$id|pre\n$tmp2\n");
    close(TMPSEQ2);
#    exit;
    system("/work/lianglab/bin/blast-2.2.26/bin/bl2seq -i $tmpreferrence -j $tmpfile -p blastn -e 10e-1 -D 1 -r 1 -q -2 -F F -m T -W 14 -o $bl2seqoutfile");
#    print STDERR "$tmpreferrence\t$tmpfile\n";
#    system ("bl2seq -i $tmpreferrence -j $tmpfile -p blastn -g F -o $bl2seqoutfile -W 9 -F F ") && die ("Problem running bl2seq: $!");
#    exit;
    my $tsdstart = 0;
    my $tsdend = length($tmp2)-1;
    my ($p1,$p2);
    my $tsd = "NNNNN";
    my $p5tsd = "NA";
    my $p3tsd = "NA";
    my $pos = $go[2];
#    my $pos = $go[4];
    $pos =~ s/,//;
    $pos =~ tr/:-/\t/;
    my ($chr,$s,$e) = split(/\t/, $pos);
    if (!$refseq{$assembly}{$chr}){
        open(REFSEQ, "/work/wt09rh/seq/$assembly/$chr.fa") or die("Can't open reference genome file");
        my $j = 0;
        $refseq{$assembly}{$chr} = "";
        my $line2;
        while (<REFSEQ>){
            $j += 1;
            if ($j !=1) {$line2 = $_; chomp($line2); $refseq{$assembly}{$chr} .= $line2;}
        }
        close(REFSEQ);
    }

    my $p5flank = lc(substr($refseq{$assembly}{$chr},$s-101,100));
    my $p3flank = lc(substr($refseq{$assembly}{$chr},$e,100));
    my $reseq = uc(substr($refseq{$assembly}{$chr},$s-1,$e-$s+1));

    my $PSEQ;

    my $j = 0;
    
    my $tmprelen = $e - $s + 1;    

    my ($i_motif,$i_motif_2);
#    open(BLASTOUT,"$bl2seqoutfile") or die "Can not open $bl2seqoutfile";
    open(BLASTOUT,"$bl2seqoutfile") or next;
    my @blocks;
    my (@lf,@rf);
#    print STDERR "ilen: $ilen\t tsdend: $tsdend\n";
#    print STDERR "$OGS\t$OGE\n";
#    print STDERR "$s\t$e\n";
    my $direction = 0;
    while (<BLASTOUT>){
#	print STDERR "There is BLASTOUT\n";
	$line = $_;
	chomp($line);
	if (!($line =~/^\#/)){
	    my @g = split(/\t/,$line); 
#	    print STDERR "$line\n";
#	    next if ($g[2] < 90 );
	    next if ($g[2] < 85 );
	    if (($g[8] >= $OGS + $e - $s + 1 || $g[9] >= $OGS + $e - $s + 1 )&& ($g[8] >= $OGS + $e - $s + 1 - 20 && $g[9] >= $OGS + $e - $s + 1 - 20)){ #backup Aug.17th
		if ($direction == 0){
		    $direction = ($g[9] - $g[8])*($g[7]-$g[6]);
		}
		if(($g[9] - $g[8])*($g[7]-$g[6])*$direction < 0){
		    next;
		}
		push(@rf,\@g);
#		print STDERR "A:rf\n";
		next;
	    }
	    if (($g[8] <= $OGS || $g[9] <= $OGS) && ($g[8] <= $OGS+20 && $g[9] <= $OGS+20)){
		if ($direction == 0){
		    $direction = ($g[9] - $g[8])*($g[7]-$g[6]);
		}
		if(($g[9] - $g[8])*($g[7]-$g[6])*$direction < 0){
		    next;
		}
		push(@lf,\@g);
#		print STDERR "B:lf\n";
		next;
	    }
#	    if ($g[8] >= $OGS + $e - $s + 1  || $g[9] >= $OGS + $e - $s + 1 ){
#		push (@rf,\@g);
#		print STDERR "C:rf\n";
#		next;
#	    }
#	    if ($g[8] <= $OGS || $g[9] <= $OGS){
#		push (@lf,\@g);
#		print STDERR "D:lf\n";
#		next;
#	    }
#	    print STDERR "discarded\n";
	}
    }
    @rf = sort {$b->[11] <=> $a->[11]} @rf;
    @lf = sort {$b->[11] <=> $a->[11]} @lf;

    my (@ttmpg1,@ttmpg2);

    if (@rf && @lf){
	my $tmpline = join("\t",@{$rf[0]});
#	print STDERR $OGS + $e - $s + 1,"\n";
#	print STDERR "rf:\n$tmpline\n";
	push(@blocks,$tmpline);
	$tmpline = join("\t",@{$lf[0]});
#	print STDERR "$OGS\n";
#	print STDERR "lf:\n$tmpline\n";
	push(@blocks,$tmpline);
	@ttmpg1 = @{$rf[0]};
	@ttmpg2 = @{$lf[0]};
    }


    my $bn = scalar(@blocks);
    my $sl = length ($seq);
#    print STDERR "$bn\n";
#    exit;

#    print STDERR ("\n\n\n$bn\n\n");
#    exit;
    close(BLASTOUT);
    my $outline;
    if ($bn < 2){
	$#lf = -1;
	$#rf = -1;
	$#blocks = -1;#clear all elements in the @lf and @rf;
	unlink ($bl2seqoutfile);
	system("/work/lianglab/bin/blast-2.2.26/bin/bl2seq -i $tmpreferrence -j $tmpfile -p blastn -q -3 -e 10e-1 -F F -D 1 -o $bl2seqoutfile");
#	exit;
	open(BLASTOUT,"$bl2seqoutfile") or die "Can not open $bl2seqoutfile";
	my $direction2 = 0;
	while (<BLASTOUT>){
	    $line = $_;
	    chomp($line);
	    if (!($line =~/^\#/)){
		my @g3 = split(/\t/,$line);
#		next if ($g3[2] < 90);
		next if ($g3[2] < 85);
#		print STDERR "$line\n";
		if (($g3[8] >= $OGS + $e - $s + 1 || $g3[9] >= $OGS + $e - $s + 1 )&& ($g3[8] >= $OGS + $e - $s + 1 - 20 && $g3[9] >= $OGS + $e - $s + 1 - 20)){
		    if ($direction2 == 0){
			$direction2 = ($g3[9] - $g3[8])*($g3[7]-$g3[6]);
		    }
		    if(($g3[9] - $g3[8])*($g3[7]-$g3[6])*$direction2 < 0){
			next;
		    }
		    push(@rf,\@g3);
#		    print STDERR "A:rf\n";
		    next;
		}
		if (($g3[8] <= $OGS || $g3[9] <= $OGS) && ($g3[8] <= $OGS+20 && $g3[9] <= $OGS+20)){
		    if ($direction2 == 0){
			$direction2 = ($g3[9] - $g3[8])*($g3[7]-$g3[6]);
		    }
		    if(($g3[9] - $g3[8])*($g3[7]-$g3[6])*$direction2 < 0){
			next;
		    }
		    push(@lf,\@g3);
#		    print STDERR "B:lf\n";
		    next;
		}
	    }
	}
	@rf = sort {$b->[11] <=> $a->[11]} @rf;
	@lf = sort {$b->[11] <=> $a->[11]} @lf;
	if (@rf && @lf){
	    my $tmpline = join("\t",@{$rf[0]});
#	    print STDERR "$tmpline\n";
	    push(@blocks,$tmpline);
	    $tmpline = join("\t",@{$lf[0]});
#	    print STDERR "$tmpline\n";
	    push(@blocks,$tmpline);
	    @ttmpg1 = @{$rf[0]};
	    @ttmpg2 = @{$lf[0]};
	}
	$bn = scalar(@blocks);
    }#if the initial megablast is too strict and couldn't properly detect the 2 blocks, use a less blastn
    if ($bn>1&&(($ttmpg1[6] - $ttmpg2[6])*($ttmpg1[7] - $ttmpg2[7]))> 0){
#	print STDERR "Got to this step\n";
#	exit;
	my %cordinates;
        my @set;
	my @tmpcor;
        foreach my $each (@blocks){
            my @g2 = split(/\t/,$each);
            $cordinates{$g2[6]}= $g2[8];
            $cordinates{$g2[7]}= $g2[9];
	    push(@tmpcor,$g2[8]);
	    push(@tmpcor,$g2[9]);
            my @g3;
            @g3 = ($g2[6]..$g2[7]);
#	    @g3 = ($g2[6],$g2[7]);
            push(@set,[@g3]);
        }
	@tmpcor = sort {$a <=> $b} @tmpcor;
	my %a = map{$_ => 1} @{$set[0]};
	my %b = map{$_ => 1} @{$set[1]};
	my @inter;
	@inter = sort {$a <=> $b} grep {$a{$_}} @{$set[1]};
#	exit;
	if (@inter) {             
	    $conditionRIMD = 1;
	}
        if (@inter && ($inter[$#inter] - $inter[0]) <= 30) {  
#	    print STDERR "Inter happens\n";
#	    exit;
	    $conditionRIMD = 1;
            $p1 = $inter[0];
            $p2 = pop(@inter);
#	    print STDERR "$p1\t$p2\n";
#	    print STDERR "$cordinates{$p1}\t$cordinates{$p2}\n";
#	    exit;
            open(GETTSD,"$tmpreferrence");
	    open(INSERTION,"$tmpfile");
            my $gettsd = <GETTSD>;
            $gettsd = <GETTSD>;
	    my $getinsertion = <INSERTION>;
	    $getinsertion = <INSERTION>;
	    chomp($getinsertion);
	    chomp($gettsd);

#	    print STDERR ("\n\n\n$p1\t$cordinates{$p1}\t$p2\t$cordinates{$p2}\t\n");
	    $tsd = substr($gettsd,$p1-1,$p2-$p1+1);
	    if ($cordinates{$p1} < $cordinates{$p2}){
		$p5tsd = substr($getinsertion,$cordinates{$p1}-length($tsd),length($tsd));
		$p3tsd = substr($getinsertion,$cordinates{$p2}-1,length($tsd));
		$tsd = sub1($tsd);
		$p5flank = lc(substr($refseq{$assembly}{$chr},$s-$OGS+$cordinates{$p1}-length($tsd)-101,100));
		$p3flank = lc(substr($refseq{$assembly}{$chr},$s-$OGS+$cordinates{$p2}+length($tsd)-2,100));
		$reseq = uc(substr($getinsertion,$cordinates{$p1},$cordinates{$p2}-$cordinates{$p1}-1));
	    }else{
                $p5tsd = substr($getinsertion,$cordinates{$p2}-length($tsd),length($tsd));
                $p3tsd = substr($getinsertion,$cordinates{$p1}-1,length($tsd));
		$p5flank = lc(substr($refseq{$assembly}{$chr},$s-$OGS+$cordinates{$p2}-length($tsd)-101,100));
		$p3flank = lc(substr($refseq{$assembly}{$chr},$s-$OGS+$cordinates{$p1}+length($tsd)-2,100));
		$reseq = uc(substr($getinsertion,$cordinates{$p2},$cordinates{$p1}-$cordinates{$p2}-1));
	    }


	    if ($go[3] eq "+"){
		$i_motif = substr($gettsd,$p1-16,30);
	    }else{
		$i_motif = substr($gettsd,$p2-16,30);
	    }

	    close(GETTSD);
#	    if ($cordinates{$p1}>$cordinates{$p2}){
#		my $tmp = $p1;
#		$p1 = $p2;
#		$p2 = $tmp;
#	    };
#	}
#	else{
#	    my @tmp1;
#	    foreach my $k (keys %cordinates){
#		push(@tmp1,$cordinates{$k});
#	    }
#	    @tmp1 = sort {$a <=> $b} @tmp1;
#	}
#    }
	    my ($p5tr,$p3tr,$p5trlen,$p3trlen) = ("NA","NA",0,0);

	    if (${$lf[0]}[9] > ${$lf[0]}[8]){
		if (${$lf[0]}[9] < $OGS){
		    my $tmp = substr($getinsertion,${$lf[0]}[9],$OGS-${$lf[0]}[9]);
		    if ($go[3] eq "+"){
			$p5tr = $tmp;
			$p5trlen = length($p5tr);
		    }else{
			$p3tr = sub1($tmp);
			$p3trlen = length($p3tr);
		    }
		}
	    }else{
		if (${$lf[0]}[8] < $OGS){
                    my $tmp = substr($getinsertion,${$lf[0]}[8],$OGS-${$lf[0]}[8]);
		    if ($go[3] eq "+"){
			$p5tr = $tmp;
			$p5trlen = length($p5tr);
		    }else{
			$p3tr = sub1($tmp);
			$p3trlen = length($p3tr);
		    }
		}
	    }
            if (${$rf[0]}[9] > ${$rf[0]}[8]){
		if (${$rf[0]}[8] > length($getinsertion) - $OGE + 1){
		    my $tmp = substr($getinsertion,length($getinsertion)-$OGE,${$rf[0]}[8]+$OGE - length($getinsertion) - 1);
		    if ($go[3] eq "+"){
			$p3tr = $tmp;
                        $p3trlen = length($p3tr);
		    }else{
			$p5tr = sub1($tmp);
                        $p5trlen = length($p5tr);
		    }
		}
	    }else{
                if (${$rf[0]}[9] > length($getinsertion) - $OGE + 1){
		    my $tmp = substr($getinsertion,length($getinsertion)-$OGE,${$rf[0]}[9]+$OGE - length($getinsertion) - 1);
                    if ($go[3] eq "+"){
                        $p3tr = $tmp;
                        $p3trlen = length($p3tr);
                    }else{
                        $p5tr = sub1($tmp);
                        $p5trlen = length($p5tr);
                    }
                }
	    }
	    $go[20] = $p5tr;
	    $go[21] = $p5trlen;
	    $go[22] = $p3tr;
	    $go[23] = $p3trlen;

	    if (length($reseq)>($e-$s+1)){$reseq = uc(substr($refseq{$assembly}{$chr},$s-1,$e-$s+1));}

	    $tsd = uc($tsd);
	    $i_motif = uc($i_motif);
	    $p5tsd = uc($p5tsd);
	    $p3tsd = uc($p3tsd);
#    print("$tsd\n\n");
       }#else{ #change Sept. 16th so every blat can search local break points for a more accurate call of TSD
#	print STDERR "Inter did not happen\n";
	my ($REchr,$REs,$REe) = $go[2] =~ /(chr.+?)\:(\d+)\-(\d+)/;
#	    print STDERR "$REchr\t$REs\t$REe\n";
#	    print STDERR "$go[2]\n";
	$go[15] = "NA";
	$go[16] = 0;
	$go[17] = "NA";
	$go[18] = "NA";
	$go[19] = "NA";
	$go[20] = "NA";
	$go[21] = 0;
	$go[22] = "NA";
	$go[23] = 0;

	my $actuallength = $tmpcor[2] - $tmpcor[1] + 1;
	my $MEIlength = $REe - $REs + 1;
#	    if ($actuallength - $MEIlength > 29){
	if ($actuallength){
	    $errorcount ++;
	    #Start processing entries to see if it's transduction
	    my ($s1,$e1,$s2,$e2);
	    my $similarity = 90;
	    my $searchspan = 0;
	    if ($go[5] eq "LTR"){$searchspan = 15}else{$searchspan = 50};
	    if ($go[8] < 100){$searchspan = 15};
	    if (@inter){
		$p1 = $inter[0];
		$p2 = pop(@inter);
#		if (!$p2 || !$p1){print STDERR "$f[1]\n"; die}
		if($p2 - $p1 + 1 < 10){
		    $searchspan = 15
		    }
	    }
#		print STDERR "$go[5]\n";
#		exit;

	    $s1 = $s - $OGS + $tmpcor[1] - $searchspan;
	    $e1 = $s - $OGS + $tmpcor[1] + $searchspan - 1;
	    $s2 = $s - $OGS + $tmpcor[2] - 2 - $searchspan + 1;
	    $e2 = $s - $OGS + $tmpcor[2] - 2 + $searchspan ;
#	    print STDERR "$s1\t$e1\n$s2\t$e2\n";
	    if (!$refseq{$assembly}{$REchr}){
		open(REFSEQ, "/work/wt09rh/seq/$assembly/$REchr.fa") or die("Can't open reference genome file");
		my $j = 0;
		$refseq{$assembly}{$REchr} = "";
		my $line2;
		while (<REFSEQ>){
		    $j += 1;
		    if ($j !=1) {$line2 = $_; chomp($line2); $refseq{$assembly}{$REchr} .= $line2;}
		}
		close(REFSEQ);
	    }
	    open (OUT1, ">$id.seq1");
	    open (OUT2, ">$id.seq2");
	    my $seq1 = substr($refseq{$assembly}{$REchr},$s1-1,$e1-$s1+1);
	    my $seq2 = substr($refseq{$assembly}{$REchr},$s2-1,$e2-$s2+1);
	    print OUT1 ">lf\n$seq1\n";
	    print OUT2 ">rf\n$seq2\n";
	    close(OUT1); close(OUT2);
	    my $secondbl2seqoutfile = "$id".".2nd.bl2seq.out";
#		system("/work/lianglab/bin/blast-2.2.26/bin/bl2seq -i $id.seq1 -j $id.seq2 -p blastn -e 10e-1 -D 1 -r 1 -q -2 -F F -W 4 -o $secondbl2seqoutfile");
	    system("/work/lianglab/bin/blast-2.2.26/bin/bl2seq -i $id.seq1 -j $id.seq2 -q -3 -e 10e-1 -W 4 -p blastn -F F -D 1 -o $secondbl2seqoutfile");
	    open(BLSTOUT2,"$secondbl2seqoutfile") or die "Can not open $secondbl2seqoutfile:$!\n";
	    my ($seq1s, $seq1e, $seq2s, $seq2e,$sim);
	    my $flag=0;
	    my @tmphash;
	    while(<BLSTOUT2>){
		my $tmp2line = $_;
		chomp($tmp2line);
#		print STDERR "$tmp2line\n";
		next if ($tmp2line=~/^\#/);
		my @tmpg = split(/\t/,$tmp2line);
		next if ($tmpg[2]<$similarity);
		next if ((($tmpg[6]-$tmpg[7])*($tmpg[8]-$tmpg[9]))<0);#ignore any reverse match;
#		next if ($go[3] eq "+" && $tmpg[7] > $searchspan * 1.1); 
#		next if ($go[3] eq "-" && $tmpg[8] < $searchspan * 0.90);#ignore any match going to the TE; allow some error room as the TSDs some time can overlap 1 - 2 bp with TE due to inaccurate call of TE by repeat masker;# we don't know how these TSD transposon move, take it out for now?

# Rewrite this section so it won't bite into it?
#		print STDERR "$s\n$e\n";
		next if ($tmpg[7] + $s1 > $searchspan * 0.1 + $s ); 
		next if ($e - $searchspan * 0.1 > $e2 + $tmpg[8]);#ignore any match going to the TE; allow some error room as the TSDs some time can overlap 1 - 2 bp with TE due to inaccurate call of TE by repeat masker;
		next if ($tmpg[3] > 30);
		push(@tmphash,\@tmpg);
		$flag=1;
#		print STDERR "Chosen\n";
#		    last;
	    }
	    close(BLSTOUT2);
	    unlink("$id.seq1");
	    unlink("$id.seq2");
	    unlink("$secondbl2seqoutfile");
	    
	    
	    @tmphash = sort {$b->[11] <=> $a->[11]} @tmphash;
#	    print STDERR "$flag\n";
	    ($seq1s, $seq1e, $seq2s, $seq2e,$sim) = (${$tmphash[0]}[6],${$tmphash[0]}[7],${$tmphash[0]}[8],${$tmphash[0]}[9],${$tmphash[0]}[2]);
	    
	    if ($flag == 1){
		my ($restart,$reend) = ($s,$e);
		if ($s1+$seq1e > $s){$restart = $s1 + $seq1e};
		if ($s2+$seq2s-2 < $e){$reend = $s2+$seq2s-2};
		if ($reend < $restart){
		    $restart = $s; 
		    $reend = $e; 
		    $flag = 0;
		    $tsd = "NNNNN";
		    $p5tsd = "NA";
		    $p3tsd = "NA";
		    $p5flank = lc(substr($refseq{$assembly}{$chr},$s-101,100));
		    $p3flank = lc(substr($refseq{$assembly}{$chr},$e,100));
		    $reseq = uc(substr($refseq{$assembly}{$chr},$s-1,$e-$s+1));
		}
	    }
	    
	    if ($flag == 1){
#		print STDERR "$seq1s\t$seq1e\t$seq2s\t$seq2e\n";
#		exit;
		$tfcondition = 1;
		my $p5tran = "NA";
		my $p5tranlen = 0;
		my $p3tran = "NA";
		my $p3tranlen = 0;
		$previousfalsenegative ++;
		$seq1=uc($seq1);
		$seq2=uc($seq2);
		my @lf_arr = split //, $seq1;
		my @rf_arr = split //, $seq2;
		my @ucltsd = @lf_arr[$seq1s-1..$seq1e-1];
		my @ucrtsd = @rf_arr[$seq2s-1..$seq2e-1];
#		    $p5flank = substr($refseq{$assembly}{$REchr},$s - $OGS + $tmpcor[1]- 100 +$seq1s-102,100);
		$p5flank = substr($refseq{$assembly}{$REchr},$s1 + $seq1s - 1 - 100 - 1,100);
#		    $p3flank = substr($refseq{$assembly}{$REchr},$s - $OGS + $tmpcor[2] - 2 -100 + $seq2e,100);
		$p3flank = substr($refseq{$assembly}{$REchr},$s2 + $seq2e - 1,100);
#		    print STDERR "$p5flank\n";
		$reseq = substr($refseq{$assembly}{$REchr},$s-1,$e-$s+1);
		$tsd = join("",@ucltsd);
		$p5tsd = join("",@ucltsd);
		$p3tsd = join("",@ucrtsd);
		my $tsdlength = length($tsd);
		if ($go[3] eq "+"){
		    if ($tsdlength >= 15){
			$i_motif = substr($refseq{$assembly}{$REchr},$seq1s - 15 + $s1 - 1 - 1,30);
		    }elsif($tsdlength <15){
			$i_motif = substr($refseq{$assembly}{$REchr},$seq1s - 15 + $s1 - 1 - 1,15+$tsdlength).substr($refseq{$assembly}{$REchr},$seq2e + 1 + $s2 - 1 - 1,15-$tsdlength);
		    }
		}else{
		    if ($tsdlength >= 15){
			$i_motif = substr($refseq{$assembly}{$REchr},$seq2e - 14 + $s2 - 1 - 1,30);
		    }elsif($tsdlength <15){
			$i_motif = substr($refseq{$assembly}{$REchr},$seq1s-(15-$tsdlength) + $s1 - 1 - 1,15-$tsdlength).substr($refseq{$assembly}{$REchr},$seq2s + $s2 - 1 - 1,15+$tsdlength);
		    }
#			$i_motif = sub1($i_motif);
#			$tsd = sub1($tsd);
#			my $tmptsd = $p5tsd;
#			$p5tsd = sub1($p3tsd);
#			$p3tsd = sub1($tmptsd);
#			my $tmpflank = $p5flank;
#			$p5flank = sub1($p3flank);
#			$p3flank = sub1($tmpflank);
#			$reseq = sub1($reseq);
		}
		my ($restart,$reend) = ($s,$e);
		
		if ($s1+$seq1e > $s){$restart = $s1 + $seq1e};
		if ($s2+$seq2s-2 < $e){$reend = $s2+$seq2s-2};
		if ($reend < $restart){$restart = $s; $reend = $e;}
		$reseq = substr($refseq{$assembly}{$REchr},$restart-1,$reend - $restart + 1);
		
		$i_motif = uc($i_motif);
		my $tmptrans;
		$p5flank = lc($p5flank);
		$p3flank = lc($p3flank);
#		    if ($s - $OGS + $tmpcor[1] - 1 < $s - 1){ #comment out;
#		    print STDERR "$OGS\n";
#		    print STDERR $s - $OGS + $tmpcor[2] - 2,"\t",$e,"\t",$seq1e,"\n";
		if($seq1e + $s1 < $s){
#			print STDERR "1st cut ",$seq1e+$s1,"\t",$s-1,"\n";
		    $tmptrans = substr($refseq{$assembly}{$REchr},$seq1e + $s1 - 1,$s -1 - ($seq1e + $s1) + 1);
#			$p5tran = $tmptrans;
#			$p5tranlen = length($p5tran);
		    if ($go[3] eq "+"){
			$p5tran = $tmptrans;
			$p5tranlen = length($p5tran);
		    }else{
			$p3tran = sub1($tmptrans);
			$p3tranlen = length($p3tran);
		    }
		}
		if ($s2 + $seq2s  - 2 > $e){
#			print STDERR "2nd cut ",$e+1,"\t",$s2+$seq2s-1,"\n";
		    $tmptrans = substr($refseq{$assembly}{$REchr},$e,$s2 + $seq2s - 2 - ($e + 1) + 1);
#			$p3tran = $tmptrans;
#			$p3tranlen = length($p3tran);
#			print STDERR "$go[3]\tcoming from here : $tmptrans\n";
		    
		    if ($go[3] eq "+"){
			$p3tran = $tmptrans;
			$p3tranlen = length($p3tran);
		    }else{
			$p5tran = sub1($tmptrans);
			$p5tranlen = length($p5tran);
		    }
		}
#		    print STDERR "5 flanking transduction : $p5tran\n 3 flanking transduction : $p3tran\n";
		
		$go[20] = $p5tran;
		$go[21] = $p5tranlen;
		$go[22] = $p3tran;
		$go[23] = $p3tranlen;
		
	    }
	    
	    
	    
	    
	    
	    #1st step: calculate break positions in the human genome only; Retrieve breakpoint sequences
	    if ((4 ** length($tsd))/4 < ($searchspan * 2)){
		$tsd = "NNNNN";
		$p5tsd = "NA";
		$p3tsd = "NA";
		$p5flank = lc(substr($refseq{$assembly}{$chr},$s-101,100));
		$p3flank = lc(substr($refseq{$assembly}{$chr},$e,100));
		$reseq = uc(substr($refseq{$assembly}{$chr},$s-1,$e-$s+1));

	    }

	}
	my @tmparray;
	foreach my $each (@blocks){
		my @g2 = split(/\t/,$each);
		push(@tmparray,$g2[6]);
		push(@tmparray,$g2[7]);
	    }
	@tmparray = sort {$a <=> $b} @tmparray;
	if ($conditionRIMD < 1){
	    $RIMD = uc(substr($tmp2,$tmparray[1],$tmparray[2]-$tmparray[1] - 1));
	    $RIMDlen = length($RIMD);
	    if ($RIMD =~ /NN/||$RIMDlen == 0){
		$RIMD = "NA";
		$RIMDlen = 0;
		}
	}
#	} #change Sept. 16th so every blat can search local break points for a more accurate call of TSD
	
	if (!($tsd =~ /NNNNN/)){
	    my $tsdlength = length($tsd);
	    my @tmph = split(/\t/,$tmpinsertionseq[0]);
	    my $motiflength = length($i_motif);
	    if ($go[3] eq "-"){
		$tsd = sub1($tsd);
		$i_motif = sub1($i_motif);
		my $tmptsd = $p5tsd;
		$p5tsd = sub1($p3tsd);
		$p3tsd = sub1($tmptsd);
		my $tmpflank = $p5flank;
		$p5flank = sub1($p3flank);
		$p3flank = sub1($tmpflank);
		$reseq = sub1($reseq);
	    }		
	    $go[15] = $tsd;
	    $go[16] = length($tsd);
	    $go[17] = $p5tsd;
	    $go[18] = $p3tsd;
	    $go[19] = $i_motif;
#	    $i_motif = uc($i_motif);
#	    print("$tmpinsertionseq[0]\t$tsd\t$i_motif\t$tsdlength\t\n");	    
	}else{
            if ($go[3] eq "-"){
                my $tmpflank = $p5flank;
                $p5flank = sub1($p3flank);
                $p3flank = sub1($tmpflank);
                $reseq = sub1($reseq);
            }
	}
#	next; #debuggging control point, comment out

	$outline = join("\t",@go[0..23],$RIMD,$RIMDlen);

    }else{
	if ($go[3] eq "-"){
	    my $tmpflank = $p5flank;
	    $p5flank = sub1($p3flank);
	    $p3flank = sub1($tmpflank);
	    $reseq = sub1($reseq);
	}
	$outline = join("\t",(@go[0..14],"NA","0","NA","NA","NA","NA","0","NA","0","NA","0"));
    }    
    if ($opt{O} == 0){
	print "$outline\n";
    }elsif($opt{O} == 1){
	$outline = join("|",@go[0..23],$RIMD,$RIMDlen);
	print ">$outline\n";
	print "$p5flank\n";
	if($go[17] eq "NA"){print "\n";}else{print"$go[17]\n"};
	if($go[20] eq "NA"){print "\n";}else{print"$go[20]\n"};
	print "$reseq\n";
	if($go[22] eq "NA"){print "\n";}else{print"$go[22]\n"};
	if ($go[18] eq "NA"){print "\n";}else{print "$go[18]\n"};
	print"$p3flank\n";
#	if ($tfcondition > 0){
#	    print STDERR ">$outline\n";
#	    print STDERR "$p5flank\n";
#	    if($p5tsd eq "NA"){print STDERR "\n";}else{print STDERR "$p5tsd\n"};
#	    if($go[20] eq "NA"){print STDERR "\n";}else{print STDERR "$go[20]\n"};
#	    print STDERR "$reseq\n";
#	    if($go[22] eq "NA"){print STDERR "\n";}else{print STDERR "$go[22]\n"};
#	    if ($p3tsd eq "NA"){print STDERR "\n";}else{print STDERR "$p3tsd\n"};
#	    print STDERR "$p3flank\n";
#	}
    }elsif($opt{O} == 2){
	if ($RIMDlen >=  30){
	    $outline = join("|",@go[0..23],$RIMD,$RIMDlen);
	    print ">$outline\n$tmp2\n";
#	    exit;
	}
    }

#    exit;
    unlink($tmpfile);
    unlink($bl2seqoutfile);
    unlink($tmpreferrence);
#    exit;
}
close(IN);#read blat output files and store the sequence;
#print STDERR "$errorcount entries may suffer from sequences divergence in the out-group genomes and may carry $previousfalsenegative entries of transduction seequences\n";
