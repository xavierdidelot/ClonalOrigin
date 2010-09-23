#!/usr/bin/perl -w
# Program to parse bzip2'ed output from a set of warg runs done gene-by-gene 
# Run from within the directory containing the .bz2 files
# Generates text summary files that can be plotted using genomeplots_log.R
# Usage: rhomulus.pl <warg output base filename> <# of genes> <file with mapping from warg run ID to gene ID> <path to alignment files> <file containing gene coordinates in the reference genome>
#
# @author Aaron Darling
# @license GPL
use strict;
my $basename = $ARGV[0];
my $maxgenes = $ARGV[1];
my $decoder = $ARGV[2];
my $genpath = $ARGV[3];
my $coords = $ARGV[4];
die unless @ARGV==5;

my @stats=("rho","theta","delta","tmrca");
# read in gene ID decoder
my %geneId;
open( DECOD, "$decoder" ) || die "unable to open decoder $decoder\n";
while( my $line = <DECOD> )
{
	chop($line);
	my @stuff = split( /:/, $line );
	$geneId{$stuff[0]} = $stuff[1];
}

open(RHO, ">rho_summary.txt");
open(THETA, ">theta_summary.txt");
open(TMRCA, ">tmrca_summary.txt");
open(DELTA, ">delta_summary.txt");
open(RHOTHETA, ">rhotheta_summary.txt");
open(RHOPERSITE, ">rhopersite_summary.txt");
open(RHODELTA, ">rhodelta_summary.txt");
open(ROVERM, ">roverm_summary.txt");

# read in gene position table
my %lends;
my %rends;
my %strands;
open( COORD, "$coords") || die "Unable to open $coords\n";
while( my $line = <COORD> )
{
	my @stuff = split( /\s+/, $line );
	$lends{$stuff[4]} = $stuff[0];
	$rends{$stuff[4]} = $stuff[1];
	$strands{$stuff[4]} = $stuff[2];
}

for( my $gI = 1; $gI < $maxgenes; $gI++ )
{
	my $bzfname = "$basename.$gI.xml.bz2";
	print STDERR "Trying $bzfname\n";
	next unless -e $bzfname;

	# look up position in the genome
	my $alnfile = $genpath."/".$geneId{$gI};
	open( ALN, $alnfile ) || die "Error opening $alnfile\n";
	my $geneid = 0;
	while( my $line = <ALN> )
	{
		next unless $line =~ /(ESCO001c01_\d+)/;
		$geneid = $1;
		last;
	}

	# read posterior file
	open( POST, "bzcat $bzfname |" ) || next;
	my @rho;
	my @theta;
	my @tmrca;
	my @delta;
	my @redgelens;
	my @curredges;
	my $prevredge = 0;	# whether the previous line had a recedge
	while( my $line = <POST> )
	{
		checkLine($line,"rho",\@rho);
		checkLine($line,"theta",\@theta);
		checkLine($line,"tmrca",\@tmrca);
		checkLine($line,"delta",\@delta);
		my $tmp = checkRecEdge($line,\@curredges);
		if($tmp == 0 && $prevredge == 1)
		{
			push( @redgelens, [ @curredges ] );
			@curredges = ();
		}
		$prevredge = $tmp;
	}
	# file is parsed, summarize stats
	print RHO summaryLine($geneid,\@rho);
	print THETA summaryLine($geneid,\@theta);
	print TMRCA summaryLine($geneid,\@tmrca);
	print DELTA summaryLine($geneid,\@delta);
	# compute rho over thetas
	my @rhotheta = @rho;
	for(my $j = 0; $j<@rhotheta;$j++){	$rhotheta[$j]/=$theta[$j];}
	print RHOTHETA summaryLine($geneid,\@rhotheta);
	my @rhopersite = @rho;
	my $L = $rends{$geneid}-$lends{$geneid};
	for(my $j = 0; $j<@rhopersite;$j++){	$rhopersite[$j]/=($delta[$j]+$L-1);}
	print RHOPERSITE summaryLine($geneid,\@rhopersite);
	my @rhodelta = @rhopersite;
	for(my $j = 0; $j<@rhopersite;$j++){	$rhodelta[$j]*=$delta[$j];}
	print RHODELTA summaryLine($geneid,\@rhodelta);
	# calculate Xavier's Y for r/m
	my @roverm = @rhodelta;
	for(my $j=0; $j<@redgelens; $j++){
		my $Y = 0;
		my $cur = $redgelens[$j];
		for( my $k=0; $k<@$cur; $k++)
		{
			$Y += 0.75*(1.0-exp(-(4.0/3.0)*$theta[$j]*@$cur[$k]/$L));
		}
		$Y /= @$cur;
		$roverm[$j] *= $Y;
	}
	print ROVERM summaryLine($geneid,\@roverm);
}


sub checkLine {
	my $linea = shift;
	my $stat = shift;
	my $statts = shift;
	if( $linea =~ /$stat\>(\d+)\.(\d+)\<\/$stat/ )
	{
		my $rv = $2;
		$rv /= 10 while( $rv > 1 );
		$rv += $1;
		push( @$statts, $rv );
	}elsif(  $linea =~ /$stat\>(\d+)\<\/$stat/ )
	{
		push( @$statts, $1 );
	}
}

sub summaryLine {
	my $gid = shift;
	my $tmp = shift;
	my @stats = @$tmp;
	my $sum = 0;
	foreach( @stats ){
		$sum += $_;
	}
	my $mean = $sum/@stats;
	my @sorted = sort { $a <=> $b } @stats;
	
	my $sd = 0;
	foreach( @stats ){
		$sd += ($mean-$_)*($mean-$_);
	}
	$sd /= @stats;
	$sd = sqrt($sd);

	my $l = $lends{$gid};
	my $r = $rends{$gid};
	my $s = $strands{$gid};
	my $liner = "$gid\t$mean\t$sd\t$l\t$r\t$s\t";
	$liner .= $sorted[0]."\t".$sorted[@sorted*0.1]."\t".$sorted[@sorted*0.25]."\t".$sorted[@sorted*.5]."\t".$sorted[@sorted*.75]."\t".$sorted[@sorted*.9]."\t".$sorted[@sorted-1]."\n";
	return $liner;
}

sub checkRecEdge {
	my $linea = shift;
	my $redgelens = shift;
	return 0 unless( $linea =~ /\<recedge\>/);
	if( $linea =~ /\<afrom\>(.+)\<\/afrom\>\<ato\>(.+)\<\/ato\>/ )
	{
		push( @$redgelens, abs($1-$2) );
	}	
	return 1;
}
