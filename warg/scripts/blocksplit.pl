#!/usr/bin/perl -w 
use strict;

my $fname = $ARGV[0];

open(INFILE, "$fname") || die "Unable to open input XMFA file $fname\n";
my $i = 1;
my $curoutfile = $fname.".1";
open(OUTFILE, ">$curoutfile");
while( my $line = <INFILE> )
{
	print OUTFILE $line;
	if($line =~ /=/)
	{
		close OUTFILE;
		$i++;
		$curoutfile = $fname.".$i";
		open(OUTFILE, ">$curoutfile");
	}
}
close OUTFILE;
# last output file is extra
`rm $curoutfile`;

