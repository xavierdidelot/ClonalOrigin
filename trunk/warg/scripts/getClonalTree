#!/usr/bin/perl -w
use strict;

if(@ARGV!=2){
	print STDERR "Usage: getClonalTree <ClonalFrame file> <output tree file>\n\n";
	exit -1;
}
my $src = $ARGV[0];
my $dest = $ARGV[1];
my $cl = "head -n2 $src | tail -n1 > $dest";
#print "$cl\n";
`$cl`;

