#!/usr/bin/perl
#This script edits the Newick string in the input file by rescaling all distances according to the factor given in the first argument

@ARGV==3 || die("Usage: rescale [SCALE] [INPUT] [OUTPUT]");

$scale=$ARGV[0];
open($file,$ARGV[1]) or die("Can't open file ",$ARGV[1],"\n");
$string=<$file>;
close($file);

while ($string =~ /:([0-9]+\.[0-9]+)/) {
$rep=$1*$scale;
$string =~ s/:([0-9]+\.[0-9]+)/!$rep/;
}

$string =~ s/!/:/g;

open($file,">",$ARGV[2]) or die("Can't open output file\n");
print $file "$string\n";
close($file);