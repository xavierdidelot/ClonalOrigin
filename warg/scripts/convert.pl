#!/usr/bin/perl
#Converts an output file from the old format (December) to the new format(January)

if ($#ARGV!=1) {die("Usage: convert.pl [input] [output]\n");}

open($file,$ARGV[0]) or die("Can't open file ",$ARGV[0],"\n");
open($fileout,">",$ARGV[1]) or die("Can't open output file\n");

while ($buf=<$file>) {
if ( $buf =~ /<Tree>/ ) {
$tree=<$file>;
$buf=<$file>;
} elsif ( $buf =~ /<Iteration>/ ) {
print $fileout $buf."<Tree>\n".$tree."</Tree>\n";
} else {
print $fileout $buf;
}
}
close($file);
close($fileout);
