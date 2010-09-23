#!/usr/bin/perl -w
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error) ;

use strict;

my $maxfile = -1;
for( my $i=0; $i < (@ARGV) - 1; $i++){
	my $fI = $ARGV[$i];
	$fI =~ s/.*\.(\d+)\.xml.bz2/$1/g;
	$maxfile = $fI if($fI>$maxfile);
}
print "Found $maxfile xml block files\n";

my $base = $ARGV[0];
$base =~ s/(.*)\.\d+\.xml.bz2/$1/g;

my $zout = new IO::Compress::Bzip2 "$base.mauveConcat.bz2" or die "bzip2 failed: $Bzip2Error\n";
$zout->print("<?xml version = '1.0' encoding = 'UTF-8'?>\n");
$zout->print("<weakArgData>\n");
for( my $i=1; $i<=$maxfile; $i++ ){
	unless( -e "$base.$i.xml.bz2" ){
		# output file doesn't exist, so create an empty block placeholder
		$zout->print("<outputFile>\n");
		$zout->print("<Blocks>\n");
		$zout->print("</Blocks>\n");
		$zout->print("<comment>Empty placeholder for missing block $i</comment>\n");
		$zout->print("</outputFile>\n");
		next;
	}
	open( CURFILE, "bzcat $base.$i.xml.bz2 |" );
	while( my $line = <CURFILE> ){
		next if substr($line, 0, 5) eq "<?xml";
		$zout->print($line);
	}
	close CURFILE;
}
$zout->print("</weakArgData>\n");
$zout->flush();
$zout->close();
