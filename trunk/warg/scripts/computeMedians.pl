#!/usr/bin/perl -w
use strict;
use IO::Uncompress::Bunzip2 qw(bunzip2 $Bunzip2Error);
use XML::Parser;

if(@ARGV==0){
	die "Usage: computeMedians.pl <ClonalOrigins XML or xml.bz2>\n";
}

my @lens;
my @meantheta;
my @meandelta;
my @meanrho;
my $itercount=0;
my $curtheta=0;
my $curdelta=0;
my $currho=0;
my $tag;

my $blockcount=scalar(@ARGV);	# assume one block per file

# extract posterior mean estimates of global parameters from each file
foreach my $f (@ARGV){
	my $fs;
	if($f =~ /\.bz2$/){
		$fs = bunzip2 $f => "tmpxml" or die "IO::Uncompress::Bunzip2 failed: $Bunzip2Error\n";
		$fs = "tmpxml";
	}else{
		$fs = $f;
	}
	my $parser = new XML::Parser();

	$parser->setHandlers(      Start => \&startElement,
                           End => \&endElement,
                           Char => \&characterData,
                           Default => \&default);

	$itercount=0;
	$curtheta=0;
	$curdelta=0;
	$currho=0;
	my $doc;
	eval{ $doc = $parser->parsefile($fs)};
	print "Unable to parse XML of $f, error $@\n" if $@;
	next if $@;
	print "parsed $f\n";
	$curtheta /= $itercount;
	$curdelta /= $itercount;
	$currho /= $itercount;
	push( @meantheta, $curtheta );
	push( @meandelta, $curdelta );
	push( @meanrho, $currho );
}

# convert to per-site values of theta and rho
for( my $i=0; $i<@meantheta; $i++){
	$meantheta[$i] /= $lens[$i];
	$meanrho[$i] /= $meandelta[$i] + $lens[$i];
}

# now compute a weighted median
my %thetalens;
my %deltalens;
my %rholens;
my $lensum=0;
for( my $i=0; $i<@meantheta; $i++){
	$thetalens{$meantheta[$i]}=$lens[$i];
	$deltalens{$meandelta[$i]}=$lens[$i];
	$rholens{$meanrho[$i]}=$lens[$i];
	$lensum += $lens[$i];
}
print "lensum is $lensum\n";


my @tsort = sort{ $a <=> $b } @meantheta;
my @dsort = sort{ $a <=> $b } @meandelta;
my @rsort = sort{ $a <=> $b } @meanrho;
my $j=0;
for(my $ttally=$thetalens{$tsort[$j]}; $ttally < $lensum/2; $ttally += $thetalens{$tsort[$j]})
{
	$j++;
}
print "Median theta: ".$tsort[$j]."\n";

$j=0;
for(my $dtally=$deltalens{$dsort[$j]}; $dtally < $lensum/2; $dtally += $deltalens{$dsort[$j]})
{
	$j++;
}
print "Median delta: ".$dsort[$j]."\n";

$j=0;
for(my $rtally=$rholens{$rsort[$j]}; $rtally < $lensum/2; $rtally += $rholens{$rsort[$j]})
{
	$j++;
}
print "Median rho: ".$rsort[$j]."\n";


exit;

sub startElement {
       my( $parseinst, $element, %attrs ) = @_;
	$tag = $element;
       SWITCH: {
              if ($element eq "Iteration") {
                     $itercount++;
                     last SWITCH;
              }
              if ($element eq "delta") {
                     last SWITCH;
              }
              if ($element eq "rho") {
                     last SWITCH;
              }
       }
}

sub endElement {
	$tag = "";
}

sub characterData {
       my( $parseinst, $data ) = @_;
	$data =~ s/\n|\t//g;
	$curtheta += $data if ($tag eq "theta");
	$curdelta += $data if ($tag eq "delta");
	$currho += $data if ($tag eq "rho");
	if($tag eq "Blocks"){
		$data =~ s/.+\,//g;
		push( @lens, $data ) if(length($data)>1);
	}
}

sub default {
}
