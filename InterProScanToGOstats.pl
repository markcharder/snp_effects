#!/usr/bin/env perl

# This script parses output from interproscan and snpEff to generate a file that is useable with the GOstats bioconductor package.

use strict;
use warnings;
use Getopt::Long;

my $interpro = "";
my $output = "";
my $help;

GetOptions('interpro-file|i=s' => \$interpro, 'prefix|p:s' => \$output);

my $usage = <<USAGE;


	$0 --interpro-file <interpro file> [--prefix <output prefix>]

	-	Options can be abbreviated to '-i' for interpro file or '-p' for prefix.
	-	Output defaults to 'out.txt'.
	-	This script creates a csv suitable for use with the GOstats package in R.


USAGE

if ($output eq ""){
	$output = "out.txt";
}

open my $FH, "<", $interpro or die ("\n$usage\n");

while (my $line = <$FH>){

	chomp $line;
	my @fields = split /\t/, $line;
	my $id = $fields[0];
	my $goterms = $fields[-1];
	my @goterms = split /\|/, $goterms;
	open my $OFH, ">>", $output or die ("\nOutput file does not seem to be writable.\n");
	if ($. == 1){
		print $OFH "GO_TERM\tEVIDENCE\tGENE_ID\n";
	}
	foreach my $term (@goterms){
		if ($term =~ /GO:\d+$/){
			print $OFH "$term\tISS\t$id\n";
		}
	}
	close $OFH;

}

close $FH;
