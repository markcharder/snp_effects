#!/usr/bin/perl

# This script was originially packaged with snpEff but has been modified (to be used with a current Perl version).
# It prints one gene with one effect per line from a vcf file output from a first pass snpEff run.

 open my $STDIN, "<", $ARGV[0];
 while( $l = <$STDIN> ) {
         chomp $l;
         if( $l !~ /^#/ ) {
                 ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info) = split /\t/, $l;
                 print "$chrom\t$pos\t$id\t$ref\t$alt\t$qual\t$filter\n";  

                 foreach $in ( split /;/, $info ) {
                         if( $in =~ /^(.*?)=(.*)/ ) {
                                 $key = $1;
                                 $values = $2; 

                                 @vals = split /,/, $values;

                                 if( $#vals > 0 ) {
                                         foreach $val ( @vals ) { print "\t\t\t\t\t\t\t$key\t$val\n"; }
                                 } else {
                                         print "\t\t\t\t\t\t\t$key\t$values\n";
                                 }
                         } else {
                                 print "\t\t\t\t\t\t\t$in\n";
                         }
                 }
         }
 }
 close $STDIN;
