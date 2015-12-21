#!/usr/bin/perl -w
# Author: Thomas Thiel
# Program name: prim_output.pl
# Description: converts the Primer3 output into an table

# Modified my Matt Gitzendanner (MAG), University of Florida
# MAG 03/18/15 Modified for new tags for Primer3
# MAG 03/18/15 Modified to not print anything for loci with no primers.
# MAG 03/18/15 Modified to print file of sequence names that can be used to pull those from the original fasta file.

open (SRC,"<$ARGV[0]") || die ("\nError: Couldn't open Primer3 results file (*.p3out) !\n\n");
my $filename = $ARGV[0];
$filename =~ s/\.p3out//;
open (IN,"<$ARGV[1]") || die ("\nError: Couldn't open source file containing MISA (*.misa) results !\n\n");
open (OUT,">$filename.results") || die ("nError: Couldn't create file !\n\n");
open (GOOD,">$filename.goodSeqs") || die ("nError: Couldn't create file !\n\n");  #MAG: Added line to print good sequences.

my ($seq_names_failed,$count,$countfailed);

print OUT "ID\tSSR nr.\tSSR type\tSSR\tsize\tstart\tend\t";
print OUT "FORWARD PRIMER1 (5'-3')\tTm(°C)\tsize\tREVERSE PRIMER1 (5'-3')\tTm(°C)\tsize\tPRODUCT1 size (bp)\tstart (bp)\tend (bp)\t";
print OUT "FORWARD PRIMER2 (5'-3')\tTm(°C)\tsize\tREVERSE PRIMER2 (5'-3')\tTm(°C)\tsize\tPRODUCT2 size (bp)\tstart (bp)\tend (bp)\t";
print OUT "FORWARD PRIMER3 (5'-3')\tTm(°C)\tsize\tREVERSE PRIMER3 (5'-3')\tTm(°C)\tsize\tPRODUCT3 size (bp)\tstart (bp)\tend (bp)\n";

undef $/;
my $in = <IN>;
study $in;

$count=0;  #MAG: Added line to count primers created.
$/ = "=\n";

while (<SRC>)
  {
  my ($id,$ssr_nr) = (/SEQUENCE_ID=(\S+)_(\d+)/);
  
  $in =~ /($id\t$ssr_nr\t.*)\n/;
  my $misa = $1;
  
  /PRIMER_LEFT_0_SEQUENCE=(.*)/ || do {$count_failed++; next}; #MAG: Cut print OUT "$misa\n"; from within {} to not print anything for failed primer designs.
  my $info = "$1\t";
  
  /PRIMER_LEFT_0_TM=(.*)/; $info .= "$1\t";  #MAG: Updated tags on this line and following for primer3 v2+
  /PRIMER_LEFT_0=\d+,(\d+)/; $info .= "$1\t";
  
  /PRIMER_RIGHT_0_SEQUENCE=(.*)/;  $info .= "$1\t";
  /PRIMER_RIGHT_0_TM=(.*)/; $info .= "$1\t";
  /PRIMER_RIGHT_0=\d+,(\d+)/; $info .= "$1\t";
  
  /PRIMER_PAIR_0_PRODUCT_SIZE=(.*)/; $info .= "$1\t";
  /PRIMER_LEFT_0=(\d+),\d+/; $info .= "$1\t";
  /PRIMER_RIGHT_0=(\d+),\d+/; $info .= "$1\t";
  
  
  /PRIMER_LEFT_1_SEQUENCE=(.*)/; $info .= "$1\t";
  /PRIMER_LEFT_1_TM=(.*)/; $info .= "$1\t";
  /PRIMER_LEFT_1=\d+,(\d+)/; $info .= "$1\t";
    
  /PRIMER_RIGHT_1_SEQUENCE=(.*)/;  $info .= "$1\t";
  /PRIMER_RIGHT_1_TM=(.*)/; $info .= "$1\t";
  /PRIMER_RIGHT_1=\d+,(\d+)/; $info .= "$1\t";
  
  /PRIMER_PRODUCT_SIZE_1=(.*)/; $info .= "$1\t";
  /PRIMER_LEFT_1=(\d+),\d+/; $info .= "$1\t";
  /PRIMER_RIGHT_1=(\d+),\d+/; $info .= "$1\t";
  
  
  /PRIMER_LEFT_2_SEQUENCE=(.*)/; $info .= "$1\t";
  /PRIMER_LEFT_2_TM=(.*)/; $info .= "$1\t";
  /PRIMER_LEFT_2=\d+,(\d+)/; $info .= "$1\t";
    
  /PRIMER_RIGHT_2_SEQUENCE=(.*)/;  $info .= "$1\t";
  /PRIMER_RIGHT_2_TM=(.*)/; $info .= "$1\t";
  /PRIMER_RIGHT_2=\d+,(\d+)/; $info .= "$1\t";
  
  /PRIMER_PRODUCT_SIZE_2=(.*)/; $info .= "$1\t";
  /PRIMER_LEFT_2=(\d+),\d+/; $info .= "$1\t";
  /PRIMER_RIGHT_2=(\d+),\d+/; $info .= "$1";
  
  $count++;
  print OUT "$misa\t$info\n";
  print GOOD "$id\n";  #MAG: Added printing of sequence name to file for subsequent steps.
  
  };

print "\nPrimer modelling was successful for $count sequences.\n";
print "Primer modelling failed for $count_failed sequences.\n";
