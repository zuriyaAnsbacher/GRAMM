#!/usr/local/bin/perl
##############################################################################
# SCRIPT NAME:	ds
# DESCRIPTION:	convert dna_ser.txt to ser_dna.txt (reads rel_ser_ser.txt)
# PARAMETERS:	none
# OUTPUT:	
# TABLES:       
#
# DATE WRITTEN: 2012-08-16
# WRITTEN BY:   Martin J. Maiers
#
# REVISION HISTORY: 
# REVISION DATE		REVISED BY	DESCRIPTION 
# ------- ----------	--------------	-------------------------------------
#
#       COPYRIGHT (C) 2011 NATIONAL MARROW DONOR PROGRAM.  
#               ALL RIGHTS RESERVED        
##############################################################################
use strict;    # always
use warnings;  # or else
# use lib "/MDP/prod/research/lib/perl";
# use Connect;

my %S;


#
# load rel_ser_ser.txt
#
my %B;
my $file_ss = "rel_ser_ser.txt";
open FILE_SS, $file_ss or die "$!: $file_ss";
while(<FILE_SS>) {
  chomp;
  tr///d;	# windoze
  my ($loc, $ser, $spl, $ass) = split /\;/;
  foreach ($spl, $ass) {
    next unless defined;
    foreach (split /\//) {
      $B{$loc}{$_}=$ser;  # assocate splits with their broads
    }
  }
}

#
# load rel_dna_ser.txt
#
my $file = "rel_dna_ser.txt";
open FILE, $file or die "$!: $file";
while (<FILE>) {
  chomp;
  tr///d;
  my ($loc, $allele, $unambig_ser, $possible_ser, $assumed_ser, $expert_ser) =
	split /\;/;
  my $serloc = $loc;
  $serloc=~tr/\*//d; # really?
  $serloc = "DR" if $serloc eq "DRB1";

  my $ct =0;
  $ct++ if defined $unambig_ser && $unambig_ser=~/\S/;
  $ct++ if defined $possible_ser && $possible_ser=~/\S/;
  $ct++ if defined $assumed_ser && $assumed_ser=~/\S/;
  print "multiple definitions for: $loc $allele" if $ct >1;  #shouldn't happen
  foreach ($unambig_ser, $possible_ser, $assumed_ser) {
    next unless defined;
    foreach (split /\//) {
      $S{$serloc}{$_}{$loc}{$allele}++;
      if (defined $B{$serloc}{$_}) {
        my $b = $B{$serloc}{$_};
        $S{$serloc}{$b}{$loc}{$allele}++;
      }
    }
  }
}

my $file_sd = "rel_ser_dna.txt";
open FILE_SD, ">$file_sd" or die "$!: $file_sd";

foreach my $serloc (sort keys %S) {
  foreach my $ser (sort keys %{$S{$serloc}}) {
    foreach my $loc (sort keys %{$S{$serloc}{$ser}}) {
      my $al = join '/', sort keys %{$S{$serloc}{$ser}{$loc}};
      print FILE_SD join ('	', $serloc, $ser, $loc, $al), "\n";
    }
  }
}
exit 0;

