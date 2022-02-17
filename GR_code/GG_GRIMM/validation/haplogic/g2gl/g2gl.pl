#!/usr/bin/env perl
##############################################################################
# SCRIPT NAME:	g2gl.pl
# DESCRIPTION:	genotype to genotype-list
#
# DATE WRITTEN: 2017-02-08
# WRITTEN BY:   Martin Maiers
#
##############################################################################
use strict;    # always
use warnings;  # or else
use lib '.';
use lib '../g2gl';
use MAC;
use Getopt::Std;
use REST::Client;
use JSON;
use JSON::Parse 'parse_json';
use MIME::Base64;
use Math::Round;

use Memoize;
memoize('expandTyp');

#
# setup neo4j client
# TODO: get config from a file
#
my $service_url = 'http://localhost:7474'; # location of the service
my $service_up = encode_base64(join (':', "neo4j", "ontological"));
my $client = REST::Client->new({ host => $service_url,});
$client->addHeader('Content-Type', 'application/json;charset=UTF-8');
$client->addHeader('Accept', 'application/json');
$client->addHeader('Authorization', "Basic $service_up");


my %opts;
getopts('i:o:f:', \%opts);
die "$0: -i inputfile -o outputfile -f freqfile" unless defined $opts{i} && defined $opts{o};



# set to absolute path of your installation
my $top = "../..";

# set to HF file
my $hf_file = "$top/graph_generator/data/wmda/freqs.txt";

# set to input genotypes file
my $gfile = $opts{i};
my $ofile = $opts{o};


my %A;  # hash table of alleles

if (defined $opts{f}) {
  # get frequencies from a file
  $hf_file = $opts{f};
  ##############################################################################
  # parse file of haplotype frequencies
  ##############################################################################
  open HF_FILE, "cat $hf_file|" or die "$!: $hf_file";
  while(<HF_FILE>) {
    chomp;
    #A*01:01g~C*01:02g~B*15:01g~DRB1*01:01~DQB1*05:01;6e-05
    my ($hap, $freq) = split /\;/;
    foreach my $allele (split /\~/, $hap) {
      $allele=~s/g//;      # remove little-g
      my $who_allele = "HLA-".$allele;
      $A{$who_allele}+=$freq;  # note that this allele exists; sum freq just because
    }
  }
  close HF_FILE;
} else {
  ##############################################################################
  # get haplotypes from Neo4j
  ##############################################################################
  my @loci = qw/A B C R Q/; # TODO: get from an info node in the graph
  foreach my $loc (@loci) {
    my $cypher = qq/match (n:$loc) return n.name, n.frequency/;
  
    foreach my $row (doCypher($client, $cypher)) {
      my ($allele, $freq) = @{$row->{row}};
      my $who_allele = "HLA-".$allele;
      $A{$who_allele}++;  # note that this allele exists 
      # this will include all alleles that exist in the frequencies
      # TODO: sum the frequency of the allele to allow a frequency threshold
    }
  }
}



##############################################################################
# parse input file
##############################################################################
open OFILE, ">$ofile" or die "$!: $ofile";
open GFILE, $gfile or die "$!: $gfile";
while(<GFILE>) {
  chomp;
  #D000001%A*11:XX+A*02:XX^C*UUUU+C*UUUU^B*13:XX+B*44:XX^DRB1*07:XX+DRB1*01:XX^DQB1*UUUU+DQB1*UUUU
  my ($id, $g) = split /\%/;
  my @gl_array = ();
  foreach my $loctyp (split /\^/, $g) {
    # within a locus

    # look for | characters
    # and if you find it, refold the glstring to be a single genotype with
    # allelic ambiguity
    my ($t1, $t2);
    if ($loctyp=~/\|/) {
      # e.g. B*40:02+B*42:01|B*40:40+B*42:02
      # expand genotype ambiguity within the locus
      my %typ1; 
      my %typ2; 
      foreach my $geno (split /\|/, $loctyp) {
        my ($p1, $p2) = split /\+/, $geno;
        $typ1{$p1}++; $typ2{$p2}++;
      }
      $t1 = join ('/', sort keys %typ1);
      $t2 = join ('/', sort keys %typ2);
    } else {
      ($t1, $t2) = split /\+/, $loctyp;
    }
    my ($l1, $a1) = split /\*/, $t1;
    my ($l2, $a2) = split /\*/, $t2;

    if ($l1 ne $l2) {
      warn "mismatching loci: $l1 vs $l2 for id: $id";
      next;
    }

    if ($a1 eq "UUUU" || $a2 eq "UUUU") {
      # if either typing is UUUU then skip this loctyp
      next; 
    }

    # expand
    my @al1 = expandTyp($id, $t1);
    my @al2 = expandTyp($id, $t2);
    
    # reduce
    my $gl = join ('+', 
      sort  join ('/', reduceTyp(@al1)), join ('/', reduceTyp(@al2)));
      
    $gl=~s/HLA\-//g;
    push @gl_array, $gl;
  }
  print OFILE join ('%', $id, join ('^', @gl_array)), "\n";
}
close OFILE;
close GFILE;

exit 0;


#
# use the alleles in the frequency set to reduce the ambiguity list
#
sub reduceTyp {
  my (@al) = @_;
  my @ret;
  
  foreach (@_) {
    push @ret, $_ if defined $A{$_};
  }
  if (@ret) {
    return @ret;
  } else {
    # at this point everything has been reduced; so pass back the original list
    return @al;
  }
}

# 
#  expand and process HLA types
#  - convert allele lists to use HLA-A instead of A
#  - convert allele lists to remove little g
#  - convert allele codes (A*02:AB, A*02:XX)
# 
sub expandTyp {
  my ($id, $t) = @_;

  my @al;
  
  if ($t=~/\//) {
    # its an allele ambiguity already
    foreach my $a (split /\//, $t) {
      $a=~s/g//;      # remove little-g
      my $who_typ = "HLA-".$a;
      push @al, $who_typ;
    }
  } else {
    # perform MAC expansion
    my $who_typ = "HLA-".$t;

    # expand
    my $rc = MAC::decode($who_typ, \@al);
    if ($rc ne 200) {
      warn "decode failed with code: $rc for $who_typ for id: $id";
      next;
    }
  }
  return @al;
}

sub doCypher {
  my ($client, $query) = @_;
  my $request = { statements => [ { statement => $query } ] };
  my $json_request = JSON::to_json($request);
  $client->POST('/db/data/transaction/commit', $json_request, {});
  if ( $client->responseCode() eq '500') {
     die $client->responseContent();
  }
  my $json_response = $client->responseContent;
  my $response = parse_json($json_response);
  return undef unless defined $$response{results}[0]{data};
  return @{$$response{results}[0]{data}};
}

