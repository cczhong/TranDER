#!/usr/bin/perl -w
use strict;

my $input = shift;    # FASTQC data file (fastqc_data.txt)
my $output = shift;   # FASTA file that can be used by trimmomatic

die "Usage: perl GetOverrepresentedSeqs [expecting fastqc_data.txt] [expecting output file]\n" if (!defined $input || !defined $output);
open my $IN, "<$input" or die "Cannot read input file $input: $!\n";
my $index = 0;

# parse the input; extract "Overrepresented sequences" and "Kmer content"
my @content;
while(<$IN>)  {
  chomp;
  push @content, $_;
}
close $IN;


# generate the sequences
open my $OUT, ">$output" or die "Cannot create output file $output: $!\n";
my $i;
my $j;
for($i = 0; $i < scalar(@content); ++ $i) {
  if($content[$i] =~ /^\>\>Overrepresented\s+sequences\s+fail/)  {
    for($j = $i + 1; $j < scalar(@content); ++ $j) {
      next if $content[$j] =~ /^\#/;
      last if $content[$j] =~ /^\>\>END\_MODULE/;
      my @decom = split /\s+/, $content[$j];
      print $OUT ">index_$index\n$decom[0]\n";
      ++ $index;
    }
    $i = $j;
  } elsif($content[$i] =~ /^\>\>Kmer\s+Content\s+fail/) {
    for($j = $i + 1; $j < scalar(@content); ++ $j) {
      next if $content[$j] =~ /^\#/;
      last if $content[$j] =~ /^\>\>END\_MODULE/;
      my @decom = split /\s+/, $content[$j];
      print $OUT ">index_$index\n$decom[0]\n";
      ++ $index;
    }
    $i = $j;
  }
}
close $OUT;
