#!/usr/bin/perl -w
use strict;

my $file = shift;	# this is expected to be the DESeq2 output
my $out = shift;	# the output file

open my $OUT, ">$out" or die "Cannot open file: $!\n";
open my $IN, "<$file" or die "Cannot open file: $!\n";
<$IN>;	# skip the header
print $OUT "Gene_ID\tMolecule_Type\tMean_normalized_read_count\tlog2FoldChange_Treatment-vs-Control\tp-value\tq-value\tComplete_Info\n";
my @content;
while(<$IN>)	{
  chomp;
  my @decom = split /\s+/, $_;
  $decom[0] =~ s/\"//g;
  my @decom2 = split /\|/, $decom[0];
  push @content, [$decom2[5], $decom2[7], $decom[1], $decom[2], $decom[5], $decom[6], $decom[0]];
}
close $IN;

# sort the array based on q-values
@content = sort {$a->[4] <=> $b->[4]} @content;
my $i; my $j;
for($i = 0; $i < scalar(@content); ++ $i)	{
  my $k = scalar(@{$content[$i]});
  for($j = 0; $j < $k - 1; ++ $j)	{
    print $OUT "$content[$i][$j]\t";
  }
  print $OUT "$content[$i][$k - 1]\n";
}
close $OUT;
