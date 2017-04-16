#!/usr/bin/perl -w
use strict;
use Getopt::Long;

# TranDER: Transcript-level Differential Expression analysis on RNA-seq data

my $manifest;   # the manifest specifying the locations of the programs and the data and the groups
my $work_dir;   # the working directory for storing temporary files and outputs

# loading the arguments

GetOptions (
  "manifest=s" => \$manifest,
  "dir=s" => \$work_dir
) or die("Error in command line arguments\n");

if(!defined $manifest || !defined $work_dir)  {
  print "TranDER: Transcript-level Differential Expression analysis on RNA-seq data\n";
  print "Usage: perl TranDER.pl --manifest=[MANIFEST_FILE] --dir=[WORK_DIRECTORY]\n";
  print " --manifest:\tthe manifest file\n";
  print " --dir:\t\tthe working directory\n";
  exit;
}

# create file directory
print "!!!  TranDER: preparing and checking run environment...\n";
mkdir "$work_dir" or die "FATAL ERROR: fail to create working directoy!" if !(-e "$work_dir");
mkdir "$work_dir/Trimmed_Reads" or die "FATAL ERROR: fail to create working directoy!" if !(-e "$work_dir/Trimmed_Reads");
mkdir "$work_dir/FASTQC_Report" or die "FATAL ERROR: fail to create working directoy!" if !(-e "$work_dir/FASTQC_Report");
mkdir "$work_dir/Intermediate" or die "FATAL ERROR: fail to create working directoy!" if !(-e "$work_dir/Intermediate");
mkdir "$work_dir/Mapping" or die "FATAL ERROR: fail to create working directoy!" if !(-e "$work_dir/Mapping");
mkdir "$work_dir/Transcript_abundance" or die "FATAL ERROR: fail to create working directoy!" if !(-e "$work_dir/Transcript_abundance");
mkdir "$work_dir/Differential_expression" or die "FATAL ERROR: fail to create working directoy!" if !(-e "$work_dir/Differential_expression");
mkdir "$work_dir/Report" or die "FATAL ERROR: fail to create working directoy!" if !(-e "$work_dir/Report");
mkdir "$work_dir/Log" or die "FATAL ERROR: fail to create working directoy!" if !(-e "$work_dir/Log");

# read the manifest file
my @treatments;       # a list of runs that are going to be treated as the treatment group
my @treatment_IDs;    # the identifier (without path and file extension) of the treatment group runs
my @controls;         # a list of runs that are going to be treated as the control group
my @control_IDs;      # the identifier (without path and file extension) of the control group runs
my $HOME;             # the home directory of TranDER
my $identifier;       # the identifier for the project
my $ref_prefix;       # the reference prefix (assume that the HISAT2 index is done)
my $ref_transcript;   # the reference transcriptome (expecting multi-FASTA file)
my $exe_fastqc;       # the executable of the FASTQC program
my $exe_trimmomatic;  # the executable of the Trimmomatic program
my $exe_hisat2;       # the executable of the HISAT2 program
my $exe_express;      # the executable of the eXpress program
my $exe_Rscript;      # the executable of the eXpress program
my $num_threads;      # the number of threads allowed
# TODO: other variables are to be added when extending the funtionality

#==================== Loading the manifest =====================================
open my $IN, "<$manifest" or die "Cannot read manifest file: $!\n";
while(<$IN>) {
  chomp;
  next if /^\#/;
  my @decom = split /\=/, $_;
  @treatments = split /\;/, $decom[1] if $decom[0] eq 'Treatment_group';
  @controls = split /\;/, $decom[1] if $decom[0] eq 'Control_group';
  $HOME = $decom[1] if $decom[0] eq 'HOME';
  $identifier = $decom[1] if $decom[0] eq 'Project_ID';
  $exe_fastqc = $decom[1] if $decom[0] eq 'exe_FASTQC';
  $exe_trimmomatic = $decom[1] if $decom[0] eq 'exe_TRIMMOMATIC';
  $exe_hisat2 = $decom[1] if $decom[0] eq 'exe_HISAT2';
  $exe_express = $decom[1] if $decom[0] eq 'exe_EXPRESS';
  $exe_Rscript = $decom[1] if $decom[0] eq 'exe_RSCRIPT';
  $ref_prefix = $decom[1] if $decom[0] eq 'Reference_prefix';
  $ref_transcript = $decom[1] if $decom[0] eq 'Reference_transcript';
  $num_threads = $decom[1] if $decom[0] eq 'Num_threads';
  # TODO: to be added
}
close $IN;

#==================== Checking whether all files exist =========================
# TODO

#==================== Parsing the IDs of the files =============================
sub GetFileStem($) {
  my $p = shift;
  if(!(-f "$p"))  {
    # does not seem to be a file
    die "FATAL ERROR: The input file \"$p\" does not exist\n";
  }
  if($p =~ /.*\/+(\S+)/)  {       # if the directory is not ended with back slash
    return $1;
  } else {                        # if no back slash exists, the stem is itself
    return $p;
  }
}

my $i; my $j; my $k;

# handling treatment groups
my $min_len = 999999;
my $trailing_index;
foreach(@treatments)  {
  my $f_stem = GetFileStem($_);
  $min_len = length($f_stem) if length($f_stem) < $min_len;
  push @treatment_IDs, $f_stem;
}
for($trailing_index = 1; $trailing_index < $min_len; ++ $trailing_index) {
  my $c = substr($treatment_IDs[0], -$trailing_index, 1);
  my $break_tag = 0;
  for($i = 1; $i < scalar(@treatment_IDs); ++ $i) {
    if(substr($treatment_IDs[$i], -$trailing_index, 1) ne $c) { $break_tag = 1; last; }
  }
  last if $break_tag;
}
for($i = 0; $i < scalar(@treatment_IDs); ++ $i) {
  $treatment_IDs[$i] = substr($treatment_IDs[$i], 0, length($treatment_IDs[$i]) - $trailing_index + 1);
}

# TODO: check for special characters in the IDs
foreach(@treatment_IDs) {
  if(/\-/)  {  
    die "FATAL ERROR: Special character \"-\" is found in the file name of $_; rename before rerunning.\n";
  }
}

# handling control groups
$min_len = 999999;
foreach(@controls)  {
  my $f_stem = GetFileStem($_);
  $min_len = length($f_stem) if length($f_stem) < $min_len;
  push @control_IDs, $f_stem;
}
for($trailing_index = 1; $trailing_index < $min_len; ++ $trailing_index) {
  my $c = substr($control_IDs[0], -$trailing_index, 1);
  my $break_tag = 0;
  for($i = 1; $i < scalar(@control_IDs); ++ $i) {
    if(substr($control_IDs[$i], -$trailing_index, 1) ne $c) { $break_tag = 1; last; }
  }
  last if $break_tag;
}
for($i = 0; $i < scalar(@control_IDs); ++ $i) {
  $control_IDs[$i] = substr($control_IDs[$i], 0, length($control_IDs[$i]) - $trailing_index + 1);
}

foreach(@control_IDs) {
  if(/\-/)  {  
    die "FATAL ERROR: Special character \"-\" is found in the file name of $_; rename before rerunning.\n";
  }
}


#==================== run FASTQC to check the data quality =====================
#==================== also generate overrepresented sequences ==================
print "!!!  TranDER: running FASTQC to check data quality...\n";

for($i = 0; $i < scalar(@treatments); ++ $i) {
  mkdir "$work_dir/FASTQC_Report/$treatment_IDs[$i]" 
    or die "FATAL ERROR: fail to create FASTQC output directory!" if !(-e "$work_dir/FASTQC_Report/$treatment_IDs[$i]");
#  system "$exe_fastqc -o $work_dir/FASTQC_Report/$treatment_IDs[$i] $treatments[$i] >&$work_dir/Log/fastqc.$treatment_IDs[$i].log";
  print "!!!  TranDER: quality assessment for $treatment_IDs[$i] done. See log file at $work_dir/Log/fastqc.$treatment_IDs[$i].log\n";
  foreach(<$work_dir/FASTQC_Report/$treatment_IDs[$i]/*.zip>) {
    my $fastqc_id = $_;
    $fastqc_id =~ s/\.zip//g;
#    system "unzip -o -d $work_dir/FASTQC_Report/$treatment_IDs[$i] $_ >&$work_dir/Log/unzip.$treatment_IDs[$i].log";
#    system "perl $HOME/scripts/GetOverrepresentedSeqs.pl $fastqc_id/fastqc_data.txt $work_dir/Intermediate/$treatment_IDs[$i].overrepresented.fasta";
  }
}
for($i = 0; $i < scalar(@controls); ++ $i) {
  mkdir "$work_dir/FASTQC_Report/$control_IDs[$i]" 
    or die "FATAL ERROR: fail to create FASTQC output directory!" if !(-e "$work_dir/FASTQC_Report/$control_IDs[$i]");
#  system "$exe_fastqc -o $work_dir/FASTQC_Report/$control_IDs[$i] $controls[$i] >&$work_dir/Log/fastqc.$control_IDs[$i].log";
  print "!!!  TranDER: quality assessment for $control_IDs[$i] done. See log file at $work_dir/Log/fastqc.$control_IDs[$i].log\n";
  foreach(<$work_dir/FASTQC_Report/$control_IDs[$i]/*.zip>) {
    my $fastqc_id = $_;
    $fastqc_id =~ s/\.zip//g;
#    system "unzip -o -d $work_dir/FASTQC_Report/$control_IDs[$i] $_ >&$work_dir/Log/unzip.$control_IDs[$i].log";
#    system "perl $HOME/scripts/GetOverrepresentedSeqs.pl $fastqc_id/fastqc_data.txt $work_dir/Intermediate/$control_IDs[$i].overrepresented.fasta";
  }
}


#==================== run Trimmomatic for each data set ========================
print "!!!  TranDER: running Trimmomatic to trim reads...\n";
for($i = 0; $i < scalar(@treatments); ++ $i) {
#  system "java -jar $exe_trimmomatic SE $treatments[$i] $work_dir/Trimmed_Reads/$treatment_IDs[$i].passed.fq.gz ILLUMINACLIP:$work_dir/Intermediate/$treatment_IDs[$i].overrepresented.fasta:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:50 >&$work_dir/Log/trimmomatic.$treatment_IDs[$i].log";
  print "!!!  TranDER: quality trimming for $treatment_IDs[$i] done. See log file at $work_dir/Log/fastqc.$treatment_IDs[$i].log\n";
}
for($i = 0; $i < scalar(@controls); ++ $i) {
#  system "java -jar $exe_trimmomatic SE $controls[$i] $work_dir/Trimmed_Reads/$control_IDs[$i].passed.fq.gz ILLUMINACLIP:$work_dir/Intermediate/$control_IDs[$i].overrepresented.fasta:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:50 >&$work_dir/Log/trimmomatic.$control_IDs[$i].log";
  print "!!!  TranDER: quality trimming for $control_IDs[$i] done. See log file at $work_dir/Log/fastqc.$control_IDs[$i].log\n";
}


#==================== run HISAT2 mapper ========================================
# check if all indexes exist for HISAT2
for($i = 1; $i <= 8; ++ $i) {
  if(!-e "$ref_prefix\.$i\.ht2" && !-e "$ref_prefix\.$i\.ht2l") {
    print "FATAL ERROR: No HISAT2 index is found from the reference index directory \"$ref_prefix\"!\n";
    print "Please specify \"Reference_prefix\" in the manifest.\n";
    die "To build the HISAT2 reference,\n\tuse \"hisat2-build <ref_genome.fasta> <$ref_prefix>\"\n";
  }
}

# start running HISAT2
print "!!!  TranDER: running HISAT2 to map the reads to the reference...\n";
for($i = 0; $i < scalar(@treatments); ++ $i) {
#  system "$exe_hisat2 -x $ref_prefix -U $work_dir/Trimmed_Reads/$treatment_IDs[$i].passed.fq.gz -p $num_threads -S $work_dir/Mapping/$treatment_IDs[$i].sam >&$work_dir/Log/hisat2.$treatment_IDs[$i].log";
  print "!!!  TranDER: read mapping for $treatment_IDs[$i] done. See log file at $work_dir/Log/fastqc.$treatment_IDs[$i].log\n";
}
for($i = 0; $i < scalar(@controls); ++ $i) {
#  system "$exe_hisat2 -x $ref_prefix -U $work_dir/Trimmed_Reads/$control_IDs[$i].passed.fq.gz -p $num_threads -S $work_dir/Mapping/$control_IDs[$i].sam >&$work_dir/Log/hiset2.$control_IDs[$i].log";
  print "!!!  TranDER: read mapping for $control_IDs[$i] done. See log file at $work_dir/Log/fastqc.$control_IDs[$i].log\n";
}


#==================== run eXpress to estimate the read count ===================
# TODO: esitmate the forward read count and reverse read count to run eXpress more accurately

print "!!!  TranDER: running eXpress to estimate the read count for each transcript..\n";
for($i = 0; $i < scalar(@treatments); ++ $i) {
#  system "$exe_express -o $work_dir/Transcript_abundance/$treatment_IDs[$i] --r-strand $ref_transcript $work_dir/Mapping/$treatment_IDs[$i].sam >&$work_dir/Log/express.$treatment_IDs[$i].log";
  print "!!!  TranDER: count estmation for $treatment_IDs[$i] done. See log file at $work_dir/Log/express.$treatment_IDs[$i].log\n";
}
for($i = 0; $i < scalar(@controls); ++ $i) {
#  system "$exe_express -o $work_dir/Transcript_abundance/$control_IDs[$i] --r-strand $ref_transcript $work_dir/Mapping/$control_IDs[$i].sam >&$work_dir/Log/express.$control_IDs[$i].log";
  print "!!!  TranDER: count estmation for $control_IDs[$i] done. See log file at $work_dir/Log/express.$control_IDs[$i].log\n";
}


#==================== summarize all count into a table =========================
print "!!!  TranDER: preparing files for differential expression analysis...\n";
# generate the count table
system "touch $work_dir/Transcript_abundance/file_list";
open my $FOUT, ">$work_dir/Transcript_abundance/file_list" or die "FATAL ERROR: Cannot create file to write abundance estimation list.\n";
for($i = 0; $i < scalar(@treatments); ++ $i) {
#  system "cp $work_dir/Transcript_abundance/$treatment_IDs[$i]/results.xprs $work_dir/Transcript_abundance/$treatment_IDs[$i]_count";
  print $FOUT "$work_dir/Transcript_abundance/$treatment_IDs[$i]_count\n";
}
for($i = 0; $i < scalar(@controls); ++ $i) {
#  system "cp $work_dir/Transcript_abundance/$control_IDs[$i]/results.xprs $work_dir/Transcript_abundance/$control_IDs[$i]_count";
  print $FOUT "$work_dir/Transcript_abundance/$control_IDs[$i]_count\n";
}
close $FOUT;
system "perl $HOME/scripts/SummarizeCountTable.pl $work_dir/Transcript_abundance/file_list >$work_dir/Differential_expression/count_table.tab";

# generate the condition table
open my $COND_OUT, ">$work_dir/Differential_expression/condition.tab" or die "FATAL ERROR: Cannot create condition table for DESeq2.\n";
print $COND_OUT "NAME\tCondition\tType\n";
for($i = 0; $i < scalar(@treatment_IDs); ++ $i) {
  print $COND_OUT "$treatment_IDs[$i]\_count\ttreatment\ttreatment_rep$i\n";
}
for($i = 0; $i < scalar(@control_IDs); ++ $i) {
  print $COND_OUT "$control_IDs[$i]\_count\tcontrol\tcontrol_rep$i\n";
}
close $COND_OUT;

#==================== run DESeq2 to compute DE genes ===========================
print "!!!  TranDER: running DESeq2 to identify differentially expressed genes...\n";
# first check whether DESeq2 is installed
system "$exe_Rscript $HOME/scripts/CheckDESeq2.R >&$work_dir/Intermediate/check_DESeq2";
open my $CK_IN, "<$work_dir/Intermediate/check_DESeq2" or die "FATAL ERROR: unable to execute R script or write to folder \"$work_dir/Intermediate\"\n";
my $ck_DESeq2 = <$CK_IN>;
die "FATAL ERROR: DESeq2 may not be installed; try install DESeq2 package in R.\n" if $ck_DESeq2 =~ /^Error/;
close $CK_IN;
# run DESeq2

system "$exe_Rscript $HOME/scripts/RunDESeq2.R $work_dir/Differential_expression/count_table.tab $work_dir/Differential_expression/condition.tab $identifier $work_dir/Differential_expression >&$work_dir/Log/DESeq2.log";
print "!!!  TranDER: DESeq2 run done. See log file at $work_dir/Log/DESeq2.log\n";
system "perl $HOME/scripts/ParseDESeq2Output.pl $work_dir/Differential_expression/diffexpr.$identifier.tab $work_dir/Report/diffexpr.$identifier.parsed.tab";
system "cp $work_dir/Differential_expression/MAplot.$identifier.pdf $work_dir/Report";

#==================== end of TranDER pipeline ==================================



