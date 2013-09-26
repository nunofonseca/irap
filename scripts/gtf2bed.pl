#!/usr/bin/env perl
# FILE: gtf2bed
# AUTHOR: William Stafford Noble
# CREATE DATE: 10 Jan 2004
# PROJECT: HS
use strict;

my $usage = "gtf2bed [-id <string>] <file>\n";

# Parse the command line.
my $id_append;
while (scalar(@ARGV) > 1) {
  my $next_arg = shift(@ARGV);
  if ($next_arg eq "-id") {
    $id_append = shift(@ARGV);
  } else {
    die("Invalid option ($next_arg).\n");
  }
}
if (scalar(@ARGV) != 1) {
  print(STDERR $usage);
  exit(1);
}
my($gtf_file) = @ARGV;

# Open the file for reading.
open(GTF_FILE, "<$gtf_file") || die("Can't open $gtf_file.");

# Read the file line by line.
while (my $line = <GTF_FILE>) {
  chomp($line);

  # Parse the line.
  my @words = split("\t", $line);
  my $chr = $words[0];
  my $start = $words[3];
  my $end = $words[4];
  (undef, my $id) = split(' ', $words[8]);
  $id = substr($id, 0, length($id) - 1);

  # Append to the ID, if requested.
  if (defined($id_append)) {
    $id .= "$id_append";
  }
  $id=~s/\"//g;
  # Print the corresponding BED line.
  printf("%s\t%d\t%d\t%s\n", $chr, $start - 1, $end, $id);
}
close(GTF_FILE);
