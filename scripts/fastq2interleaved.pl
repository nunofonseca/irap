#!/usr/bin/env perl
# Author:      Nuno A. Fonseca <nf-at-ibmc-dot-up-dot-pt>

use strict;

use Bio::SeqIO;

my $input_file1 = shift;
my $input_file2 = shift;
my $output_file = shift;

my $seq_in1  = Bio::SeqIO->new( -format => 'fastq',
                               -file => $input_file1);
my $seq_in2  = Bio::SeqIO->new( -format => 'fastq',
                               -file => $input_file2);

my $out1 = Bio::SeqIO->new('-file' => ">$output_file",
                                       '-format' =>"fastq");

while (my $inseq1 = $seq_in1->next_seq) {
    my $inseq2 = $seq_in2->next_seq;
    $out1->write_seq($inseq1);
    $out1->write_seq($inseq2);
}
exit(0);
