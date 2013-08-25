#!/usr/bin/env perl
# Author:      Nuno A. Fonseca <nuno-dot-fonseca-at-gmail-dot-pt>
use Bio::SeqIO;
use strict;
my $nfiles=1;

my $input_file = shift;
my $output_dir = shift;


my $seq_in  = Bio::SeqIO->new( -format => 'fasta',
			       -file => $input_file);
`mkdir -p $output_dir`;
`rm -f $output_dir.lst; touch $output_dir.lst`;
my $output_file_prefix="$output_dir/";

# write each entry in the input file to the output file
while (my $inseq = $seq_in->next_seq) {
    print "Creating file $nfiles\n";
    my $ofile=">$output_file_prefix".$inseq->id().".fa";
    my $seq_out = Bio::SeqIO->new('-file' => $ofile,
			       '-format' =>"fasta");
    ++$nfiles;	
    `echo $ofile >> $output_dir.lst`;
    $seq_out->write_seq($inseq);
    $seq_out->close;
}
