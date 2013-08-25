#!/usr/bin/env perl
# Author:      Nuno A. Fonseca <nuno-dot-fonseca-at-gmail-dot-pt>
use Bio::SeqIO;
use strict;
my $input_file = shift;

my $seq_in  = Bio::SeqIO->new( -format => 'fasta',
			       -file => $input_file);
my $n=0;

# write each entry in the input file to the output file
while (my $inseq = $seq_in->next_seq) {
    my $l=length($inseq->seq);
    print $inseq->id()."\t".$l."\n";
}
exit 0;
