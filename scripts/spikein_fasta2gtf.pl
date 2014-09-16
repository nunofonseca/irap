#!/usr/bin/env perl
# Author:      Nuno A. Fonseca <nf-at-ibmc-dot-up-dot-pt>

use strict;

use Bio::SeqIO;

my $input_file1 = shift;

die("Unable to find $input_file1") if ! -e $input_file1;
die("$input_file1 should be a file") if ! -f $input_file1;
die("Unable to open $input_file1") if ! -r $input_file1;


my $seq_in  = Bio::SeqIO->new( -format => 'fasta',
                               -file => $input_file1);

while (my $inseq = $seq_in->next_seq) {
    my $len=length($inseq->seq);
    my $id=$inseq->id;
    #print "$id\tprotein_coding\texon\t1\t$len\t.\t+\t.\tgene_id \"$id\"; transcript_id \"t$id\"; exon_number \"1\"; gene_biotype \"spikein\" \n";
    print "$id\tspikein\texon\t1\t$len\t.\t+\t.\tgene_id \"$id\"; transcript_id \"t$id\"; exon_number \"1\"; gene_biotype \"spikein\" \n";
    print "$id\tCDS\texon\t1\t$len\t.\t+\t.\tgene_id \"$id\"; transcript_id \"t$id\"; exon_number \"1\"; gene_biotype \"spikein\"; transcript_name \"t$id\"; \n";
}

exit(0);
