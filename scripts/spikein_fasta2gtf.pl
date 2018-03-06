#!/usr/bin/env perl
# =========================================================
# Copyright 2012-2018,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This file is part of iRAP.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================
use strict;

use Bio::SeqIO;

my $input_file1 = shift;
## usage: spikein_fasta2gtf.pl input_fasta_file
die("Unable to find $input_file1") if ! -e $input_file1;
##die("$input_file1 should be a file") if ! -f $input_file1;
die("Unable to open $input_file1") if ! -r $input_file1;


my $seq_in  = Bio::SeqIO->new( -format => 'fasta',
                               -file => $input_file1);

while (my $inseq = $seq_in->next_seq) {
    my $len=length($inseq->seq);
    my $id=$inseq->id;
    #print "$id\tprotein_coding\texon\t1\t$len\t.\t+\t.\tgene_id \"$id\"; transcript_id \"t$id\"; exon_number \"1\"; gene_bio_type \"spikein\" \n";
    print "$id\tspikein\texon\t1\t$len\t.\t+\t.\tgene_id \"$id\"; transcript_id \"$id\"; exon_number \"1\"; gene_name \"$id\"; exon_id \"e${id}_1\"; gene_biotype \"protein_coding\"; gene_type \"protein_coding\";\n";
    #print "$id\tspikein\tCDS\t1\t$len\t.\t+\t.\tgene_id \"$id\"; transcript_id \"t$id\"; gene_bio_type \"protein_coding\"; transcript_name \"t$id\";\n";
    print "$id\tspikein\ttranscript\t1\t$len\t.\t+\t.\tgene_id \"$id\"; transcript_id \"$id\"; gene_biotype \"protein_coding\"; transcript_biotype \"protein_coding\"; transcript_type \"protein_coding\"; transcript_name \"t$id\";\n";
    print "$id\tspikein\tgene\t1\t$len\t.\t+\t.\tgene_id \"$id\"; gene_biotype \"protein_coding\"; gene_type \"protein_coding\"; gene_name \"$id\";\n";
}

exit(0);
