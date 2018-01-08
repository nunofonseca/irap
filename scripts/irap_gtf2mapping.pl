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
use Getopt::Long;
use strict;

sub byebye {
    my ($msg) = @_;
    print STDERR $msg;
    exit(1);
}

my $usage="irap_gtf_to_mapping --gtf gtf_file --feature exon|transcript";

my $gtf_file = "";
my $feature = "";
# types of transcripts to keep (if more than one then separate them by |
my $out_file = "";
GetOptions ("gtf|g=s"   => \$gtf_file,      # string
	    "feature|f=s"   => \$feature)      # string	    
  or byebye("Error:".$usage."\n");

# TODO: check values
byebye("Error:".$usage."\n") if ($gtf_file eq "" || $feature eq "" );

byebye("Error:".$gtf_file." not found.\n") if ! -e $gtf_file;

print STDERR "gtf=$gtf_file\n";
print STDERR "feature=$feature\n";
## check if gtf_to_fasta is in the path

if ( $feature!="exon" || $feature!="transcript" ) {
    print STDERR "Warning: expected exon or transcript as feature (got $feature)"
}
print STDERR "Processing GTF...\n";

my $fid=$feature."_id";
print STDERR $fid."\n";
my $feat_pat=$feature;
my $fid2="";
my $feat_trans_id="transcript_id";
if ( $feature eq "transcript" ) {
    $feat_pat="(CDS|transcript)";
} else {
    $fid2="exon_number" if ( $feature eq "exon" );
}

## load the mapping between the transcripts and the strand
open(GTF,$gtf_file) || byebye "Could not open $gtf_file";


my $linenum=0;
print "ID\tgene_id\tseqid\tstart\tend\tstrand\tsource\n";
while(<GTF>) {
    chomp($_);
    ++$linenum;    
    my @line = split(/\t/,$_);
    $_=$line[2];
    next if ( ! m/$feat_pat/ );
    ##print STDERR "---$line[8]---\n";
    my $s1=$line[8];
    ## gene_id
    my $gid=$s1;
    $gid=~ s/^.*\s*gene_id\s*\"?//g;
    $gid=~ s/"?;.*//g;
    ## feature
    my $id=$line[8];
    $_=$id;
    if ( ! m/.*\s*$fid\s*\"/ ) {
	byebye "no $fid not found in line $linenum\n" if ( $fid2 eq "");
	byebye "no $fid or $fid2 not found in line $linenum\n" if (! m/.*\s*$fid2\s*\"/ );
	$id=~ s/.*\s*$fid2\s*\"?//g;
	$id=~ s/"?\s*;.*//g;
	## transcript id
	my $trans_id=$line[8];
	$_=$trans_id;	
	byebye "no $feat_trans_id  found in line $linenum\n" if (! m/.*\s*$feat_trans_id\s*\"/ );
	$trans_id=~ s/.*\s*$feat_trans_id\s*\"?//g;
	$trans_id=~ s/"?\s*;.*//g;
	$id=$gid.".".$trans_id.".".$id;
    } else {
	$id=~ s/.*\s*$fid\s*\"?//g;
	$id=~ s/"?\s*;.*//g;
    }
    print "$id\t$gid\t".$line[0]."\t".$line[3]."\t".$line[4]."\t".$line[6]."\t".$line[1]."\n";
}
exit(0);
