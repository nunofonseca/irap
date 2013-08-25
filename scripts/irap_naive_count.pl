#!/usr/bin/env perl
# POD documentation - main docs before the code
=pod
=head1 NAME
  irap_naive_count.pl
  Nuno A. Fonseca (first name.last name@gmail.com) - 07-10-2010

=head1 SYNOPSIS
  Compute (naively) the number of reads per gene or exon from a BED or BAM file.
  
=head1 DESCRIPTION
  Compute (naively) the number of reads per gene or exon from a BAM file and a BED/GTF file.

=head1 OPTIONS
  None

=head1 EXAMPLES
  irap_naive_count.pl bed_file bam_file exons.tsv genes.tsv
  irap_naive_count.pl bed_file exons.tsv genes.tsv

=cut

use List::MoreUtils 'uniq';

use strict;
use warnings;

die "ERROR\nUsage: irap_naive_count.pl gtf/bed_file [bam_file] countsbyexon_file countsbygene_file\n" if ( $#ARGV!=3 && $#ARGV!=2 );

my $debug=0;
#my $countby = shift; # gene|exon
my $gtf_file;
my $bam_file;
my $exons_counts_file;
my $genes_counts_file;

($gtf_file,$bam_file,$exons_counts_file,$genes_counts_file)=($ARGV[0],$ARGV[1],$ARGV[2],$ARGV[3]) if ( $#ARGV==3 );
($gtf_file,$bam_file,$exons_counts_file,$genes_counts_file)=($ARGV[0],$ARGV[0],$ARGV[1],$ARGV[2]) if ( $#ARGV==2 );

print STDERR "irap_naive_count.pl $gtf_file $bam_file $exons_counts_file $genes_counts_file\n" if ($#ARGV==3 );
print STDERR "irap_naive_count.pl $gtf_file  $exons_counts_file $genes_counts_file\n" if ($#ARGV==2 );
my $tmp_file = "$gtf_file.tmp";
# my $tmp_file = "filt.bed.cov";
# my $tmp_file = "chr19.cov.bed";
# $countby="exon";
my $countby="exon";

for my $f ($gtf_file,$bam_file) {
    die "File $f not found\n" if ( ! -e $f ) 
}

# Output format
# geneid(exonloc) count len

my $output_dir = shift;

open(EXONS_FILE, ">$exons_counts_file") or die("Can't open $exons_counts_file");
open(GENES_FILE, ">$genes_counts_file") or die("Can't open $genes_counts_file");

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

my %genes;
my %exons;
my %coords;
my %trans;

if (  $gtf_file ne $bam_file ) {
    # TODO: create a function!!!
    open (FILE, $gtf_file) or die($!);
    print STDERR "Loading ".$gtf_file."\n";
    while (<FILE>) {
	chomp;
	# TODO: handle comments
	my ($seqname,$source,$feature,$start,$end,$score,$strand,$frame,$attr) = split("\t");
	my $count=0;
	my %vars=("gene_id"=>"","transcript_id"=>"","exon_number"=>"","gene_name"=>"","gene_bio_type"=>"");
	for my $attv (split(/;/,$attr)) {
	    my ($name,$val)=split(/ /,trim($attv));
	    $vars{$name}=$val;	
	    $vars{$name}=~s/"//g;
	}
	if ( $feature eq "exon" ) {
	    my $key=$vars{"gene_id"}."-".$vars{"transcript_id"}."-".$vars{"exon_number"};
	    print STDERR "WARNING: repeated exon for ".$vars{"gene_id"}.".Using the largest value.\n" if defined $exons{$key};
	    $exons{$key}=0;
	    $genes{$vars{"gene_id"}}=0;
	    $trans{$vars{"gene_id"}}=0;
	    $coords{$key}="$start $end";
	}
    }
    close(FILE);
    print STDERR "Generating BED file $tmp_file...this may take a while\n";
    `bedtools coverage -abam $bam_file -b $gtf_file -split > $tmp_file`;

} else { 
    $tmp_file=$gtf_file
}

#<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
open (FILE, "$tmp_file");
# features: "CDS", "start_codon", "stop_codon". The features "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS" and "exon" are optional. 
while (<FILE>) {
    chomp;
    # TODO: handle comments
    my ($seqname,$source,$feature,$start,$end,$score,$strand,$frame,$attr,$count) = split("\t");
    #print STDERR "$seqname, $source, $feature, $start, $end, $score, $strand,$frame, $attr, $count\n";
    # get the gene id, transcript, exon, gene_name and gene_biotype from the attributes
    my %vars=("gene_id"=>"","transcript_id"=>"","exon_number"=>"","gene_name"=>"","gene_bio_type"=>"");
    for my $attv (split(/;/,$attr)) {
	my ($name,$val)=split(/ /,trim($attv));
	$vars{$name}=$val;	
	$vars{$name}=~s/"//g;
    }
    if ( $feature eq "exon" ) {
	my $key=$vars{"gene_id"}."-".$vars{"transcript_id"}."-".$vars{"exon_number"};
	print STDERR "WARNING: repeated exon for ".$vars{"gene_id"}.".Using the largest value.\n" if defined $exons{$key};
	if ( defined($count) && defined($exons{$key}) ) {
	    if ($count > $exons{$key} ) {
		$exons{$key}=$count;
		$genes{$vars{"gene_id"}}="" if ! defined $genes{$vars{"gene_id"}};
		$trans{$vars{"gene_id"}}="" if ! defined $trans{$vars{"gene_id"}};
		$coords{$key}="$start $end";
		$genes{$vars{"gene_id"}}=$genes{$vars{"gene_id"}}." ".$vars{"exon_number"};
		$trans{$vars{"gene_id"}}=$trans{$vars{"gene_id"}}." ".$vars{"transcript_id"};
	    }
	}
    }
}
close(FILE);

print STDERR "Found ".(keys %exons)." exons\n";
print STDERR "Found ".(keys %genes)." genes\n";
for my $gene (keys %genes) {
    print STDERR "Processing $gene\n" if ($debug);
    my @transcripts=uniq(split(/ /,trim($trans{$gene})));
    print STDERR "#transcripts:".@transcripts."\n" if ($debug);
    for my $t (@transcripts) {
	my @transcript_exons=grep(/$gene-$t-.*/,keys(%exons));
	print STDERR "Processing transcript $t: ".@transcript_exons." exons\n" if ($debug);	
    }
    # get all exons    
    my @exons=sort(grep(/$gene-.*/,keys(%exons)));
    ############################
    # exclude 'duplicated' exons
    my %exons2consider=();
    for my $e (@exons) {
	#print STDERR "$e\n";
	my $discard=0;
	my $e_coords=$coords{$e};
	for my $e2 (keys %exons2consider) {
	    # check if it is the same exon
	    my $e2_coords=$coords{$e2};
	    if ( $e2_coords eq $e_coords ) {
		print STDERR "Discarding $e (==$e2)\n" if ($debug);
		$discard=1;
		last;
	    } else {
		# check if there is an overlap
		my ($start1,$end1)=split(/ /,$e_coords);
		my ($start2,$end2)=split(/ /,$e2_coords);
		if ( $start1 > $end1 ) { print STDERR "COORDS\n"; exit 1;}
		if ( $start2 > $end2 ) { print STDERR "COORDS\n"; exit 1;}
		if ( $end2<$start1 || $end1<$start2) {
		    #print STDERR "No overlap\n";
		} else {
		    my $len1=$end1-$start1;
		    my $len2=$end2-$start2;
		    print STDERR "Exons overlap: $e_coords ($len1) $e2_coords ($len2)\n" if ($debug);
		    # pick the larger exon
		    if ( $len1 > $len2) {
			$discard=0;
			# remove 
			delete $exons2consider{$e2};
			print STDERR "Exons overlap: removing $e2 \n" if ($debug);
			last;
		    } else {
			$discard=1;
		    }
		}
	    }
	}	    
	if ( $discard == 0 ) {
	    print STDERR "Keeping $e \n" if ($debug);
	    $exons2consider{$e}=$exons{$e}
	} else {
	    print STDERR "Exons overlap: discarding $e \n" if ($debug);
	}
    }
    print STDERR "Exons considered ".keys(%exons2consider)."\n" if ($debug);
    #if ($countby eq "gene") {
	my $tot_count=0;
	my $tot_len=0;
	for my $ec (keys %exons2consider) {
	    $tot_count+=$exons2consider{$ec};
	    my ($start,$end)=split(/ /,$coords{$ec});
	    $tot_len+=abs($end-$start);
	}
	print GENES_FILE "$gene\t$tot_count\t$tot_len\n";
    #} else {
	for my $ec (keys %exons2consider) {
	    my $count=$exons2consider{$ec};
	    my ($start,$end)=split(/ /,$coords{$ec});
	    my $len+=abs($end-$start);
	    print EXONS_FILE "$gene.$start.$end\t$count\t$len\n";
	}		
    #}
    
}
if ( $#ARGV==3 ) {
    `rm -f $tmp_file`;
}
close(EXONS_FILE);
close(GENES_FILE);

exit;

# #read from a file
# my $in  = Bio::FeatureIO->new(-file => "$tmp_file" , -format => 'GTF');
# my $tag;
# while ( my $feat = $in->next_feature() ) {
#     my $chr=$feat->get_tag_values("seq_id");
#     my $strand=$feat->strand;
#     my $tag=$feat->primary_tag;
#     my $start=$feat->start;
#     my $end=$feat->end;
#     print STDERR "CHR=$chr strand=$strand tag=$tag   $start->$end\n";
#     my $str = $feat->gff_string;
#     print STDERR "$str\n";
#     print STDERR "Feature from ", , "to ",
#                $feat->end, " Primary tag  ", $feat->primary_tag,
#                   ", produced by ", $feat->source_tag(), "\n";

#        print STDERR "feature location is ",$feat->start, "..",
#           $feat->end, " on strand ", $feat->strand, "\n";
#        print STDERR "easy utility to print STDERR locations in GenBank/EMBL way ",
#           $feat->location->to_FTstring(), "\n";

#        foreach $tag ( $feat->get_all_tags() ) {
#                     print STDERR "Feature has tag ", $tag, " with values, ",
#                       join(' ',$feat->get_tag_values($tag)), "\n";
#        }
#        print STDERR "new feature\n" if $feat->has_tag('new');
#        # features can have sub features
#        my @subfeat = $feat->get_SeqFeatures();
# }
#use Bio::FeatureIO;
