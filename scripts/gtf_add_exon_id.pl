#!/usr/bin/env perl
# POD documentation - main docs before the code
=pod
=head1 NAME
  gtf_add_exon_id.pl
  Nuno A. Fonseca (first name.last name@gmail.com) - 21-12-2012

=head1 SYNOPSIS
 Adds unique exon_id field to a given GTF file.
  
=head1 DESCRIPTION
 Adds unique exon_id field base on gene name and start_pos to a given GTF file.

=head1 OPTIONS
  None

=head1 EXAMPLES
  gtf_add_exon_id.pl gtf

=cut
# $Id: 0.1.3 Nuno Fonseca Fri Dec 21 16:07:17 2012$
use strict;
use warnings;

die "ERROR\nUsage: gtf_add_exon_id.pl gtf/bed_file \n" if ( $#ARGV!=0 );

my $debug=0;
#my $countby = shift; # gene|exon
my $gtf_file;

$gtf_file=$ARGV[0];


sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}
# 
open (FILE, $gtf_file) or die($!);
print STDERR "Loading ".$gtf_file."\n";
while (<FILE>) {
    chomp;
    if ( !  /^#/ ) {
    # TODO: handle comments
	my ($seqname,$source,$feature,$start,$end,$score,$strand,$frame,$attr) = split("\t");
	my $count=0;
	my %vars=("gene_id"=>"","transcript_id"=>"","exon_number"=>"","gene_name"=>"","gene_bio_type"=>"");
	for my $attv (split(/;/,$attr)) {
	    my ($name,$val)=split(/ /,trim($attv));
	    $vars{$name}=$val;	
	    $vars{$name}=~s/"//g;
	}
	#$vars{"exon_id"}=$vars{"gene_id"}.".".$vars{"transcript_id"}.".".$vars{"exon_number"};
	$vars{"exon_id"}=$vars{"gene_id"}.".".$vars{"exon_number"};
	$attr=$attr." exon_id \"".$vars{"exon_id"}."\";";
	print join("\t",($seqname,$source,$feature,$start,$end,$score,$strand,$frame,$attr))."\n";
    } else {
	print $_."\n";
    }
}
close(FILE);
