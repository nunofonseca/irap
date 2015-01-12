#!/usr/bin/env perl
# convert gtf file from Ensembl to gff3 file usable by tophat
# tested on ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh37.56.gtf.gz

use strict;


my @v;
my ($g,$tr,$last_tr,$last_g);
my @exons=();
my @trs=();
my $counter = 0;
while (<>) {
        next if ( /^#/ ); # ignore comments
	@v = split /\t/;
	next unless $v[2] eq 'exon';
	#$v[0] = 'chr' . $v[0] unless index($v[0], 'chr') == 0;
	($g,$tr) = $v[8] =~ /gene_id \"?([^;\"]+)\"?;.* transcript_id \"?([^;\"]+)\"?/;
	if ($tr ne $last_tr) {
		push @trs, [@exons];
		if ($g ne $last_g) {
			process(@trs) if defined $last_g;
			@trs=();
			$last_g = $g;
			print STDERR '.' if $counter % 1000 == 0;
		}
		$last_tr = $tr;
		@exons = [@v];
	} else {
		push @exons, [@v];
	}
}
push @trs, [@exons];
process(@trs);

sub process {
	my @trs = @_;
	return unless scalar @trs;
	my @all_starts = ();
	my @all_ends = ();
	my @out = ();
	my @tr_here;
	my ($gid, $tid, $exn, $gname);
	my (@starts, @ends, $start, $end);
	my (@transcript, @gene);
	my ($exon_info_stringA, $exon_info_stringB);


	($gid) = $trs[0]->[0]->[8] =~ /gene_id \"?([^;\"]+)\"?/;
	($gname) = $trs[0]->[0]->[8] =~ /gene_name \"?([^;\"]+)\"?/;
	$gname =~ tr/;/./;

	foreach my $exons (@trs) {
		@starts = sort {$a<=>$b} map { $_->[3] } @$exons;
		@ends = sort {$b<=>$a} map { $_->[4] } @$exons;
		$start = $starts[0];
		$end = $ends[0];

		($tid) = $exons->[0]->[8] =~ /transcript_id \"?([^;\"]+)\"?/;
		$exon_info_stringA = "ID=$tid.";
		$exon_info_stringB = ";Name=$gname;Parent=$tid\n";

		@transcript = @{$exons->[0]};
		$transcript[2] = 'mRNA';
		$transcript[3] = $start;
		$transcript[4] = $end;
		$transcript[8] =  "ID=$tid;Name=$gname;Parent=$gid\n";

		push @all_starts, $start;
		push @all_ends, $end;
		foreach (reverse @$exons) {
			($exn) = $_->[8] =~ /exon_number "(\d+)"/;
			$_->[8] = $exon_info_stringA . $exn. $exon_info_stringB;
			push @out, join("\t", @$_);
		}
		push @out, join("\t", @transcript);
	}
	@all_starts = sort {$a<=>$b} @all_starts;
	@all_ends = sort {$b<=>$a} @all_ends;
	@gene = @transcript;
	$gene[2] = 'gene';
	$gene[3] = $all_starts[0];
	$gene[4] = $all_ends[0];
	$gene[8] = "ID=$gid;Name=$gname\n";
	push @out, join("\t", @gene);
	print reverse @out;
}
