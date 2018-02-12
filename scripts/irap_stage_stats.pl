#!/usr/bin/env perl
# Author:      Nuno A. Fonseca <nuno-dot-fonseca-at-gmail-dot-pt>
use strict;

if ( $ARGV[6] eq "" ) { 
    die("ERROR: Usage: irap_stage_stats stage name qc mapper quant de file1 file2 ...");
}
my $stage=shift @ARGV;
my $name=shift @ARGV;
my $qc=shift @ARGV;
my $mapper=shift @ARGV;
my $quant=shift @ARGV;
my $de=shift @ARGV;
my @files=@ARGV;

my $nfiles=scalar(@files);
my $nfiles_ok=0;


sub get_ext {
    my $file=shift;    
    $file=~/.([a-z-A-Z0-9]+)$/;
    return $1
}
sub get_file_type {
    
    my @files=@_;
    my $ext="";

    for my $f (@files) {
	my $new_ext=get_ext($f);
	if ($ext ne "" && $ext ne $new_ext ) {
	    die("get_file_type:: $f file extension inconsistent (expected $ext and got $new_ext)");
	}
	$ext=$new_ext;
	print STDERR "$f $new_ext\n";
    }
    return $ext;
}

sub valid_file {
    my $f=shift;
    my $t=shift;

    return(0) if ( ! -e $f );
    if ( $t eq "bam" ) {
	my $l=`samtools view -c $f` or return 0;
	chomp $l;
	return(0) if ( $l == 0 );
	print STDERR "[INFO] $f - $l\n";
    }
    if ( $t eq "fastq" ) {
	my $l=`wc -l $f | cut -f 1 -d\\ ` or return 0;
	chomp $l;
	return(0) if ( $l == 0 );
	print STDERR "[INFO] $f - $l\n";
    }
    if ( $t eq "tsv" ) {
	my $l=`wc -l $f | cut -f 1 -d\\  ` or return 0;
	chomp $l;
	return(0) if ( $l == 0 );
	# TODO: check the number of columns
	my $h=`head -n 1 $f` or return 0;
	chomp $h;
	my $n=scalar(split(/\t/,$h));
	return(0) if ( $n < 1 );
	    
	print STDERR "[INFO] $f - $l  $n\n";
    }
    return(1);
}
#my $ext=get_file_type @files;
my $ext="";
for my $f (@files) {
    if ( valid_file($f,$ext) ) {
	$nfiles_ok++;
    } else {
	print STDERR "[WARNING] file $f not ok\n";
    }
}
my $p_done=int($nfiles_ok/$nfiles*100);
print "$name,$stage,$qc,$mapper,$quant,$de,$nfiles,$nfiles_ok,$p_done\n";
exit;
