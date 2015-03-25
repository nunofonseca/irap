#!/bin/sh

file1=$1
threads=$2
paired=$3
outfile=$4
file2=$5


if [ "$outfile-" = "-" ]; then
    echo "Usage: irap_gen_fusionmap_conf.sh file1 threads paired outfile [file2]" >&2
    exit 1
fi

PAIRED=False
if [ "$paired" == "pe" ] ; then
    PAIRED=True
fi
# always false
PAIRED=FALSE
outdir=`dirname $outfile`
outfile=`basename $outfile`

cat <<EOF
<Files>
$file1
$file2

<Options>
ThreadNumber=$threads
PairedEnd=$PAIRED
RnaMode=True			//Detect fusion results 
FileFormat=BAM		//Possible values: FASTQ, QSEQ, FASTA, BAM. Default value=FASTQ
OutputFusionReads=True		//Out put Fusion reads as BAM files for genome browser. Default value=False
QualityEncoding=Automatic	//Auto detect quality coding in the fastq file or specify with Illumina or Sanger
AutoPenalty=True		//Set alignment penalty cutoff to automatic based on read length: Max (2,(read length-31)/15)
FixedPenalty=2			//If AutoPenalty=False, Fixed Penalty will be used
FilterUnlikelyFusionReads=True	//Enable filtering step
FullLengthPenaltyProportion=8	//Filtering normal reads allowing 8% of alignment mismatches of the reads
MinimalFusionAlignmentLength=0	//Default (alpha in the paper) value=0 and the program will automatically set minimal Seed Read end length to Min(20, Max(17,floor(ReadLength/3))). The program will use the specified value if user sets any > 0.
FusionReportCutoff=1		//# of allowed multiple hits of read ends; Possible values: 1-5. Default value=1 (beta in paper); 
NonCanonicalSpliceJunctionPenalty=4 //Possible values: 0-10. Default value= 2 (G); 
MinimalHit=4			//Minimal distinct fusion read; Possible values: 1-10000, Default value=2 
MinimalRescuedReadNumber=1	//Minimal rescued read number. Default value= 1
MinimalFusionSpan=5000		//Minimal distance (bp) between two fusion breakpoints
RealignToGenome=FALSE		//If True, seed read ends are re-aligned to genome to see if it is <= FusionReportCutoff in RNA-Seq.
OutputFusionReads=True		//Out put Fusion reads as BAM files for genome browser. Default value=True

<Output>
TempPath=$TEMPDIR
OutputPath=$outdir
OutputName=$outfile
EOF
exit 0
