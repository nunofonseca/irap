#!/bin/bash
set -e

if [ "$1-" = "-" ]; then
    echo "Usage: irapBAM2stats2.sh bam_file_name" > /dev/stderr
    exit 1
fi

if [ ! -e $1.ridx ]; then
    echo "Index not found: $1.ridx" > /dev/stderr
    exit 1
fi

time sqlite3 $1.ridx <<EOF
.mode line
.separator ", "


select count(*) as All_entries from bam_index;
select count(*) as Valid_entries from bam_index where isDuplicate=0 and isNotPassingQualityControls=0;
select count(*) as Duplicate from bam_index where isDuplicate=1 ;
select count(*) as Alignments from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 ;
select count(*) as Spliced from bam_index where nSpliced>=1 and isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1;
select count(*) as ReadsSpliced from bam_index where nSpliced>=1 and isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and isPrimary=1;
select count(*) as Paired from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isPaired=1 ;
select count(*) as Paired_mapped from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and isPaired=1 and isProperPair=1  ;
select count(*) as Paired_mapped_mate1 from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and isPaired=1 and isSecondMateRead=0 ;
select count(*) as Paired_mapped_mate2 from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and isPaired=1 and isSecondMateRead=1 ;
select count(*) as ReadsUmapped from bam_index where isMapped=0 and  isDuplicate=0 and isNotPassingQualityControls=0;
select count(*) as Reads from bam_index where  isPrimary=1;
select count(*) as NM0 from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and nm=0 ;
select count(*) as NM1 from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and nm=1 ;
select count(*) as NM2 from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and nm=2 ;
select count(*) as NMGE3 from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and nm>=3 ;
select count(*) as plus_strand from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and xs='+' ;
select count(*) as minus_strand from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and xs='-' ;
select count(*) as NH1 from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and nh=1 ;
select count(*) as NHGE2 from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and nh>=2 ;
select count(*) as NHGE10 from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and nh>=10 ;
select count(*) as NHGE20 from bam_index where isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1 and nh>=20 ;
select count(*) as Spliced1 from bam_index where nSpliced=1 and isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1;
select count(*) as Spliced2 from bam_index where nSpliced=2 and isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1;
select count(*) as Spliced3 from bam_index where nSpliced=3 and isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1;
select count(*) as SplicedGE4 from bam_index where nSpliced>3 and isDuplicate=0 and isNotPassingQualityControls=0  and isMapped=1;
EOF
exit 0
