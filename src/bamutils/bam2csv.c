/*
 * =========================================================
 * Copyright 2012-2018,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
 *
 * This file is part of iRAP.
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *    $Id$
 * =========================================================
 */
#include <stdio.h>  
#include <math.h>
#include <bam.h>
#include <sam.h>
#include <kstring.h>      

#define BOOLEAN int
#define FALSE 0
#define TRUE 1

#define BACKLINE "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
#define PRINT_ALNS_PROCESSED(c) { if (c%100000==0) { fprintf(stderr,"%s%ld",BACKLINE,c);fflush(stderr); }}


int main(int argc, char *argv[])  
{  
  bamFile in; 
  
  if (argc != 2) {  
    fprintf(stderr, "Usage: bam2csv <in.bam>\n");  
    return 1;  
  }  
  
  // Open file and exit if error
  //in = strcmp(argv[1], "-")? bam_open(argv[1], "rb") : bam_dopen(fileno(stdin), "rb");
  in = bam_open(argv[1], "rb");
  if (in == 0 ) {  
    fprintf(stderr, "ERROR: Fail to open BAM file %s\n", argv[1]);  
    return 1;  
  }  
  // ***********
  // Read header
  bam_header_t *header;
  header = bam_header_read(in);
  // sorted by name?
  // Should not rely on the value in SO 
  bam1_t *aln=bam_init1();
  unsigned long num_alns=0;

  while(bam_read1(in,aln)>=0) { // read alignment
    //aln->core.tid < 0 ? 
    uint8_t *nh = bam_aux_get(aln, "NH");
    uint8_t *nm = bam_aux_get(aln, "NM");
    uint8_t *xs = bam_aux_get(aln, "XS");
    
    BOOLEAN isPrimary;
    BOOLEAN isMapped;
    BOOLEAN isDuplicate;
    BOOLEAN isNotPassingQualityControls;
    BOOLEAN isPaired;
    BOOLEAN isSecondMateRead,isProperPair;
    //secondary alignment
    isMapped=((aln->core.flag & BAM_FMUNMAP) || aln->core.mtid <0) ? TRUE: FALSE;
    isPrimary= (aln->core.flag & BAM_FSECONDARY) ? FALSE:TRUE;
    isProperPair=(aln->core.flag & BAM_FPROPER_PAIR) ? TRUE:FALSE;
    isPaired=(aln->core.flag & BAM_FPAIRED ) ? TRUE:FALSE;
    isSecondMateRead=(aln->core.flag & BAM_FREAD2 ) ? TRUE: FALSE;
    isNotPassingQualityControls=(aln->core.flag & BAM_FQCFAIL ) ? TRUE:FALSE;
    isDuplicate=(aln->core.flag & BAM_FDUP) ? TRUE: FALSE;

    BOOLEAN isSpliced=FALSE;
    BOOLEAN hasSimpleCigar=TRUE;
    int i;
    if (aln->core.n_cigar != 0) {
      for (i = 0; i < aln->core.n_cigar; ++i) {
	char l="MIDNSHP=X"[bam1_cigar(aln)[i]&BAM_CIGAR_MASK];
	if ( l == 'S' ) { isSpliced=TRUE; hasSimpleCigar=FALSE;}	  
	if ( l != 'M' && l!='=' ) { isSpliced=TRUE; hasSimpleCigar=FALSE;}	  
      }
    }   
    // isDuplicate,isNotPassingQualityControls,
    // isSpliced,isPAired,isPrimary,hasSimpleCigar,isSecondMateRead,isProperPair,nh,nm,xs
    printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%c\n",
	   isDuplicate,isNotPassingQualityControls,
	   isSpliced,isPaired,isPrimary,isMapped,hasSimpleCigar,isSecondMateRead,isProperPair,
	   (nh==0?0:bam_aux2i(nh)),(nm==0?0:bam_aux2i(nm)),
	   (xs==0?' ':bam_aux2A(xs)));
    ++num_alns;
    PRINT_ALNS_PROCESSED(num_alns);
  }
  bam_close(in);  
  return 0;  
}  
