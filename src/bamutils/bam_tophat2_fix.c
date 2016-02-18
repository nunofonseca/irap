/*
 * =========================================================
 * Copyright 2012-2016,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
 *    $Id: src/bamutils/bam_pe_insert.c 0.1.1 Nuno Fonseca Fri Dec 21 01:07:37 2012$
 * =========================================================
 */
#include <stdio.h>  
#include <math.h>
#include <bam.h>
#include <sam.h>
#include <kstring.h>      

#define VERSION "0.8.0p1"
#define MAX_INSERT_SIZE 1000
#define BIN_SIZE 10
     
// FIX MRNM and unaligned reads
int main(int argc, char *argv[])  
{  
  bamFile in; 
  long num_unmapped=0;
  long num_alns_pe=0;

  if (argc != 3) {  
    fprintf(stderr, "Usage: bam_tophat2_pe_fix <in.bam> <out.bam>\n");  
    return 1;  
  }  

  in = strcmp(argv[1], "-")? bam_open(argv[1], "rb") : bam_dopen(fileno(stdin), "rb"); 
  if (in == 0) {  
    fprintf(stderr, "ERROR: Fail to open input BAM file %s\n", argv[1]);  
    return 1;  
  }  
  

  int ref;  
  unsigned long num_alns=0;
  // counts
  unsigned long unalign_mapq_fix=0;
  unsigned long mtid_fix=0;
  unsigned long mpos_fix=0;
  bamFile out; 
  bam_header_t *header;
  header = bam_header_read(in);
  bam1_t *aln=bam_init1();

  out = strcmp(argv[2], "-")? bam_open(argv[2], "w") : bam_dopen(fileno(stdout), "w"); 
  if (out == 0) {  
    fprintf(stderr, "ERROR: Fail to open BAM file %s\n", argv[2]);  
    return 1;  
  }  
  bam_header_write(out,header);

  while(bam_read1(in,aln)>=0) {
    ++num_alns;
    if (aln->core.tid < 0) { 
      // unaligned reads
      if ( aln->core.qual!=0 ) {
	//fprintf(stderr, "ERROR: Unaligned read with quality > 0 in line %lu\n",num_alns);
	aln->core.qual=0;
	unalign_mapq_fix++;
      }
    }
    //fprintf(stderr,"%s %c %d %d\n",bam1_qname(aln),(aln->core.tid<0?'U':'M'),aln->core.mtid,aln->core.mpos);
    if ( aln->core.flag & BAM_FPAIRED ) {
      //fprintf(stderr,"paired %d\n",(aln->core.flag & BAM_FMUNMAP));
      // paired
      if ( aln->core.mtid <0 && !(aln->core.flag & BAM_FMUNMAP) ) {
	aln->core.flag |= BAM_FMUNMAP;
	aln->core.mpos=-1;
	mtid_fix++;
      }
      if ( aln->core.mpos <0 && !(aln->core.flag & BAM_FMUNMAP) ) {
	aln->core.flag |= BAM_FMUNMAP;
	aln->core.mtid=-1;
	mpos_fix++;
      }
    }
    
    bam_write1(out,aln);

  }
  bam_destroy1(aln);
  bam_close(in);  
  bam_close(out);  
  // 
  fprintf(stderr,"unaligned MAPQ fixes: %lu\n",unalign_mapq_fix);
  fprintf(stderr,"unaligned mtid fixes: %lu\n",mtid_fix);
  fprintf(stderr,"unaligned mpos fixes: %lu\n",mpos_fix);
  return 0;  
}  
