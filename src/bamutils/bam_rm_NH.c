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
 *    $Id$
 * =========================================================
 */
#include <stdio.h>  
#include <math.h>
#include <bam.h>
#include <sam.h>
#include <kstring.h>      

#include "../fastq_utils/hash.h"
#include <string.h>
#include <assert.h>

// Same version as IRAP
#define VERSION "0.8.0d4"


#define BACKLINE "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
#define PRINT_ALNS_PROCESSED(c) { if (c%100000==0) { printf("%s%ld",BACKLINE,c);fflush(stdout); }}


int main(int argc, char *argv[])  
{  
  short out2stdout=0;
  bamFile in,in2; 
  bamFile out; 
  int paired;//1 if not paired or pair read 1, 2 otherwise

  if (argc != 3) {  
    fprintf(stderr, "Usage: bam_rm_NH <in.bam> <out.bam or - for stdout>\n");  
    return 1;  
  }  
  // Open file and exit if error
  in = bam_open(argv[1], "rb");
  out2stdout = strcmp(argv[2], "-")? 0 : 1; 
  out = strcmp(argv[2], "-")? bam_open(argv[2], "w") : bam_dopen(fileno(stdout), "w"); 
  if (in == 0 ) {  
    fprintf(stderr, "ERROR: Fail to open BAM file %s\n", argv[1]);  
    return 1;  
  }  
  if (out == 0) {  
    fprintf(stderr, "ERROR: Fail to open BAM file %s\n", argv[2]);  
    return 1;  
  }  

  unsigned long num_alns=0;
  unsigned long num_nh=0;
  int ref;  

  // ***********
  // Copy header
  bam_header_t *header;
  header = bam_header_read(in);
  bam_header_write(out,header);

  // sorted by name?
  // Should not rely on the value in SO 
  bam1_t *aln=bam_init1();
  bam1_t *prev=bam_init1();

  if (!out2stdout) {
    printf("bam_fix_NH version %s\n",VERSION);
    printf("Processing %s\n",argv[1]);
  }

  // reopen
  in2 = bam_open(argv[1], "rb");
  if (in2 == 0 ) {  
    fprintf(stderr, "ERROR: Fail to open BAM file %s\n", argv[1]);  
    return 1;  
  }  

  header = bam_header_read(in2);
  num_alns=0;
  num_nh=0;
  while(bam_read1(in2,aln)>=0) { // read alignment
    paired=1;
    if (aln->core.tid < 0) continue;//ignore unaligned reads
    if (aln->core.flag & BAM_FUNMAP) continue;
    if (aln->core.flag & BAM_FREAD2) paired=2;
    ++num_alns;

    // update the NH field
    uint8_t *old_nh = bam_aux_get(aln, "NH");    
    if (old_nh) {
      ++num_nh;
      bam_aux_del(aln, old_nh);
#ifdef DEBUG
      //      printf("!>%s %d\n",bam1_qname(aln),r->ctr);
#endif
    }
    bam_write1(out,aln);
    if(!out2stdout) PRINT_ALNS_PROCESSED(num_alns);
  }
  // 
  bam_destroy1(aln);
  bam_close(in2);  
  bam_close(out);  
  if(!out2stdout) {
    printf("%s%lu %lu\n",BACKLINE,num_alns,num_nh);
    printf("Done.\n");
  }
  return 0;  
}  
