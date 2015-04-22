/*
 * =========================================================
 * Copyright 2012-2015,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
 *    $Id: src/bamutils/bam_fix_NH.c 0.1.1 Nuno Fonseca Fri Dec 21 01:07:37 2012$
 * =========================================================
 */
#include <stdio.h>  
#include <math.h>
#include <bam.h>
#include <sam.h>
#include <kstring.h>      

#include <string.h>
#include <assert.h>

#define VERSION "0.6.1p15"


#define BACKLINE "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
#define PRINT_ALNS_PROCESSED(c) { if (c%100000==0) { fprintf(stderr,"%s%ld",BACKLINE,c);fflush(stderr); }}


int main(int argc, char *argv[])  
{  
  short out2stdout=0;
  bamFile in,in2; 
  bamFile out; 


  if (argc != 3) {  
    fprintf(stderr, "Usage: bam_fix_se_flag <in.bam> <out.bam or - for stdout>\n");  
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
    fprintf(stderr,"bam_fix_se_flag version %s\n",VERSION);
    fprintf(stderr,"Processing %s\n",argv[1]);
  }

  // reopen
  in2 = bam_open(argv[1], "rb");
  if (in2 == 0 ) {  
    fprintf(stderr, "ERROR: Fail to open BAM file %s\n", argv[1]);  
    return 1;  
  }  

  header = bam_header_read(in2);
  num_alns=0;
  while(bam_read1(in2,aln)>=0) { // read alignment
    if (aln->core.tid < 0) continue;//ignore unaligned reads
    if (aln->core.flag & BAM_FUNMAP) continue;
    if (aln->core.flag & BAM_FPAIRED ) { // PAIRED

    } else { //SE 
      //turn off the other pair related flags
      aln->core.flag&=~BAM_FPROPER_PAIR;
      aln->core.flag&=~BAM_FMUNMAP;
      aln->core.flag&=~BAM_FREAD1;
      aln->core.flag&=~BAM_FREAD2;
      fprintf(stderr, ".");  
    }
    bam_write1(out,aln);
    if(!out2stdout) PRINT_ALNS_PROCESSED(num_alns);
    ++num_alns;
  }
  // 
  bam_destroy1(aln);
  bam_close(in2);  
  bam_close(out);  
  if(!out2stdout) {
    fprintf(stderr,"%s%lu\n",BACKLINE,num_alns);
    fprintf(stderr,"Done.\n");
  }
  return 0;  
}  
