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
 *    $Id: src/bamutils/bam_pe_insert.c 0.1.1 Nuno Fonseca Fri Dec 21 01:07:37 2012$
 * =========================================================
 */
#include <stdio.h>  
#include <math.h>
#include <bam.h>
#include <sam.h>
#include <kstring.h>      

#define VERSION "0.6.2d14"
#define MAX_INSERT_SIZE 1000
#define BIN_SIZE 10
     
int main(int argc, char *argv[])  
{  
  bamFile in; 

  if (argc == 1) {  
    fprintf(stderr, "Usage: bam_pe_insert <in.bam>\n");  
    return 1;  
  }  
  int k;
  long num_negative_vals=0;
  long num_alns=0;
  long num_alns_pe=0;
  unsigned int bins[MAX_INSERT_SIZE/BIN_SIZE+1];

  for (k=0;k<=MAX_INSERT_SIZE/BIN_SIZE;++k) bins[k]=0;
  in = bam_open(argv[1], "rb");  
  if (in == 0) {  
    fprintf(stderr, "Fail to open BAM file %s\n", argv[1]);  
    return 1;  
  }  
  int ref;  
  bam_header_t *header;
  header = bam_header_read(in);
  bam1_t *aln=bam_init1();
  while(bam_read1(in,aln)>=0) {
    ++num_alns;
    if (aln->core.tid < 0) continue;
    if (aln->core.flag & BAM_FSECONDARY) continue; // skip secondary alignments
    if (aln->core.flag & ! BAM_FPAIRED) continue; // not paired
    if (aln->core.flag & ! BAM_FPROPER_PAIR) continue; // not a proper pair
    if (aln->core.flag & ! BAM_FMUNMAP) continue; // the mate is mapped
    if (aln->core.flag & BAM_FSECONDARY) continue; // secundary read
    if (aln->core.flag & BAM_FREAD2) continue; // only count each pair once
    if (aln->core.isize<=0) ++num_negative_vals;
    else {
      //printf("%d %d bin=%d pos=%lu matepos=%lu\n",aln->core.flag,aln->core.isize,aln->core.isize/BIN_SIZE,     aln->core.pos,aln->core.mpos);
      if (aln->core.isize>MAX_INSERT_SIZE) {
	++bins[MAX_INSERT_SIZE/BIN_SIZE];
      } else ++bins[aln->core.isize/BIN_SIZE];
      ++num_alns_pe;
    }
  }
  // 
  printf("#Alignments: %lu    #PE alignments: %lu\n",num_alns,num_alns_pe*2);
  printf("#Insert size<0: %lu\n",num_negative_vals);
  if ( !num_alns_pe ) {
    exit(1);
  }
  unsigned long long tot=0;
  for (k=0;k<=MAX_INSERT_SIZE/BIN_SIZE;++k) {
    printf("%d\t%d\n",k*BIN_SIZE,bins[k]);
    tot+=k*BIN_SIZE*bins[k];
  }
  unsigned int avg=tot/num_alns_pe;
  tot=0;
  for (k=0;k<=MAX_INSERT_SIZE/BIN_SIZE;++k) {    
    tot+=((k*BIN_SIZE)-avg)*((k*BIN_SIZE)-avg)*bins[k];
    //if (bins[k]) printf("%ld %lu %lu %lu %ld\n",k*BIN_SIZE,avg,(k*BIN_SIZE),((k*BIN_SIZE)-avg)*((k*BIN_SIZE)-avg),bins[k]);
  }
  unsigned long var=tot/(num_alns_pe-1);
  long int stdev=sqrt(var);
  printf("Avg. insert size: %lu  stdev:%ld\n",avg,stdev);
  bam_destroy1(aln);
  bam_close(in);  
  return 0;  
}  
