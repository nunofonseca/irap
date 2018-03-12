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
#include <time.h>
#include <assert.h>

#define BOOLEAN int
#define FALSE 0
#define TRUE 1

#define BACKLINE "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
#define PRINT_ALNS_PROCESSED(c) { if (c%100000==0) { fprintf(stderr,"%s%ld",BACKLINE,c);fflush(stderr); }}


int main(int argc, char *argv[])  
{  
  bamFile in; 
  char * sErrMsg = NULL;
  char * tail = 0;
  unsigned long nentries=1000000;
  int nRetCode;
  clock_t startClock,startClock2;

  if ( argc != 2 && argc!=3 ) {  
    fprintf(stderr, "Usage: bamStats <in.bam> [alignmentsToInspect@100000] \n");  
    return 1;  
  }  
  
  
  in = bam_open(argv[1], "rb");
  if (in == 0 ) {  
    fprintf(stderr, "ERROR: Fail to open BAM file %s\n", argv[1]);  
    return 1;  
  }  
  startClock = clock();
  
  // ***********
  // Read header
  bam_header_t *header;
  header = bam_header_read(in);
  // sorted by name?
  // Should not rely on the value in SO 
  bam1_t *aln=bam_init1();
  unsigned long num_alns=0;
  unsigned long uniq_alns=0;
  unsigned long mapped_read1=0;
  unsigned long reverse_read1=0;
  unsigned long reverse_read2=0;
  unsigned long forward_read1=0;
  unsigned long forward_read2=0;
  unsigned long plus_plus=0;
  unsigned long plus_minus=0;
  unsigned long minus_minus=0;
  unsigned long minus_plus=0;
  unsigned long SE_map=0;
  unsigned long SE_forward=0;
  unsigned long SE_reverse=0;


  while(bam_read1(in,aln)>=0 && num_alns<nentries) { // read alignment
    //aln->core.tid < 0 ? 
    uint8_t *nh = bam_aux_get(aln, "NH");
    uint8_t *xs = bam_aux_get(aln, "XS");    
    BOOLEAN isPrimary;
    BOOLEAN isDuplicate;
    BOOLEAN isNotPassingQualityControls;
    BOOLEAN isPaired;
    BOOLEAN isProperPair;
    // read is mapped to the reverse strand
    isPrimary= (aln->core.flag & BAM_FSECONDARY) ? FALSE:TRUE;
    isProperPair=(aln->core.flag & BAM_FPROPER_PAIR) ? TRUE:FALSE;
    isPaired=(aln->core.flag & BAM_FPAIRED ) ? TRUE:FALSE;
    // consider only the uniquely mapped reads
    if ( bam_aux2i(nh) <2) {                                                                                                                                                                  
      ++uniq_alns;                                                                                                                                                                            
      if ( isPaired && isProperPair ) {      // PE                                                                                                                                            
        if ( (aln->core.flag & 0x0040 ) != 0) {     //1st read                                                                                                                                
          if ( (aln->core.flag & 0x0004 ) == 0) {                                                                                                                                             
            ++mapped_read1;                                                                                                                                                                   
            if ( (aln->core.flag & 0x0010) != 0) ++reverse_read1;                                                                                                                             
            if ( (aln->core.flag & 0x0010 )== 0) ++forward_read1;                                                                                                                             
          }                                                                                                                                                                                   
        }                                                                                                                                                                                     
        if ((aln->core.flag & 0x0080) != 0) { // 2nd read                                                                                                                                     
          if ((aln->core.flag & 0x0004) == 0) {                                                                                                                                               
            if ( (aln->core.flag & 0x0010) != 0) ++reverse_read2;                                                                                                                             
                 if ( (aln->core.flag & 0x0010) == 0) ++forward_read2;                                                                                                                        
          }                                                                                                                                                                                   
        }                                                                                                                                                                                     
                                                                                                                                                                                              
        if (  (aln->core.flag & 0x0010 ) != 0 && ( aln->core.flag & 0x0020) != 0) minus_minus++;                                                                                              
        if ( (aln->core.flag & 0x0010 )!= 0 && (aln->core.flag & 0x0020) == 0) minus_plus++;                                                                                                  
        if ((aln->core.flag & 0x0010 ) == 0 && (aln->core.flag & 0x0020 ) != 0) plus_minus++;                                                                                                 
        if ( (aln->core.flag & 0x0010 ) == 0 && (aln->core.flag & 0x0020 ) == 0) plus_plus++;                                                                                                 
      }                                                                                                                                                                                       
      if (!isPaired) { // SE                                                                                                                                                                  
        if ( (aln->core.flag & 0x0004) == 0)  {                                                                                                                                               
          ++SE_map;                                                                                                                                                                           
          if ( ( aln->core.flag & 0x0010 ) != 0) ++SE_reverse;                                                                                                                                
          if ( ( aln->core.flag & 0x0010 ) == 0) ++SE_forward;                                                                                                                                
        }                                                                                                                                                                                     
      }                                                                                                                                                                                       
    }                                                                                                                                                                                         
    ++num_alns;                            
  }
   
  bam_close(in);  
  printf("\nProcessed %d records in %4.2f seconds\n", num_alns, ( (double) (clock() - startClock))/CLOCKS_PER_SEC);
  //
  printf("Alignments: %ld\n",num_alns);
  printf("Unique alignments: %ld (%4.2f)\n",num_alns,((double)uniq_alns/num_alns));
  printf("Mapped read 1 forward: %ld (%4.2f)\n",forward_read1,((double)forward_read1/(forward_read1+reverse_read1)));
  printf("Mapped read 1 reverse: %ld (%4.2f)\n",reverse_read1,((double)forward_read1/(forward_read1+reverse_read1)));

  printf("Mapped read 2 forward: %ld (%4.2f)\n",forward_read2,((double)forward_read1/(forward_read2+reverse_read2)));
  printf("Mapped read 2 reverse: %ld (%4.2f)\n",reverse_read2,((double)forward_read1/(forward_read2+reverse_read2)));


  printf("Mapped to +,-: %ld (%4.2f)\n",plus_minus,((double)plus_minus/(plus_plus+minus_plus+minus_minus+plus_minus)));
  printf("Mapped to -,+: %ld (%4.2f)\n",minus_plus,((double)minus_plus/(plus_plus+minus_plus+minus_minus+plus_minus)));
  printf("Mapped to +,+: %ld (%4.2f)\n",plus_plus,((double)plus_plus/(plus_plus+minus_plus+minus_minus+plus_minus)));
  printf("Mapped to -,-: %ld (%4.2f)\n",minus_minus,((double)minus_minus/(plus_plus+minus_plus+minus_minus+plus_minus)));

  printf("SE Mapped: %ld\n",SE_map);
  printf("Mapped to +: %ld (%4.2f)\n",SE_forward,((double)SE_forward/(SE_forward+SE_reverse)));
  printf("Mapped to -: %ld (%4.2f)\n",SE_reverse,((double)SE_forward/(SE_forward+SE_reverse)));
  return 0;  
}  
