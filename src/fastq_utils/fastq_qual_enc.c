/*
 * =========================================================
 * Copyright 2012-2013,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>

#define PRINT_READS_PROCESSED(c) { if (c%1000000==0) { printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%ld",cline/4);fflush(stdout); }}

#define MAX_READ_LENGTH 1024000
char read_buffer[MAX_READ_LENGTH];
#define READ_LINE(fd) fgets(&read_buffer[0],MAX_READ_LENGTH,fd)

#define MAX_RANGE 127

 /* Source: http://en.wikipedia.org/wiki/FASTQ_format */
 /*  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS..................................................... */
 /*  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX...................... */
 /*  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII...................... */
 /*  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ...................... */
 /*  LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL.................................................... */
 /*  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~ */
 /*  |                         |    |        |                              |                     | */
 /* 33                        59   64       73                            104                   126 */

 /* S - Sanger        Phred+33,  raw reads typically (0, 40) */
 /* X - Solexa        Solexa+64, raw reads typically (-5, 40) */
 /* I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40) */
 /* J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40) */
 /*    with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)  */
 /*    (Note: See discussion above). */
 /* L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41) */


inline FILE* open_fastq(char* filename) {

  FILE *fd1=fopen(filename,"r");
  if (fd1==NULL) {
    fprintf(stderr,"Unable to open %s\n",filename);
    exit(1);
  }
  return(fd1);
}

int main(int argc, char **argv ) {

  //printf("%d",sizeof(struct index_entry)); 


  if (argc!=2) {
    fprintf(stderr,"Usage: fastq_qual_enc fastq1\n");
    exit(1);
  }
  FILE *fd1=open_fastq(argv[1]);
  // ************************************************************
  unsigned long cline=1;
  unsigned long cur_read=0;
  long num_reads=1000;
  short range[MAX_RANGE]; 
  for (cline=0; cline<MAX_RANGE; cline++) {
    range[cline]=0;
  }
  cline=1;
  // read the entry using another fd
  cline=1;
  while(!feof(fd1) && cur_read<num_reads) {
    char *hdr=READ_LINE(fd1);
    if ( hdr==NULL) break;
    if ( hdr[0]!='@' ) {
      fprintf(stderr,"line %ul: error in header %s",cline,hdr);
      return 1;
    }
    //SEQ
    READ_LINE(fd1);
    //@
    READ_LINE(fd1);
    // qual
    READ_LINE(fd1);
    // process qual
    long i=0;
    printf("%s",read_buffer);
    while ( read_buffer[i]!='\0' ) {
      short c=read_buffer[i];
      range[c]++;
      if ( range[c]<=0 ) { // overflow
	range[c]-=1;
      }
      printf("%c=%d\n",read_buffer[i],c);
      ++i;
    }
    cline+=4;
    cur_read++;
    
    PRINT_READS_PROCESSED(cline/4);
  }
  fclose(fd1);
  // 
  short min1,max1,min2=0,max2=0;
  for (cline=0; cline<MAX_RANGE; cline++) {
    if (range[cline]) 
      printf("%d|%d ",cline,range[cline]);
  }
  printf("\n");
  exit(0);
}
