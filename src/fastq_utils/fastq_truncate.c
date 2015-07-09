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

#define PRINT_READS_PROCESSED(c) { if (c%1000000==0) { fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%lu",cline/4);fflush(stderr); }}


#define MAX_READ_LENGTH 1024000
char read_buffer[MAX_READ_LENGTH];
#define READ_LINE(fd) fgets(&read_buffer[0],MAX_READ_LENGTH,fd)

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
  
  if (argc!=3) {
    fprintf(stderr,"Usage: fastq_truncate fastq1 num_reads\n");
    exit(1);
  }
  FILE *fd1=open_fastq(argv[1]);
  long num_reads=atol(argv[2]);
  // ************************************************************
  unsigned long cline=1;
  unsigned long cur_read=0;

  // read the entry using another fd
  cline=1;
  while(!feof(fd1) && cur_read<num_reads) {
    char *hdr=READ_LINE(fd1);
    if ( hdr==NULL) break;
    if ( hdr[0]!='@' ) {
      fprintf(stderr,"line %lu: error in header %s",cline,hdr);
      return 1;
    }
    READ_LINE(fd1);
    READ_LINE(fd1);
    READ_LINE(fd1);
    cline+=4;
    cur_read++;
    
    PRINT_READS_PROCESSED(cline/4);
  }
  if ( cur_read==num_reads) {
    long pos=ftell(fd1);
    fclose(fd1);
    truncate(argv[1],pos);
  } else {
    fclose(fd1);
  }
  exit(0);
}
