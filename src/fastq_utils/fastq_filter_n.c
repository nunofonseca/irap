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

#define MAX_READ_LENGTH 1024000
char read_buffer[5][MAX_READ_LENGTH];
#define READ_LINE(fd,pos) fgets(&read_buffer[pos][0],MAX_READ_LENGTH,fd)
#define WRITE_READ(fd) {fputs(read_buffer[1],fd);fputs(read_buffer[2],fd);fputs(read_buffer[3],fd);fputs(read_buffer[4],fd);}


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
    fprintf(stderr,"Usage: fastq_filter_n fastq1\n");
    exit(1);
  }
  FILE *fd1=open_fastq(argv[1]);
  // ************************************************************
  unsigned long cline=1;
  unsigned long cur_read=0;
  char tmp_buffer[MAX_READ_LENGTH];
  // read the entry using another fd
  cline=1;
  while(!feof(fd1)) {
    char *hdr=READ_LINE(fd1,1);
    if ( hdr==NULL) break;
    if ( hdr[0]!='@' ) {
      fprintf(stderr,"line %ul: error in header %s",cline,hdr);
      return 1;
    }
    //
    char *seq=READ_LINE(fd1,2);
    READ_LINE(fd1,3);
    READ_LINE(fd1,4);
    
    short n_found=0;
    int k;
    for ( k=0;k<MAX_READ_LENGTH;k++) {
      if (seq[k]=='\n') break;
      if (seq[k]=='N' || seq[k]=='n' ) {
	n_found=1; break;
      }
    }
    if ( ! n_found ) 
      WRITE_READ(stdout);
    cline+=4;
    cur_read++;    
  }
  fclose(fd1);
  exit(0);
}
