/*
 * =========================================================
 * Copyright 2012-2017,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
 * =========================================================
 */
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>

#include "fastq.h"
#define VERSION "0.9.0a"


int main(int argc, char **argv ) {

  int nopt=0; 
  int c;
  opterr = 0;

  fprintf(stderr,"Version iRAP %s\n",VERSION); 
  // add an option -n N
  //if (optopt == 'c')
  char *cvalue = NULL;
  unsigned max_n=0;
  while ((c = getopt (argc, argv, "n:")) != -1)
    switch (c) 
      {
      case 'n':
	cvalue = optarg;
	max_n=atoi(cvalue);
	if ( max_n > 100 ) max_n=100;
  	nopt+=2;
        break;
      default:
  	++nopt;
        fprintf(stderr,"ERROR: Option -%c invalid\n",optopt);
  	exit(1);
      }
  
  if (argc-nopt<2 || argc-nopt>3) { 
    fprintf(stderr,"Usage: fastq_filter_n [ -n 0 ] fastq1\n");
    exit(1);
  }
  
  fprintf(stderr,"Discard reads with more than %d%% of Ns\n",max_n);
  
  FASTQ_FILE *fd1=fastq_new(argv[nopt+1],FALSE,"r");

  unsigned num_n,max_num_n;
  FASTQ_ENTRY *m1=fastq_new_entry();

  while(!gzeof(fd1->fd)) {
    num_n=0;
    if (fastq_read_entry(fd1,m1)==0) break;
    
    int k;
    max_num_n=m1->read_len*max_n/100;
    for ( k=0;k<MAX_READ_LENGTH;k++) {
      if (m1->seq[k]=='\n' || m1->seq[k]=='\0') break;
      if (m1->seq[k]=='N' || m1->seq[k]=='n'  ) {
	++num_n;
	if ( num_n > max_num_n ) break;
      }
    }
    // 
    if ( num_n <= max_num_n  ) {
      fastq_write_entry2stdout(m1);
    }
    PRINT_READS_PROCESSED(fd1->cline,100000);
  }
  fastq_destroy(fd1);
  exit(0);
}
