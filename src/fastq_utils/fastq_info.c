/*
# =========================================================
# Copyright 2012-2017,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
# This file is part of iRAP.
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with iRAP.  If not, see <http://www.gnu.org/licenses/>.
#
#
#    $Id: 0.1.1$
# =========================================================
*/
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <regex.h> 
#include <zlib.h> 


#include "hash.h"
// 1MB
// disable this option if disk access is fast (local disk)
// enable it for network disks
#define VERSION "0.8.5.p8"


#define SEQDISKACCESS 1

#define MAX_READ_LENGTH 1024000
#define MIN_READ_LENGTH 1

#define HASHSIZE 19000001
#define TRUE 1
#define FALSE 0
#define UNDEF -1

#define DEFAULT  0
#define CASAVA18 1
#define INTEGERNAME 2

//#define HASHSIZE 5

#define READ_LINE(fd) gzgets(fd,&read_buffer[0],MAX_READ_LENGTH)
#define READ_LINE_HDR(fd) gzgets(fd,&read_buffer_hdr1[0],MAX_READ_LENGTH)
#define READ_LINE_HDR2(fd) gzgets(fd,&read_buffer_hdr2[0],MAX_READ_LENGTH)
#define READ_LINE_SEQ(fd) gzgets(fd,&read_buffer_seq[0],MAX_READ_LENGTH)
#define READ_LINE_QUAL(fd) gzgets(fd,&read_buffer_qual[0],MAX_READ_LENGTH)

#define READ_LINE_HDR2_1(fd) gzgets(fd,&read_buffer2_hdr1[0],MAX_READ_LENGTH)
#define READ_LINE_HDR2_2(fd) gzgets(fd,&read_buffer2_hdr2[0],MAX_READ_LENGTH)
#define READ_LINE_SEQ2(fd) gzgets(fd,&read_buffer2_seq[0],MAX_READ_LENGTH)
#define READ_LINE_QUAL2(fd) gzgets(fd,&read_buffer2_qual[0],MAX_READ_LENGTH)

#define PRINT_READS_PROCESSED(c) { if (c%500000==0) { fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%ld",c);fflush(stderr); }}

struct index_entry {
  
  // file offset: start entry
  // file offset: end entry
  // chat hdr(40)
  char *hdr;
  off_t entry_start;
  //unsigned  int  nbytes;
};
typedef struct index_entry INDEX_ENTRY;

// four buffers per file
char read_buffer_hdr1[MAX_READ_LENGTH];
char read_buffer_hdr2[MAX_READ_LENGTH];
char read_buffer_seq[MAX_READ_LENGTH];
char read_buffer_qual[MAX_READ_LENGTH];

char read_buffer2_hdr1[MAX_READ_LENGTH];
char read_buffer2_hdr2[MAX_READ_LENGTH];
char read_buffer2_seq[MAX_READ_LENGTH];
char read_buffer2_qual[MAX_READ_LENGTH];

char hdr[MAX_READ_LENGTH];

char read_buffer[MAX_READ_LENGTH];
long index_mem=0;

int readname_format=UNDEF;

#define VERSION "0.9.0a"


// approx. median read length
inline unsigned int median_rl(FASTQ_FILE* fd1,FASTQ_FILE* fd2) {
  unsigned long long ctr=0;
  unsigned int crl=1;
  unsigned long nreads=fd1->num_rds;
  
  if ( fd1->num_rds==1 && fd2==NULL) return(fd1->min_rl);
  if ( fd2!=NULL) nreads+=fd2->num_rds;
  while ( crl < MAX_READ_LENGTH ) {    
    ctr+=fd1->rdlen_ctr[crl];
    if (fd2!=NULL) ctr+=fd2->rdlen_ctr[crl];
    //printf("%d-%lu\n",crl,rdlen_ctr[crl]);
    if ( fd1->num_rds>1 && ctr>nreads/2) break;
    ++crl;
  }

  return(crl);
}
int has_gz_extension(const char *s) {
  regex_t regex;
  int reti;
  reti = regcomp(&regex,".gz$",0);  
  if ( reti ) { 
    fprintf(stderr, "\nInternal error: Could not compile regex\n"); 
    exit(2); 
  }
  /* Execute regular expression */
  //fprintf(stderr,"%s\n",hdr);
  reti = regexec(&regex, s, 0, NULL, 0);
  regfree(&regex);
  return(!reti);
}

FASTQ_FILE* validate_interleaved(char *f) {
  unsigned long cline=1;

  fprintf(stderr,"Paired-end interleaved\n");

  FASTQ_FILE* fd1=fastq_new(f,FALSE,"r");
  fastq_is_pe(fd1);
  FASTQ_ENTRY *m1=fastq_new_entry(),
    *m2=fastq_new_entry();

  char rname1[MAX_LABEL_LENGTH];
  char rname2[MAX_LABEL_LENGTH];
  unsigned long nreads1=0;
  unsigned long len=0;    

  while(!gzeof(fd1->fd)) {
    // read 1
    if (fastq_read_entry(fd1,m1)==0) break;
    // read 2
    if (fastq_read_entry(fd1,m2)==0) {
      fprintf(stderr,"\nError in file %s, line %lu: file truncated?\n",f,fd1->cline);
      exit(1);
    }
    // match
    char *readname1=fastq_get_readname(fd1,m1,&rname1[0],&len,TRUE);
    char *readname2=fastq_get_readname(fd1,m2,&rname2[0],&len,TRUE);

    // TODO
    // replace_dots(start_pos,seq1,hdr1,hdr1_2,qual1,fdf);    
    // replace_dots(start_pos,seq2,hdr2,hdr2_2,qual2,fdf);    

    if ( strcmp(readname1,readname2) ) {
      fprintf(stderr,"\nError in file %s, line %lu: unpaired read - %s\n",f,fd1->cline,readname1);
      exit(1);
    } 

    if (fastq_validate_entry(fd1,m1)) {
      exit(1);
    }
    if (fastq_validate_entry(fd1,m2)) {
      exit(1);
    }
    PRINT_READS_PROCESSED(cline/4,100000);
    nreads1+=2;
  }
  printf("\n");
  //close_fastq(fdf); ???
  //fastq_destroy(fd1);
  return(fd1);
}


int main(int argc, char **argv ) {
  //long paired=0;
  unsigned long num_reads1=0,
    num_reads2=0;

  unsigned long max_rl; // maximum read length
  unsigned long min_rl;  // minimum read length
  unsigned long min_qual; // minimum quality
  unsigned long max_qual;   // maximum quality
 
  int is_paired_data=FALSE;
  int is_interleaved=FALSE;
  //int fix_dot=FALSE;
  
  int nopt=0;
  int c;
  opterr = 0;

  fprintf(stderr,"Version iRAP %s\n",VERSION);
  
  while ((c = getopt (argc, argv, "f")) != -1)
    switch (c)
      {
      case 'f':
        //fix_dot = TRUE;
	fprintf(stderr,"Fixing (-f) enabled: Replacing . by N (creating .fix.gz files)\n");
	fprintf(stderr,"ERROR: -f option is no longer valid.\n");
	exit(1);
	++nopt;
        break;
      default:
	++nopt;
        fprintf(stderr,"ERROR: Option -%c invalid\n",optopt);
	exit(1);
      }
  
  if (argc-nopt<2 || argc-nopt>3) {
    fprintf(stderr,"Usage: fastq_info fastq1 [fastq2 file|pe]\n");
    //fprintf(stderr,"%d",argc);
    exit(1);
  }


  if (argc-nopt ==3) {
    is_paired_data=TRUE;
    //fprintf(stderr,"%d %d %d %s\n",argc,nopt,argc-nopt,argv[2+nopt]);
    if ( strncmp(argv[2+nopt],"pe",2) ) {
      is_interleaved=FALSE;
    } else {
      //fprintf(stderr,"Expecting interleaved reads...\n");
      is_interleaved=TRUE;
    }
  }

  FASTQ_FILE* fd1=NULL;
  FASTQ_FILE* fd2=NULL;
  hashtable index;
  // ************************************************************
  if ( is_interleaved ) {
    // interleaved    
    fd1=validate_interleaved(argv[1+nopt]);
  } else {
    // single or pair of fastq file(s)
    //unsigned long cline=1;
    fd1=fastq_new(argv[1+nopt],FALSE,"r");
    if ( is_paired_data) fastq_is_pe(fd1);

    fprintf(stderr,"HASHSIZE=%lu\n",(long unsigned int)HASHSIZE);
    index=new_hashtable(HASHSIZE);
    index_mem+=sizeof(hashtable);
    fprintf(stderr,"Scanning and indexing all reads from %s\n",fd1->filename);
    fastq_index_readnames(fd1,index,0,FALSE);
    fprintf(stderr,"Scanning complete.\n");

    num_reads1=index->n_entries;
    fprintf(stderr,"\n");
    // print some info
    fprintf(stderr,"Reads processed: %ld\n",index->n_entries);    
    fprintf(stderr,"Memory used in indexing: ~%ld MB\n",index_mem/1024/1024);
  }
  min_rl=fd1->min_rl;
  max_rl=fd1->max_rl;
  min_qual=fd1->min_qual;
  max_qual=fd1->max_qual;

  // pair-end
  if (argc-nopt ==3 && !is_interleaved ) {
    fprintf(stderr,"File %s processed\n",argv[1+nopt]);  
    fprintf(stderr,"Next file %s\n",argv[2+nopt]);  
    // validate the second file and check if all reads are paired
    fd2=fastq_new(argv[2+nopt],FALSE,"r");
    fastq_is_pe(fd2);
    
    unsigned long len;
    FASTQ_ENTRY *m2=fastq_new_entry();
    char rname[MAX_LABEL_LENGTH];
    // 
    while(!gzeof(fd2->fd)) {
      // read entry
      if (fastq_read_entry(fd2,m2)==0) break;
      char *readname=fastq_get_readname(fd2,m2,&rname[0],&len,TRUE);
      INDEX_ENTRY* e=fastq_index_lookup_header(index,readname);
      if (e==NULL) {
	// complain and exit if not found
	fprintf(stderr,"\nError in file %s, line %lu: unpaired read - %s\n",argv[2+nopt],fd2->cline,readname);
	exit(1);
      }
      fastq_index_delete(readname,index);
      // TODO
      //replace_dots(start_pos,seq,hdr,hdr2,qual,fdf);
      PRINT_READS_PROCESSED(fd2->cline/4,100000);
    }
    printf("\n");
    //fastq_destroy(fdf);//???
    if (index->n_entries>0 ) {
      fprintf(stderr,"\nError in file %s: found %lu unpaired reads\n",argv[1+nopt],index->n_entries);
      exit(1);
    }
    // stats
    min_rl=min(fd2->min_rl,min_rl);
    max_rl=max(fd2->max_rl,max_rl);
    min_qual=min(fd2->min_qual,min_qual);
    max_qual=max(fd2->max_qual,max_qual);    
  }

  // stats
  // min qual/max qual/read len
  FILE* out;  
  out=stderr;

  fprintf(out,"------------------------------------\n");
  if ( num_reads2>0 ) {
    fprintf(out,"Number of reads: %lu %lu\n",num_reads1,num_reads2);
  } else {
    fprintf(out,"Number of reads: %lu\n",num_reads1);
  }
  fprintf(out,"Quality encoding range: %lu %lu\n",min_qual,max_qual);
  char *enc=fastq_qualRange2enc(min_qual,max_qual);
  if ( enc == NULL ) {
    fprintf(stderr,"\nERROR: Unable to determine quality encoding - unknown range [%lu,%lu]\n",min_qual,max_qual);
    exit(1);
  }
  fprintf(out,"Quality encoding: %s\n",enc);
  fprintf(out,"Read length: %lu %lu %u\n",min_rl,max_rl,median_rl(fd1,fd2));
  fprintf(out,"OK\n"); 
  exit(0);
}

