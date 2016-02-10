/*
# =========================================================
# Copyright 2012-2016,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
#    $Id: scripts/irap Nuno Fonseca Wed Jan 2 10:16:41 2013$
# =========================================================
*/
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <time.h>


#define PRINT_READS_PROCESSED(c) { if (c%1000000==0) { fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%lu",cline/4);fflush(stderr); }}

#define GET_RAND(max)  (max>RAND_MAX?rand()*1000000+rand():rand())

#define MAX_READ_LENGTH 1024000
char read_buffer[MAX_READ_LENGTH];
#define READ_LINE(fd) fgets(&read_buffer[0],MAX_READ_LENGTH,fd)

#ifndef DEBUG
#define DEBUG 0
#endif

#define PAIRED(void) (globalArgs.fastq_filename2!=NULL)

/* ************************************************************************** */
/*****************************************************************************/
// Bits operations
#define BITMAP_empty(b)                ((b) == 0)
#define BITMAP_member(b,n)             (((b) & (1<<(n))) != 0)
#define BITMAP_alone(b,n)              ((b) == (1<<(n)))
#define BITMAP_subset(b1,b2)           (((b1) & (b2)) == b2)
#define BITMAP_same(b1,b2)             ((b1) == (b2))

#define BITMAP_on_all(b)                  ((b) = 255)

#define BITMAP_clear(b)                ((b) = 0)
#define BITMAP_and(b1,b2)              ((b1) &= (b2))
#define BITMAP_minus(b1,b2)            ((b1) &= ~(b2))
#define BITMAP_insert(b,n)             ((b) |= (1<<(n)))
#define BITMAP_delete(b,n)             ((b) &= (~(1<<(n))))
#define BITMAP_copy(b1,b2)             ((b1) = (b2))
#define BITMAP_intersection(b1,b2,b3)  ((b1) = ((b2) & (b3)))
#define BITMAP_difference(b1,b2,b3)    ((b1) = ((b2) & (~(b3))))

#define BMAP_BYTE(number) (number-1)/8
#define BMAP_BIT(number)  (number-1)%8
#define BMAP_SIZE(number) ((number)/8+(number%8>0?1:0))
#define VALID_NUMBER(bmap,number) (number>0 && number<=bm->max)


#ifndef NUM
#define NUM unsigned long
#endif

#ifdef TRUE
typedef  short BOOLEAN;
#else
typedef enum { FALSE=0, TRUE=1 } BOOLEAN;
#endif
typedef  enum { IN=1, OUT=0} STATUS;

typedef struct {
  NUM  max;      // max number <=Size
  NUM  size;     // array size
  char *array;
} BMAP;

BMAP* new_bmap(const NUM max_size) {
  NUM bmap_size=0L; // size of the last bitmap created
  BMAP *new;
  char *p;

  if ( max_size < 1  ) return NULL;

  bmap_size=BMAP_SIZE(max_size);
  //if (max_size%8>0) bmap_size++;
  new=(BMAP*)malloc(sizeof(BMAP)+bmap_size);
  if(new==NULL)
    return NULL;
  p=(char*)new;
  new->array=(char*)(p+sizeof(BMAP));

  new->size=bmap_size;
  new->max=max_size;
  memset(new->array,OUT,new->size); // empty bitmap

  return new;
}

inline
void free_bmap(BMAP *bm) {
  free(bm);
}

BOOLEAN set_in_bmap(BMAP *bm, const NUM number, const STATUS status) {

  NUM   byte=BMAP_BYTE(number);
  short bit=BMAP_BIT(number);
  char  *storage=&bm->array[byte];
  
  if ( ! VALID_NUMBER(bm,number) )  return FALSE;

  if ( status==IN )
    BITMAP_insert(*storage,bit);
  else
    BITMAP_delete(*storage,bit);
  return TRUE;
}

BOOLEAN  in_bmap(const BMAP *bm, const NUM number) {
  NUM  byte=BMAP_BYTE(number);
  short bit=BMAP_BIT(number);

  if ( ! VALID_NUMBER(bm,number) ) return FALSE;
  return BITMAP_member(bm->array[byte],bit);
}

/*****************************************************************************/

char c_time_string[100];
inline char* timestamp(void) {
  time_t current_time;
 
  /* Obtain current time as seconds elapsed since the Epoch. */
  current_time = time(NULL);
  //c_time_string = ctime(&current_time);
  struct tm *tm=localtime(&current_time);
  sprintf(&c_time_string[0],"%d/%d %d:%d:%d",tm->tm_mday,tm->tm_mon+1,tm->tm_hour,tm->tm_min,tm->tm_sec);
  return &c_time_string[0];
}

inline FILE* open_fastq(char* filename) {

  FILE *fd1=fopen(filename,"r");
  if (fd1==NULL) {
    fprintf(stderr,"[ERROR] Unable to open %s\n",filename);
    exit(1);
  }
  return(fd1);
}

long get_read_length(const char *r) {
  long n=strlen(r);
  //fprintf(stderr,">>>%c<<<<<\n",r[n-1]);
  if ( r[n-1]=='\n' ) {
    return(n-1);
  }
  return(n);
}

// 
long print_read(FILE *from,FILE* to) {
  long rl;
  fputs(READ_LINE(from),to);
  fputs(READ_LINE(from),to);
  rl=get_read_length(read_buffer);
  fputs(READ_LINE(from),to);
  fputs(READ_LINE(from),to);
  return rl;
}

// count the number of reads in a fastq file
inline long num_reads(char *fastq_file) {

  FILE *fd1=open_fastq(fastq_file);
  long n_reads=0;
  long cline=1;
  char *hdr;
  while(!feof(fd1)) {
    hdr=READ_LINE(fd1);
    if ( hdr==NULL) break;
    if ( hdr[0]!='@' ) {
      fprintf(stderr,"line %lu: error in header %s",cline,hdr);
      return 1;
    }
    READ_LINE(fd1);
    READ_LINE(fd1);
    READ_LINE(fd1);
    cline+=4;
    ++n_reads;
  }
  return n_reads;
}

struct globalArgs_t {
  long sample_size;                /* -n option */ 
  int  seed;                       /* -s option */
  char *fastq_filename1;           /* f1 option */
  char *fastq_filename2;           /* f2 option */
  short extra_info;                /* -x option */
  char *outfile_prefix;            /* -o option */
  char out1[1024]; 
  char out2[1024]; 
} globalArgs;


void display_usage() {
   fprintf(stderr,"Usage: fastq_sample -1 fastq_file [-2 fastq_file] [-s seed] [-n num_reads_in_sample] [-o file_prefix] [-h] [-x]\n");
   fprintf(stderr,"\n");
} 

// fastq_randsample -n num_reads [-s seed] -1 fastq -2 fastq -h -x
int main(int argc, char **argv ) {

  static const char *optString = "o:n:s:1:2:xh?";
  static const struct option longOpts[] = {
    { "out", required_argument, NULL, 'o' },
    { "extra_info", no_argument, NULL, 'x' },
    { "num_reads", required_argument, NULL, 'n' },
    { "seed", required_argument, NULL, 's' },
    { "fastq_file1", required_argument, NULL, '1' },
    { "fastq_file2", required_argument, NULL, '2' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
  };
  int opt = 0;
  int longIndex = 0;

  /* Initialize globalArgs before we get to work. */
  globalArgs.sample_size = 1000000;
  globalArgs.extra_info = 0; /* false */
  globalArgs.seed = (int)time(NULL); 
  globalArgs.fastq_filename1=NULL;
  globalArgs.fastq_filename2=NULL;
  globalArgs.outfile_prefix=NULL;
  globalArgs.out1[0]='\0';
  globalArgs.out2[0]='\0';

  //  while( (opt = getopt( argc, argv, optString )) != -1 ) {
  while( (opt = getopt_long( argc, argv, optString, longOpts, &longIndex )) != -1 ) {
    switch( opt ) {
    case 'o':
      globalArgs.outfile_prefix=optarg;
      break;

    case 'x':
      globalArgs.extra_info++;
      break;

    case 'n':
      globalArgs.sample_size = atol(optarg); 
      break;
      
    case 's':
      globalArgs.seed = atoi(optarg);
      break;
      
    case '1':
      globalArgs.fastq_filename1 = optarg;
      break;

    case '2':
      globalArgs.fastq_filename2 = optarg;
      break;
      
    case 'h':   /* fall-through is intentional */
    case '?':
      display_usage();
      exit(0);
      break;
      
    default:
      break;
    }        
    
  }
    
  if ( globalArgs.fastq_filename1==NULL ) {
    display_usage();
    exit(1);
  }
  // *****************************************
  // print settings
  fprintf(stderr," [INFO %s] Sample_size=%lu\n",timestamp(),globalArgs.sample_size);
  fprintf(stderr," [INFO %s] Seed=%d\n",timestamp(),globalArgs.seed);  
  fprintf(stderr," [INFO %s] FASTQ file=%s\n",timestamp(),globalArgs.fastq_filename1);

  if ( PAIRED() ) {
    fprintf(stderr," [INFO %s] Assuming PE data with fastq ordered by name\n",timestamp());
    fprintf(stderr," [INFO %s] FASTQ (2) file=%s\n",timestamp(),globalArgs.fastq_filename2);
    if (globalArgs.outfile_prefix!=NULL) {
      sprintf(&globalArgs.out1[0],"%s_1.fastq",globalArgs.outfile_prefix);
      sprintf(&globalArgs.out2[0],"%s_2.fastq",globalArgs.outfile_prefix);
      fprintf(stderr," [INFO %s] writing sample to files %s and %s\n",timestamp(),globalArgs.out1,globalArgs.out2);
    }
  } else {
    if (globalArgs.outfile_prefix!=NULL) {
      sprintf(&globalArgs.out1[0],"%s.fastq",globalArgs.outfile_prefix);
      fprintf(stderr," [INFO %s] writing sample to file %s\n",timestamp(),globalArgs.out1);
    }
  }
  
  srandom(globalArgs.seed);

  fprintf(stderr," [INFO %s] Reading %s...\n",timestamp(),globalArgs.fastq_filename1);
  long n_reads=num_reads(globalArgs.fastq_filename1);
  fprintf(stderr," [INFO %s] Reading %s...done (%ld reads).\n",timestamp(),globalArgs.fastq_filename1,n_reads);

  if (PAIRED()) {
    fprintf(stderr," [INFO %s] Reading %s...\n",timestamp(),globalArgs.fastq_filename2);
    long n_reads2=num_reads(globalArgs.fastq_filename2);
    fprintf(stderr," [INFO %s] Reading %s...done (%ld reads).\n",timestamp(),globalArgs.fastq_filename2,n_reads2);
    if ( n_reads2!=n_reads ) {
      fprintf(stderr," [ERROR] Number of reads are not the same on both files (%ld!= %ld )\n",n_reads,n_reads2);
      exit(1);
    }
  }

  if ( n_reads<globalArgs.sample_size ) {
    fprintf(stderr," [ERROR] Sample size greater than the number of reads  (%ld > %ld )\n",globalArgs.sample_size,n_reads);
    exit(1);
  }

  fprintf(stderr," [INFO %s] Starting random selection of %ld reads...\n",timestamp(),globalArgs.sample_size);
  double p=globalArgs.sample_size*1.0/n_reads; // probability of selecting a read  
  fprintf(stderr," [INFO %s] p=%f \n",timestamp(),p);
  BMAP *bit=new_bmap(n_reads);
  unsigned long selected=globalArgs.sample_size,
    ctr=1;
  //
  while (selected ) {
    double r=random()*1.0/RAND_MAX;  
    if (r <= p ) {
      set_in_bmap(bit,ctr,IN);
      --selected;
    }
    ctr++;
    if (ctr>n_reads) {
      ctr=1;
      p*=1.01; 
      fprintf(stderr,"|");
      fprintf(stderr,"%f%%",1-selected*1.0/globalArgs.sample_size);
    }
  }

  fprintf(stderr," [INFO %s] Starting random selection of %ld reads complete.\n",timestamp(),globalArgs.sample_size);
  // min, max and avg. read length
  FILE *fd1=open_fastq(globalArgs.fastq_filename1);
  FILE *fd2=NULL;
  
  FILE *fout1=NULL;
  FILE *fout2=NULL;
  
  if ( globalArgs.outfile_prefix!=NULL) {
      fout1=fopen(globalArgs.out1,"w");
      if (fout1==NULL) {
	fprintf(stderr,"[ERROR] Unable to open %s\n",globalArgs.out1);
	exit(1);
      }
      if (PAIRED()) {
	fout2=fopen(globalArgs.out2,"w");
	if (fout2==NULL) {
	  fprintf(stderr,"[ERROR] Unable to open %s\n",globalArgs.out2);
	  exit(1);
	}
      }
  } else {
    fout1=stdout;
    fout2=stdout;    
  }
  
  long key;  
  long cur_read2,cur_read=1;
  long cline2,cline=1;
  // stats
  long total_read_length=0, 
    min_read_length=1000000, 
    max_read_length=0, 
    info_num_reads=0;

  if (PAIRED()) {
    fd2=open_fastq(globalArgs.fastq_filename2);
    cur_read2=1;
    cline2=1;
  }
  
  ctr=1;
  //
  while (ctr <= n_reads ) {
    if ( in_bmap(bit,ctr) ) {
    key=ctr;
    if (DEBUG) fprintf(stderr,"picked:%ld\n",key);
    // 
    if (PAIRED()) {
      while(!feof(fd2) && cur_read2<key) { 
	char *hdr=READ_LINE(fd2);
	if ( hdr==NULL) break;
	if ( hdr[0]!='@' ) {
	  fprintf(stderr,"line %lu: error in header %s",cline2,hdr);
	  return 1;
	}
	char *r=READ_LINE(fd2); // read
	if ( globalArgs.extra_info ) {
	  long rl=get_read_length(r);
	  //fprintf(stderr, ">%ld\n",rl);
	  total_read_length+=rl;
	  if ( rl > max_read_length) max_read_length=rl;
	  if ( rl < min_read_length) min_read_length=rl;
	  ++info_num_reads;
	}
	READ_LINE(fd2); // quality
	READ_LINE(fd2); // quaity
	cline2+=4;
	cur_read2++;
	PRINT_READS_PROCESSED(cline2/4); 
      }      
    }
    while(!feof(fd1) && cur_read<key) { 
      char *hdr=READ_LINE(fd1);
      if ( hdr==NULL) break;
      if ( hdr[0]!='@' ) {
	fprintf(stderr,"line %lu: error in header %s",cline,hdr);
	return 1;
      }
      char *r=READ_LINE(fd1); // read
      if ( globalArgs.extra_info ) {
	long rl=get_read_length(r);
	//fprintf(stderr, ">%ld\n",rl);
	total_read_length+=rl;
	if ( rl > max_read_length) max_read_length=rl;
	if ( rl < min_read_length) min_read_length=rl;
	++info_num_reads;
      }
      READ_LINE(fd1); // quality
      READ_LINE(fd1); // quaity
      cline+=4;
      cur_read++;
      PRINT_READS_PROCESSED(cline/4); 
    }
    if(cur_read==key) {
      long rl=print_read(fd1,fout1);
      if ( globalArgs.extra_info ) {
	if (DEBUG) fprintf(stderr, ">%ld\n",rl);
	total_read_length+=rl;
	if ( rl > max_read_length) max_read_length=rl;
	if ( rl < min_read_length) min_read_length=rl;
	++info_num_reads;
      }
      cur_read++;
      cline+=4;
    }
    // PE
    if (PAIRED()) {
      if(cur_read2==key) {
	long rl=print_read(fd2,fout2);
	if ( globalArgs.extra_info ) {
	  if (DEBUG) fprintf(stderr, ">%ld\n",rl);
	  total_read_length+=rl;
	  if ( rl > max_read_length) max_read_length=rl;
	  if ( rl < min_read_length) min_read_length=rl;
	  ++info_num_reads;
	}
	cur_read2++;
	cline2+=4;
      }
    }
    }
    ++ctr;
  }  
  if ( globalArgs.extra_info ) {
    fprintf(stderr,"Read length\n");
    fprintf(stderr,"Max: %ld\n",max_read_length);
    fprintf(stderr,"Min: %ld\n",min_read_length);
    fprintf(stderr,"Avg: %f\n",total_read_length*1.0/info_num_reads);
  }
  exit(0);
}
