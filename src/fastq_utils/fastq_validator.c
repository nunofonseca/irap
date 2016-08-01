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


#include "hash.h"
// 1MB
// disable this option if disk access is fast (local disk)
// enable it for network disks
#define VERSION "0.8.1d8"

#define SEQDISKACCESS 1

#define MAX_READ_LENGTH 1024000
#define MIN_READ_LENGTH 1

#define HASHSIZE 19000001
//#define HASHSIZE 5

#define READ_LINE(fd) fgets(&read_buffer[0],MAX_READ_LENGTH,fd)
#define READ_LINE_HDR(fd) fgets(&read_buffer_hdr1[0],MAX_READ_LENGTH,fd)
#define READ_LINE_HDR2(fd) fgets(&read_buffer_hdr2[0],MAX_READ_LENGTH,fd)
#define READ_LINE_SEQ(fd) fgets(&read_buffer_seq[0],MAX_READ_LENGTH,fd)
#define READ_LINE_QUAL(fd) fgets(&read_buffer_qual[0],MAX_READ_LENGTH,fd)

#define READ_LINE_HDR2_1(fd) fgets(&read_buffer2_hdr1[0],MAX_READ_LENGTH,fd)
#define READ_LINE_HDR2_2(fd) fgets(&read_buffer2_hdr2[0],MAX_READ_LENGTH,fd)
#define READ_LINE_SEQ2(fd) fgets(&read_buffer2_seq[0],MAX_READ_LENGTH,fd)
#define READ_LINE_QUAL2(fd) fgets(&read_buffer2_qual[0],MAX_READ_LENGTH,fd)

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

// four buffers
char read_buffer_hdr1[MAX_READ_LENGTH];
char read_buffer_hdr2[MAX_READ_LENGTH];
char read_buffer_seq[MAX_READ_LENGTH];
char read_buffer_qual[MAX_READ_LENGTH];

char read_buffer2_hdr1[MAX_READ_LENGTH];
char read_buffer2_hdr2[MAX_READ_LENGTH];
char read_buffer2_seq[MAX_READ_LENGTH];
char read_buffer2_qual[MAX_READ_LENGTH];

char read_buffer[MAX_READ_LENGTH];
long index_mem=0;
int is_casava_18;
int is_paired_data;
int is_interleaved;

/*
sdbm
this algorithm was created for sdbm (a public-domain reimplementation of ndbm) database library.
it was found to do well in scrambling bits, causing better distribution of the keys and fewer splits.
it also happens to be a good general hashing function with good distribution.
the actual function is hash(i) = hash(i - 1) * 65599 + str[i]; what is included below is the faster version used in gawk.
[there is even a faster, duff-device version] the magic constant 65599 was picked out of thin air while experimenting with
different constants, and turns out to be a prime. this is one of the algorithms used in berkeley db (see sleepycat) and elsewhere.
*/
static ulong hashit(char *str) {

  ulong hash = 0;
  int c;
  
  while ((c = *str++))
    hash = c + (hash << 6) + (hash << 16) - hash;
  
  return hash;
}

// return 1 if the headers are the same...0 otherwise
static inline int compare_headers(const char *hdr1,const char *hdr2) {

  unsigned int slen=0;
  //char c;
  while ( hdr1[slen]!='\0' && hdr2[slen]!='\0' ) {
    if ( hdr1[slen]!=hdr2[slen] ) {
      return 0;
    }
    slen++;
  }
  if ( hdr1[slen]!='\r' && hdr1[slen]!='\0') {
    return 1;
  }
  if ( hdr2[slen]!='\r' && hdr2[slen]!='\0') {
    return 1;
  }
  return 0;
}

static INDEX_ENTRY* lookup_header(hashtable sn_index,char *hdr) {
  // lookup hdr in sn_index
  ulong key=hashit(hdr);
  //printf("looking for %s: key=%lu\n",hdr,key);
  INDEX_ENTRY* e=(INDEX_ENTRY*)get_object(sn_index,key);
  while (e!=NULL) {      // confirm that hdr are equal
    if ( !strcmp(hdr,e->hdr)) break;
    e=(INDEX_ENTRY*)get_next_object(sn_index,key);
  }
  return e;
}

//long collisions[HASHSIZE+1];
INDEX_ENTRY* new_indexentry(hashtable ht,char*hdr,int len,long start_pos) {
  
  // Memory chunck: |[index_entry]len bytes+1|
  char *mem_block=(char*)malloc(sizeof(INDEX_ENTRY)+len+1);
  if (mem_block==NULL) { return(NULL);}  
  INDEX_ENTRY *e=(INDEX_ENTRY*)&mem_block[0];

  e->hdr=(char*)&mem_block[sizeof(INDEX_ENTRY)];
  e->entry_start=start_pos;
  
  strncpy(e->hdr,hdr,len);
  e->hdr[len]='\0';
  // add to hash table
  ulong key=hashit(e->hdr);
  //collisions[key%HASHSIZE]++;
  if(insere(ht,key,e)<0) {
    fprintf(stderr,"Error adding %s to index\n",hdr);
    return(NULL);
  }
  index_mem+=sizeof(INDEX_ENTRY)+len+1+sizeof(hashnode);
  return(e);
}

void free_indexentry(INDEX_ENTRY *e) {
  free(e);
  // remove entry from hash table
  return;
}

// return 0 on sucess, 1 otherwise
inline int validate_entry(char *hdr,char *hdr2,char *seq,char *qual,unsigned long linenum,const char* filename) {
  
  // Sequence identifier
  if ( hdr[0]!='@' ) {
    fprintf(stderr,"Error in file %s, line %lu: sequence identifier should start with an @ - %s\n",filename,linenum,hdr);
    return 1;
  }  
  if ( hdr[1]=='\0' || hdr[1]=='\n' || hdr[1]=='\r') {
    fprintf(stderr,"Error in file %s, line %lu: sequence identifier should be longer than 1\n",filename,linenum);
    return 1;
  }
  // sequence
  unsigned int slen=0;
  while ( seq[slen]!='\0' && seq[slen]!='\n' && seq[slen]!='\r' ) {
    // check content: ACGT acgt nN 0123....include the .?
    if ( seq[slen]!='A' && seq[slen]!='C' && seq[slen]!='G' && seq[slen]!='T' &&
	 seq[slen]!='a' && seq[slen]!='c' && seq[slen]!='g' && seq[slen]!='t' &&
	 seq[slen]!='0' && seq[slen]!='1' && seq[slen]!='2' && seq[slen]!='3' &&
	 seq[slen]!='n' && seq[slen]!='N' && seq[slen]!='.' ) {
      fprintf(stderr,"Error in file %s, line %lu: invalid character '%c' (hex. code:'%x'), expected ACGTacgt0123nN.\n",filename,linenum+1,seq[slen],seq[slen]);
      return 1;
    }
    slen++;
  }
  // check len
  if (slen < MIN_READ_LENGTH ) {
    fprintf(stderr,"Error in file %s, line %lu: read length too small - %u\n",filename,linenum+1,slen);
    return 1;
  }
  // hdr2=+
  // be tolerant
  //if (hdr2[1]!='\0' && hdr2[1]!='\n' && hdr2[1]!='\r') {
  //  fprintf(stderr,"Error in file %s, line %lu:  header2 wrong. The line should contain only '+' followed by a newline.\n",filename,linenum+2);
  //  return 1;
  //}  
  if (hdr2[0]!='+') {
    fprintf(stderr,"Error in file %s, line %lu:  header2 wrong. The line should contain only '+' followed by a newline or read name (header1).\n",filename,linenum+2);
    return 1;
  }
  // length of hdr2 should be 1 or be the same has the hdr1
  // ignore the + sign
  hdr2=&hdr2[1];
  hdr=&hdr[1];
  if (hdr2[0]!='\0' && hdr2[0]!='\r' ) {
    if ( !compare_headers(&hdr[1],&hdr2[1]) ) {
      fprintf(stderr,"Error in file %s, line %lu:  header2 differs from header1\nheader 1 \"%s\"\nheader 2 \"%s\"\n",filename,linenum,hdr,hdr2);
      return 1;
    }
  }
  // qual length==slen
  unsigned int qlen=0;
  while ( qual[qlen]!='\0' && qual[qlen]!='\n' && qual[qlen]!='\r') {
    qlen++;    
  }  
  if ( qlen!=slen ) {
    fprintf(stderr,"Error in file %s, line %lu: sequence and quality don't have the same length %u!=%u\n",filename,linenum+3,slen,qlen);
    return 1;
  }
  return 0;
}

char* get_readname(char *s,int *len_p,unsigned long cline,const char *filename) {
  int len;
  if (s[0]!='@' ) {
    fprintf(stderr,"Error in file %s, line %lu: wrong header %s\n",filename,cline,s);
    exit(1);
  } 
  s=&s[1]; // ignore/discard @
  if (is_casava_18) {
    len=0;
    while (s[len]!=' ' && s[len]!='\0') ++len;
    s[len]='\0';
    if  ( s[len-2] == '/' ) {
      // discard /[12]
      s[len-2]='\0';
      len=len-2;
    }
  } else {
    // discard last character if PE && not casava 1.8
    len=strlen(s);
    if (is_paired_data) 
      len--;
    s[len-1]='\0'; 
  }
  *len_p=len;
  //fprintf(stderr,"read=%s=\n",s);
  return(s);
}
// TODO: optimize this
void copy_read(long offset,FILE *from,FILE* to) {
  fseek(from,offset,SEEK_SET);
  fputs(READ_LINE(from),to);
  fputs(READ_LINE(from),to);
  fputs(READ_LINE(from),to);
  fputs(READ_LINE(from),to);
}

void index_file(char *filename,hashtable sn_index,long start_offset,long length) {
  FILE *fd1=fopen(filename,"r");  
  if (fd1==NULL) {
    fprintf(stderr,"Error: Unable to open %s\n",filename);
    exit(1);
  }
  // move to the right position
  if(length>0) {
    fprintf(stderr, "Internal error: Not implemented\n");
    exit(2);
  }
  long cline=1;
  // sn_index creation could be done in parallel
  while(!feof(fd1)) {
    long start_pos=ftell(fd1);
    char *hdr=READ_LINE_HDR(fd1);

    if ( hdr==NULL) break;
    int len;
    //fprintf(stderr,"sn_index: =%s=\n",readname);
    // get seq
    //printf("cline=%ld\nLEN=%ld  hdr=%s\n",cline,len,hdr);
    char *seq=READ_LINE_SEQ(fd1);
    char *hdr2=READ_LINE_HDR2(fd1);
    char *qual=READ_LINE_QUAL(fd1);
    char* readname=get_readname(hdr,&len,cline,filename);
    if (seq==NULL || hdr2==NULL || qual==NULL ) {
      fprintf(stderr,"Error in file %s, line %lu: file truncated?\n",filename,cline);
      exit(1);
    }
    if (validate_entry(hdr,hdr2,seq,qual,cline,filename)!=0) {
      exit(1);
    }
    // check for duplicates
    if ( lookup_header(sn_index,readname)!=NULL ) {
      fprintf(stderr,"Error in file %s, line %lu: duplicated sequence %s\n",filename,cline,readname);
      exit(1);
    }
    if ( new_indexentry(sn_index,readname,len,start_pos)==NULL) {
      fprintf(stderr,"Error in file %s, line %lu: malloc failed?",filename,cline);
      exit(1);
    }
    
    PRINT_READS_PROCESSED(cline/4);
    //
    cline+=4;
  }
  
  fclose(fd1);
  return;
}



inline FILE* open_fastq(char* filename) {

  FILE *fd1=fopen(filename,"r");
  if (fd1==NULL) {
    fprintf(stderr,"Error: Unable to open %s\n",filename);
    exit(1);
  }
  return(fd1);
}
// check if the read name format was generated by casava 1.8
int is_casava_1_8(char *f) {
  regex_t regex;
  int reti;
  int is_casava_1_8=0;
  reti = regcomp(&regex,"[A-Z0-9:]* [12]:[YN]:[0-9]*:.*",0);  
  if ( reti ) { 
    fprintf(stderr, "Internal error: Could not compile regex\n"); 
    exit(2); 
  }
  FILE *fd1=open_fastq(f);
  char *hdr=READ_LINE(fd1);
  fclose(fd1);
  /* Execute regular expression */
  //fprintf(stderr,"%s\n",hdr);
  reti = regexec(&regex, hdr, 0, NULL, 0);
  if ( !reti ) {    // match
    is_casava_1_8=1;
  } 
  /* else{
    char msgbuf[100];
    regerror(reti, &regex, msgbuf, sizeof(msgbuf));
    //fprintf(stderr, "Regex match failed: %s\n", msgbuf);
    } */
  regfree(&regex);
  return is_casava_1_8;
}

int validate_interleaved(char *f) {
  unsigned long cline=1;
  FILE *fd1=NULL;
  
  fprintf(stderr,"Paired-end interleaved\n");
  fd1=open_fastq(f);
  while(!feof(fd1)) {
    //long start_pos=ftell(fd1);
    // Read 1
    char *hdr1=READ_LINE_HDR(fd1);
    if ( hdr1==NULL) break;
    int len;
    char *seq1=READ_LINE_SEQ(fd1);
    char *hdr1_2=READ_LINE_HDR2(fd1);
    char *qual1=READ_LINE_QUAL(fd1);
    // Read 2
    char *hdr2=READ_LINE_HDR2_1(fd1);
    char *seq2=READ_LINE_SEQ2(fd1);
    char *hdr2_2=READ_LINE_HDR2_2(fd1);
    char *qual2=READ_LINE_QUAL2(fd1);
    
    if ( seq1==NULL || hdr1_2==NULL || qual1==NULL ||
	 hdr2==NULL || seq2==NULL || hdr2_2==NULL || qual2==NULL ) {
      fprintf(stderr,"Error in file %s, line %lu: file truncated?\n",f,cline);
      return(1);
    }
    if (validate_entry(hdr1,hdr1_2,seq1,qual1,cline,f)!=0) {
      return(1);
    }
    if (validate_entry(hdr2,hdr2_2,seq2,qual2,cline+4,f)!=0) {
      return(1);
    }
    char* readname1=get_readname(hdr1,&len,cline,f);
    char* readname2=get_readname(hdr2,&len,cline+4,f);
    if ( strcmp(readname1,readname2) ) {
      fprintf(stderr,"Error in file %s, line %lu: unpaired read - %s\n",f,cline,readname1);
      return(1);
    } 
    PRINT_READS_PROCESSED(cline/4);
    //
    cline+=8;
    }
  printf("\n");
  fclose(fd1);
  return 0;
}

int main(int argc, char **argv ) {
  //long paired=0;
  is_paired_data=0;
  is_interleaved=0;
  printf("Version iRAP %s\n",VERSION);
  if (argc<2 || argc>3) {
    fprintf(stderr,"Usage: fastq_validator fastq1 [fastq2 file|pe]\n");
    //fprintf(stderr,"%d",argc);
    exit(1);
  }

  FILE *fd1=NULL;
  FILE *fd2=NULL;
  // open & close
  fd1=open_fastq(argv[1]);
  fclose(fd1);
  //fprintf(stderr,"%d\n",argc);
  //bin/fprintf(stderr,"%s\n",argv[0]);
  if (argc ==3) {
    is_paired_data=1;
    if ( !strncmp(argv[2],"pe",2) ) {
      is_interleaved=1;
    } else  {
      fd2=open_fastq(argv[2]);
      fclose(fd2);
    }
  }
  // ************************************************************
  // casava 1.8?
  is_casava_18=is_casava_1_8(argv[1]);
  if (is_casava_18) fprintf(stderr,"CASAVA=1.8\n");
  // ************************************************************
  //off_t cur_offset=1;
  // interleaved
  if ( is_interleaved ) {
    exit(validate_interleaved(argv[1]));
  }
  unsigned long cline=1;
  fprintf(stderr,"HASHSIZE=%lu\n",(long unsigned int)HASHSIZE);
  //memset(&collisions[0],0,HASHSIZE+1);
  hashtable sn_index=new_hashtable(HASHSIZE);
  index_mem+=sizeof(hashtable);

  index_file(argv[1],sn_index,0,-1);
  fprintf(stderr,"\n");
  // print some info
  fprintf(stderr,"Reads processed: %ld\n",sn_index->n_entries);
  fprintf(stderr,"Memory used in indexing: ~%ld MB\n",index_mem/1024/1024);  
  // pair-end
  if (argc ==3 ) {
    fprintf(stderr,"File %s processed\n",argv[1]);  
    fprintf(stderr,"Next file %s\n",argv[2]);  
    // validate the second file and check if all reads are paired
    fd2=open_fastq(argv[2]);
    INDEX_ENTRY* e;
    // read the entry using another fd
    cline=1;
    // TODO: improve code - mostly duplicated:(
    while(!feof(fd2)) {
      //long start_pos=ftell(fd2);
      char *hdr=READ_LINE_HDR(fd2);
      if ( hdr==NULL) break;
      int len;
      char *seq=READ_LINE_SEQ(fd2);
      char *hdr2=READ_LINE_HDR2(fd2);
      char *qual=READ_LINE_QUAL(fd2);
      char* readname=get_readname(hdr,&len,cline,argv[2]);
      if (seq==NULL || hdr2==NULL || qual==NULL ) {
	fprintf(stderr,"Error in file %s, line %lu: file truncated?\n",argv[2],cline);
	exit(1);
      }
      if (validate_entry(hdr,hdr2,seq,qual,cline,argv[2])!=0) {
	exit(1);
      }
      //fprintf(stderr,"Reads processed: %ld\n",sn_index->n_entries);
      // check for duplicates
      if ( (e=lookup_header(sn_index,readname))==NULL ) {
	fprintf(stderr,"Error in file %s, line %lu: unpaired read - %s\n",argv[2],cline,readname);
	exit(1);
      } else {
	ulong key=hashit(readname);
	// remove entry from sn_index
	if (delete(sn_index,key,e)!=e) {
	  fprintf(stderr,"Error in file %s, line %lu: unable to delete entry from sn_index - %s\n",argv[2],cline,readname);
	  exit(1);
	}
	free_indexentry(e);
      }
      PRINT_READS_PROCESSED(cline/4);
      //
      cline+=4;
    }
    printf("\n");
    if (sn_index->n_entries>0 ) {
      fprintf(stderr,"Error in file %s: found %lu unpaired reads\n",argv[1],sn_index->n_entries);
      exit(1);
    }
  }
  printf("OK\n");  
  exit(0);
}
