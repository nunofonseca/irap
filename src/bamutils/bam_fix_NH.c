/*
 * =========================================================
 * Copyright 2012-2014,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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

#include "../fastq_utils/hash.h"
#include <string.h>
#include <assert.h>

#define HASHSIZE 7000001
// Same version as IRAP
#define VERSION "0.5.1.d1d"
struct read {
  uint8_t ctr; // how many times a read appears in one alignment
  char *name;// read name
};
typedef struct read READ_ALN;



#define BACKLINE "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
#define PRINT_ALNS_PROCESSED(c) { if (c%100000==0) { fprintf(stderr,"%s%ld",BACKLINE,c);fflush(stderr); }}

unsigned long index_mem=0;
static unsigned long hashit(unsigned char *str);

READ_ALN* get_read_aln(hashtable ht,char*read_name) {
  ulong key=hashit(read_name); 
  READ_ALN* r=(READ_ALN*)get_object(ht,key);
  while (r!=NULL) {      // confirm that hdr are equal
    if (!strcmp(read_name,r->name)) break;
    r=(READ_ALN*)get_next_object(ht,key);
  }
  return r;
}

#define MAX_READ_NAME_LENGTH 1024
char tmp_read_name[MAX_READ_NAME_LENGTH];
char* fix_read_name(char *read_name,int paired) {
  if (paired==2) {
    int len=strlen(read_name);
    assert(len<MAX_READ_NAME_LENGTH-2);
    strncpy(tmp_read_name,read_name,strlen(read_name)+1);
    // add info about the read being read2 of the pair
    tmp_read_name[len]='<';
    tmp_read_name[len+1]='\0';
    return &tmp_read_name[0];
  } 
  return read_name;
}

READ_ALN* new_read_aln(hashtable ht,char*read_name) {

  // look for read in HT
  // add to hash table
  READ_ALN* r;
  if((r=get_read_aln(ht,read_name))==NULL) {
    // Memory chunck: |[index_entry]len bytes+1|
    int len=strlen(read_name);
    ulong key=hashit(read_name);
    char *mem_block=(char*)malloc(sizeof(READ_ALN)+len+1);
    index_mem+=sizeof(READ_ALN)+len+1+sizeof(hashnode);
    if (mem_block==NULL) { 
      fprintf(stderr,"Error allocating memory\n");
      exit(1);
    }
    r=(READ_ALN*)&mem_block[0];    
    r->name=(char*)&mem_block[sizeof(READ_ALN)];
    r->name[len]='\0';
    strncpy(r->name,read_name,len+1);
    r->ctr=0;
    if(insere(ht,key,r)<0){
      fprintf(stderr,"Error adding %s to index\n",read_name);
      exit(1);
    }
  } 
  r->ctr++;
  return(r);
}

void free_read_aln(READ_ALN *r) {
  free(r);
  return;
}


/*
sdbm
this algorithm was created for sdbm (a public-domain reimplementation of ndbm) database library.
it was found to do well in scrambling bits, causing better distribution of the keys and fewer splits.
it also happens to be a good general hashing function with good distribution.
the actual function is hash(i) = hash(i - 1) * 65599 + str[i]; what is included below is the faster version used in gawk.
[there is even a faster, duff-device version] the magic constant 65599 was picked out of thin air while experimenting with
different constants, and turns out to be a prime. this is one of the algorithms used in berkeley db (see sleepycat) and elsewhere.
*/
static unsigned long hashit(unsigned char *str) {

  unsigned long hash = 0;
  int c;
  
  while (c = *str++)
    hash = c + (hash << 6) + (hash << 16) - hash;
  
  return hash;
}


int main(int argc, char *argv[])  
{  
  short out2stdout=0;
  hashtable ht=new_hashtable(HASHSIZE);
  bamFile in,in2; 
  bamFile out; 
  int paired;//1 if not paired or pair read 1, 2 otherwise
  index_mem=sizeof(hashtable)*sizeof(hashnode**)*HASHSIZE*2;

  if (argc != 3) {  
    fprintf(stderr, "Usage: bam_fix_NH <in.bam> <out.bam or - for stdout>\n");  
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
    fprintf(stderr,"bam_fix_NH version %s\n",VERSION);
    fprintf(stderr,"Processing %s\n",argv[1]);
    fprintf(stderr,"Hashing...\n");fflush(stderr);
  }

  while(bam_read1(in,aln)>=0) { // read alignment
    if (aln->core.tid < 0) continue;//ignore unaligned reads
    if (aln->core.flag & BAM_FUNMAP) continue;
    if (aln->core.flag & BAM_FREAD2) paired=2;
    else paired=1;
    ++num_alns;
    new_read_aln(ht,fix_read_name(bam1_qname(aln),paired));
    if(!out2stdout) PRINT_ALNS_PROCESSED(num_alns);
  }
  bam_close(in);  
  if(!out2stdout) {
    fprintf(stderr,"%s%lu\n",BACKLINE,num_alns);
    fprintf(stderr,"Hashing complete (%lu alignments)\n",num_alns);
    fprintf(stderr,"Memory used: %ld MB\n",index_mem/1024/1024);  
    fprintf(stderr,"Updating entries with NH and printing BAM...\n");
    fflush(stderr);
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
    paired=1;
    if (aln->core.tid < 0) continue;//ignore unaligned reads
    if (aln->core.flag & BAM_FUNMAP) continue;
    if (aln->core.flag & BAM_FREAD2) paired=2;
    ++num_alns;
    READ_ALN *r=get_read_aln(ht,fix_read_name(bam1_qname(aln),paired));

    assert(r!=NULL);
    // update the NH field
    uint8_t *old_nh = bam_aux_get(aln, "NH");    
    int32_t nh=r->ctr;
    if (old_nh) {
      if (nh!=bam_aux2i(old_nh)) {
	fprintf(stderr,"warning: value mismatch! replacing>%s %d->%d\n",bam1_qname(aln),bam_aux2i(old_nh),nh);
      }
      bam_aux_del(aln, old_nh);
      bam_aux_append(aln, "NH", 'i', 4, (uint8_t*)&nh);
#ifdef DEBUG
      //      printf("!>%s %d\n",bam1_qname(aln),r->ctr);
#endif
    }
    if (!old_nh) { // add NH  
      bam_aux_append(aln, "NH", 'i', 4, (uint8_t*)&nh);
#ifdef DEBUG
      fprintf(stderr,"!>%s %d\n",bam1_qname(aln),bam_aux2i(old_nh));
#endif
    }
    bam_write1(out,aln);
    if(!out2stdout) PRINT_ALNS_PROCESSED(num_alns);
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
