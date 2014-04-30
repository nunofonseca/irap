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
 *    $Id$
 * =========================================================
 */
#include <stdio.h>  
#include <math.h>
#include <bam.h>
#include <sam.h>
#include <kstring.h>      

#include "../fastq_utils/hash.h"
#include <string.h>
//#include <assert.h>

#define HASHSIZE 2000000

// keep a record of how many times a read appears in one alignment
struct read {
  short ctr;
  char *name;
};
typedef struct read READ_ALN;


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
READ_ALN* new_read_aln(hashtable ht,char*read_name) {

  // look for read in HT
  // add to hash table
  READ_ALN* r;
  if((r=get_read_aln(ht,read_name))==NULL) {
    // Memory chunck: |[index_entry]len bytes+1|
    int len=strlen(read_name);
    ulong key=hashit(read_name);
    char *mem_block=(char*)malloc(sizeof(READ_ALN)+len+1);
    index_mem+=sizeof(READ_ALN)+len+1;
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
    //printf("New %s\n",r->name);
  } 
  //else {
  //  printf("HIT %d\n",r->ctr);
  //}
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
  hashtable ht=new_hashtable(HASHSIZE);
  bamFile in,in2; 
  bamFile out; 
  
  if (argc != 3) {  
    fprintf(stderr, "Usage: bam_fix_NH <in.bam> <out.bam>\n");  
    return 1;  
  }  
  
  // Open file and exit if error
  //in = strcmp(argv[1], "-")? bam_open(argv[1], "rb") : bam_dopen(fileno(stdin), "rb");
  in = bam_open(argv[1], "rb");
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

  printf("Hashing...\n");flush(stdout);
  while(bam_read1(in,aln)>=0) { // read alignment
    if (aln->core.tid < 0) continue;//ignore unaligned reads
    ++num_alns;
    new_read_aln(ht,bam1_qname(aln));
  }
  bam_close(in);  
  printf("Hashing complete (%lu alignments)\n",num_alns);
  printf("Memory used in the hash: %ld MB\n",index_mem/1024/1024);  
  flush(stdout);
  // reopen
  in2 = bam_open(argv[1], "rb");
  if (in2 == 0 ) {  
    fprintf(stderr, "ERROR: Fail to open BAM file %s\n", argv[1]);  
    return 1;  
  }  

  header = bam_header_read(in2);
  
  while(bam_read1(in2,aln)>=0) { // read alignment
    if (aln->core.tid < 0) continue;//ignore unaligned reads
    ++num_alns;
    
    READ_ALN *r=get_read_aln(ht,bam1_qname(aln));

    //assert(r!=NULL);
    // update the NH field
    uint8_t *old_nh = bam_aux_get(aln, "NH");    
    uint8_t nh=r->ctr;
    if (old_nh) {
      if (nh!=bam_aux2i(old_nh)) {
	fprintf(stderr,"warning: value mismatch! replacing>%s %d->%d\n",bam1_qname(aln),bam_aux2i(old_nh),nh);
      }
      bam_aux_del(aln, old_nh);
      bam_aux_append(aln, "NH", 'i', 4, (uint8_t*)&nh);
    }
    if (!old_nh) { // add NH  
      bam_aux_append(aln, "NH", 'i', 4, (uint8_t*)&nh);
#ifdef DEBUG
      printf("!>%s %d\n",bam1_qname(aln),bam_aux2i(old_nh));
#endif
    }
    // in->header
    // Also fix the XS:A tag
    // BAM_FREAD1
    // BAM_FREAD2
    // BAM_FREVERSE the read is mapped to the reverse strand 
    //bam1_cigar(b) 
      //BAM_CREF_SKIP 3 CIGAR skip on the reference (e.g. spliced alignment)
      //BAM_FREVERSE 16 the read is mapped to the reverse strand
    if (aln->core.flag & BAM_FSECONDARY) continue; // skip secondary alignments
    if (aln->core.flag & ! BAM_FPAIRED) continue; // not paired
    if (aln->core.flag & ! BAM_FPROPER_PAIR) continue; // not a proper pair
    if (aln->core.flag & ! BAM_FMUNMAP) continue; // the mate is mapped
    if (aln->core.flag & BAM_FSECONDARY) continue; // secundary read
    if (aln->core.flag & BAM_FREAD2) continue; // only count each pair once
    // core.strand == 0 (f/+) 1 r/-
    // flag
    // bam1_qname(b)
    bam_write1(out,aln);
  }
  // 
  bam_destroy1(aln);
  bam_close(in2);  
  bam_close(out);  
  return 0;  
/*
uint8_t *old_nm = bam_aux_get(b, "NM");
90 	if (c->flag & BAM_FUNMAP) return;
91 	if (old_nm) old_nm_i = bam_aux2i(old_nm);
92 	if (!old_nm) bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm);
93 	else if (nm != old_nm_i) {
94 	fprintf(stderr, "[bam_fillmd1] different NM for read '%s': %d -> %d\n", bam1_qname(b), old_nm_i, nm);
95 	bam_aux_del(b, old_nm);
96 	bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm);
97 	}
*/
}  
