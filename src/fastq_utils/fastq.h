

#define DEFAULT  0
#define CASAVA18 1
#define INTEGERNAME 2

#ifndef MAX_READ_LENGTH
// 5Mb should cover most of the cases for now :)
#define MAX_READ_LENGTH 5000000
#endif

#ifndef MAX_LABEL_LENGTH
#define MAX_LABEL_LENGTH 10000
#endif

#ifndef MAX_FILENAME_LENGTH
#define MAX_FILENAME_LENGTH 5000
#endif


#define MIN_READ_LENGTH 1
#define TRUE 1
#define FALSE 0
#define UNDEF -1

#ifndef HASHSIZE
#define HASHSIZE 19000001
#endif

#include "hash.h"
#include <zlib.h> 


#define min(a,b) (a<b?a:b)
#define max(a,b) (a>b?a:b)
#define PRINT_READS_PROCESSED(c,n) { if (c%n==0) { fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%lu",c);fflush(stderr); }}

extern unsigned long index_mem;
extern char* encodings[];

struct index_entry {  
  // file offset: start entry
  // file offset: end entry
  // chat hdr(40)
  char *hdr;
  off_t entry_start;
  //unsigned  int  nbytes;
};
typedef struct index_entry INDEX_ENTRY;

struct fastq_entry {  
  // file offset: start entry
  // file offset: end entry
  // chat hdr(40)
  char hdr1[MAX_LABEL_LENGTH];
  char hdr2[MAX_LABEL_LENGTH];
  char seq[MAX_READ_LENGTH];
  char qual[MAX_READ_LENGTH];
  unsigned long read_len;
  long long offset;
};
typedef struct fastq_entry FASTQ_ENTRY;

struct fastq_file {
  gzFile fd;
  long long cur_offset;
  unsigned long cline;
  char filename[MAX_FILENAME_LENGTH];

  unsigned long max_rl; // maximum read length
  unsigned long last_rl; // read length of the last read
  unsigned long min_rl;  // minimum read length
  unsigned long min_qual; // minimum quality
  unsigned long max_qual;   // maximum quality
  unsigned long num_rds; // number of reads
  unsigned long rdlen_ctr[MAX_READ_LENGTH]; // keep a tally on how many reads we observed per length
  
  int fix_dot;
  int fixed_dot;
  int is_pe;
  int readname_format;
  int is_casava_18;
};
typedef struct fastq_file  FASTQ_FILE;

FASTQ_ENTRY* fastq_new_entry(void);
void fastq_write_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e);

void fastq_index_delete(char *rname,hashtable index);
INDEX_ENTRY* fastq_index_lookup_header(hashtable sn_index,char *hdr);
char* fastq_get_readname(FASTQ_FILE*, FASTQ_ENTRY *,char* rn,unsigned long*,int is_header1);
int fastq_read_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e);
void fastq_new_entry_stats(FASTQ_FILE *, FASTQ_ENTRY* );
int fastq_validate_entry(FASTQ_FILE *fd,FASTQ_ENTRY *e);
int fastq_read_next_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e);

FASTQ_FILE* fastq_new(const char* filename, const int fix_dot,const char *mode);
void fastq_destroy(FASTQ_FILE*);
void fastq_is_pe(FASTQ_FILE* fd);
void fastq_index_readnames(FASTQ_FILE *,hashtable,long long,int);
void fastq_write_entry(FASTQ_FILE* fd,FASTQ_ENTRY *e);
void fastq_write_entry2stdout(FASTQ_ENTRY *e);
void fastq_seek_copy_read(long offset,FASTQ_FILE* from,FASTQ_FILE *to);
char* fastq_qualRange2enc(int min_qual,int max_qual);
void fastq_rewind(FASTQ_FILE* fd);
