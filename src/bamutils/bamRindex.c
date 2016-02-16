/*
 * =========================================================
 * Copyright 2012-2016,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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
#include "sqlite3.h"
#include <time.h>
#include <assert.h>

#define BOOLEAN int
#define FALSE 0
#define TRUE 1

#define BACKLINE "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
#define PRINT_ALNS_PROCESSED(c) { if (c%100000==0) { fprintf(stderr,"%s%ld",BACKLINE,c);fflush(stderr); }}


#define TABLE "CREATE TABLE bam_index(isDuplicate Boolean, isNotPassingQualityControls Boolean,  nSpliced smallint,  isPaired Boolean,  isPrimary Boolean,  isMapped Boolean,  hasSimpleCigar Boolean,  isSecondMateRead Boolean,  isProperPair Boolean,  nh smallint,  nm smallint,  mapq smallint, xs character(1));"

#define INDEXES "create index idx1 on bam_index (isMapped); CREATE INDEX idx2 on bam_index ( isPrimary); create index idx3 on bam_index (isProperPair ); create index idx4 on bam_index ( nSpliced ); create index idx5 on bam_index ( nh ); create index idx6 on bam_index ( nm ); create index idx7 on bam_index ( xs );"

#define BUFFER_SIZE 1024

#define SQLITE_CHECK_ERROR() {if (sErrMsg!=NULL) {fprintf(stderr,"ERROR: %s\n",sErrMsg);  exit(1);}}

int main(int argc, char *argv[])  
{  
  bamFile in; 
  sqlite3 * db;
  sqlite3_stmt * stmt;
  char * sErrMsg = NULL;
  char * tail = 0;
  int nRetCode;
  char sSQL [BUFFER_SIZE] = "\0";
  char database[BUFFER_SIZE];
  clock_t startClock,startClock2;

  if ( argc != 2 && argc!=3 ) {  
    fprintf(stderr, "Usage: bamRindex <in.bam> [out.index]\n");  
    return 1;  
  }  
  
  if ( argc == 2 ) {  
    assert(strcpy(database,argv[1])!=NULL);
  } else {
    assert(strcpy(database,argv[2])!=NULL);
  }
  // always add the extension
  assert(strcat(database,".ridx")!=NULL);
  remove(database);
  // Open file and exit if error
  //in = strcmp(argv[1], "-")? bam_open(argv[1], "rb") : bam_dopen(fileno(stdin), "rb");
  //fprintf(stderr,"Options ok\n");
  in = bam_open(argv[1], "rb");
  if (in == 0 ) {  
    fprintf(stderr, "ERROR: Fail to open BAM file %s\n", argv[1]);  
    return 1;  
  }  
  //fprintf(stderr,"BAM opened\n");

  // ***********
  // Read header
  bam_header_t *header;
  header = bam_header_read(in);
  // sorted by name?
  // Should not rely on the value in SO 
  bam1_t *aln=bam_init1();
  unsigned long num_alns=0;

  /*********************************************/
  /* Open the Database and create the Schema */
  // TODO: check the errors
  sqlite3_open(database, &db);
  sqlite3_exec(db, TABLE, NULL, NULL, &sErrMsg); // create the table
  SQLITE_CHECK_ERROR();
  startClock = clock();
  sqlite3_exec(db, "PRAGMA synchronous = 0;", NULL, NULL, &sErrMsg);
  SQLITE_CHECK_ERROR();
  sqlite3_exec(db, "PRAGMA journal_mode = OFF;", NULL, NULL, &sErrMsg);
  SQLITE_CHECK_ERROR();
  // Use up to 8GB of memory
  sqlite3_exec(db, "PRAGMA cache_size = -8000000;", NULL, NULL, &sErrMsg);
  SQLITE_CHECK_ERROR();
  sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, &sErrMsg);
  SQLITE_CHECK_ERROR();
  while(bam_read1(in,aln)>=0) { // read alignment
    //aln->core.tid < 0 ? 
    uint8_t *nh = bam_aux_get(aln, "NH");
    uint8_t *nm = bam_aux_get(aln, "NM");
    uint8_t *xs = bam_aux_get(aln, "XS");
    
    BOOLEAN isPrimary;
    BOOLEAN isMapped;
    BOOLEAN notMapped;
    BOOLEAN isDuplicate;
    BOOLEAN isNotPassingQualityControls;
    BOOLEAN isPaired;
    BOOLEAN isSecondMateRead,isProperPair;
    //secondary alignment
    notMapped=(aln->core.flag & BAM_FUNMAP) ? TRUE: FALSE;
    //notMapped=((aln->core.flag & BAM_FUNMAP) || (aln->core.mtid ==0)) ? TRUE: FALSE;
    isMapped=!notMapped;
    isPrimary= (aln->core.flag & BAM_FSECONDARY) ? FALSE:TRUE;
    isProperPair=(aln->core.flag & BAM_FPROPER_PAIR) ? TRUE:FALSE;
    isPaired=(aln->core.flag & BAM_FPAIRED ) ? TRUE:FALSE;
    isSecondMateRead=(aln->core.flag & BAM_FREAD2 ) ? TRUE: FALSE;
    isNotPassingQualityControls=(aln->core.flag & BAM_FQCFAIL ) ? TRUE:FALSE;
    isDuplicate=(aln->core.flag & BAM_FDUP) ? TRUE: FALSE;

    BOOLEAN isSpliced=FALSE;
    BOOLEAN hasSimpleCigar=TRUE;
    int nSpliced=0;
    int i;
    if (aln->core.n_cigar != 0) {
      for (i = 0; i < aln->core.n_cigar; ++i) {
	char l="MIDNSHP=X"[bam1_cigar(aln)[i]&BAM_CIGAR_MASK];
	//fprintf(stderr,"%c",l);
	if ( l == 'N' ) { isSpliced=TRUE; hasSimpleCigar=FALSE;++nSpliced;}	  
	if ( l != 'M' && l!='=' ) {  hasSimpleCigar=FALSE;}	  
      }
    } 
    //fprintf(stderr,"read %ld\n",num_alns);
    // isDuplicate,isNotPassingQualityControls,
    // isSpliced,isPAired,isPrimary,hasSimpleCigar,isSecondMateRead,isProperPair,nh,nm,qual/mapq,xs
    sprintf(sSQL,"INSERT into bam_index values (%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,'%c')",
	   isDuplicate,isNotPassingQualityControls,
	   nSpliced,isPaired,isPrimary,isMapped,hasSimpleCigar,isSecondMateRead,isProperPair,
	   (nh==0?0:bam_aux2i(nh)),(nm==0?0:bam_aux2i(nm)),
	    aln->core.qual,
	    (xs==0?' ':(bam_aux2A(xs)==0?' ':bam_aux2A(xs))));
    sqlite3_exec(db, sSQL, NULL, NULL, &sErrMsg);
    SQLITE_CHECK_ERROR();
    ++num_alns;
    PRINT_ALNS_PROCESSED(num_alns);
  }
  bam_close(in);  
  sqlite3_exec(db, "END TRANSACTION;", NULL, NULL, &sErrMsg);
  SQLITE_CHECK_ERROR();
  printf("\nImported %d records in %4.2f seconds\n", num_alns, ( (double) (clock() - startClock))/CLOCKS_PER_SEC);
  // Create the indexes
  startClock2 = clock();
  // generating the indexes does not pay off
  //sqlite3_exec(db, INDEXES, NULL, NULL, &sErrMsg);
  //printf("Indexed %d records in %4.2f seconds\n", num_alns, ( (double) (clock() - startClock2))/CLOCKS_PER_SEC);
  printf("Total time: %4.2f seconds\n", ((double)(clock() - startClock))/CLOCKS_PER_SEC);
  sqlite3_close(db);
  return 0;  
}  
