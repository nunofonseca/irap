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

int main(int argc, char **argv ) {

  if (argc!=2) {
    fprintf(stderr,"Usage: fastq_num_reads fastq_file\n");
    exit(1);
  }

  FASTQ_FILE *fd1=fastq_new(argv[1],FALSE,"r");
  FASTQ_ENTRY *m1=fastq_new_entry();

  while(!gzeof(fd1->fd)) {
    if (fastq_read_next_entry(fd1,m1)==0) break;
  }
  printf("%lu\n",fd1->num_rds);
  fastq_destroy(fd1);
  exit(0);
}

