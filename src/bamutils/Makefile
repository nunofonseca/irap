CFLAGS=-O1
#CFLAGS=-g -DDEBUG

OFILES=bam_fix_NH  bam_pe_insert bam_rm_NH bam_fix_se_flag

all: $(OFILES)

clean: 	
	rm -f *.o $(OFILES) bam_fix_se_flag

bam_pe_insert: bam_pe_insert.c
	gcc  $<   -L $(IRAP_DIR)/lib -lbam -lm -lz -o $@  -I $(IRAP_DIR)/include/bam/ 


bam_fix_se_flag: bam_fix_se_flag.c
	gcc  $(CFLAGS) $^   -L $(IRAP_DIR)/lib -lbam -lm -lz -o $@  -I $(IRAP_DIR)/include/bam/ 

bam_fix_NH: bam_fix_NH.o ../fastq_utils/hash.o
	gcc  $(CFLAGS) $^   -L $(IRAP_DIR)/lib -lbam -lm -lz -o $@  -I $(IRAP_DIR)/include/bam/ 

bam_rm_NH: bam_rm_NH.o 
	gcc  $(CFLAGS) $^   -L $(IRAP_DIR)/lib -lbam -lm -lz -o $@  -I $(IRAP_DIR)/include/bam/ 

bam_fix_NH.o: bam_fix_NH.c
	gcc $(CFLAGS) -c  $<   -I $(IRAP_DIR)/include/bam/ 

bam_rm_NH.o: bam_rm_NH.c
	gcc $(CFLAGS) -c  $<   -I $(IRAP_DIR)/include/bam/ 

../fastq_utils/hash.o: ../fastq_utils/hash.c ../fastq_utils/hash.h
	gcc $(CFLAGS) -c $<

tests: run_bam_fix_tests

run_bam_fix_tests: tests/bam_fix_NH_test1.bam tests/bam_fix_NH_test2.bam tests/bam_fix_NH_test3.bam
	for f in $^; do ./bam_fix_NH $$f /dev/null >/dev/null; done && echo "bam_fix_nh OK"
