CFLAGS=-O1
#CFLAGS=-g -DDEBUG
IDIRS=$(IRAP_DIR)/bin/samtools0.x/include/
LIBDIRS=$(IRAP_DIR)/bin/samtools0.x/lib

OFILES=bam_fix_NH  bam_pe_insert bam_rm_NH bam_fix_se_flag bam_tophat2_fix bamRindex

all: $(OFILES)

clean: 	
	rm -f *.o $(OFILES) bam_fix_se_flag

bam_pe_insert: bam_pe_insert.c
	gcc  $<   -L $(LIBDIRS) -lpthread -lbam -lm -lz -pthread -o $@  -I $(IDIRS) 

bam_tophat2_fix: bam_tophat2_fix.c
	gcc  $<   -L $(LIBDIRS) -lpthread -lbam -lm -lz  -pthread -o $@  -I $(IDIRS) 

bamStrand: bamStrand.c
	gcc  $<   -L $(LIBDIRS) -lpthread -lbam -lm -lz  -pthread -o $@  -I $(IDIRS) 


# -lrt
bamRindex: bamRindex.c
	gcc  $<   -L $(LIBDIRS)  -lpthread  -lbam -lm -lz -lsqlite3  -pthread  -o $@  -I $(IDIRS) 


bam_fix_se_flag: bam_fix_se_flag.c
	gcc  $(CFLAGS) $^   -L $(LIBDIRS)  -lpthread  -lbam -lm -lz  -pthread -o $@  -I $(IDIRS) 

bam_fix_NH: bam_fix_NH.o hash.o
	gcc  $(CFLAGS) $^   -L $(LIBDIRS)  -lpthread  -lbam -lm -lz  -pthread -o $@  -I $(IDIRS) 

bam_rm_NH: bam_rm_NH.o 
	gcc  $(CFLAGS) $^   -L $(LIBDIRS)  -lpthread -lbam -lm -lz  -pthread -o $@  -I $(IDIRS) 

bam_fix_NH.o: bam_fix_NH.c
	gcc $(CFLAGS) -c  $<   -I $(IDIRS) 

bam_rm_NH.o: bam_rm_NH.c
	gcc $(CFLAGS) -c  $<   -I $(IDIRS) 

hash.o: hash.c hash.h
	gcc $(CFLAGS) -c $<

tests: run_bam_fix_tests run_bamRindex_tests

run_bam_fix_tests: tests/bam_fix_NH_test1.bam tests/bam_fix_NH_test2.bam tests/bam_fix_NH_test3.bam
	for f in $^; do ./bam_fix_NH $$f /dev/null >/dev/null; done && echo "bam_fix_nh OK"

# Initial tests
# TODO: improve this
run_bamRindex_tests: tests/bam_fix_NH_test1.bam tests/bam_fix_NH_test2.bam tests/bam_fix_NH_test3.bam
	for f in $^; do ./bamRindex $$f; done && echo "bamRindex OK"
	for f in $^; do ./bamRindex $$f test; done && echo "bamRindex OK"
