# avoid having boost in  -I : it may result on including the wrong regexp.h
CFLAGS= -L$(IRAP_DIR)/lib -Wall -O3
#-g -Wunused-result 
#
SHELL=bash

ifdef DEBUG
CFLAGS+= -g -DDEBUG=1
endif

TARGETS=fastq_truncate fastq_filterpair fastq_randsample fastq_filter_n fastq_validator fastq_num_reads fastq_info

############################################################################
all: $(TARGETS)

fastq_filterpair: hash.o fastq_filterpair.o fastq.o 
	gcc  $(CFLAGS) $^ -o $@ -lz

fastq_info:  hash.o fastq_info.o fastq.o
	gcc  $(CFLAGS) $^ -o $@ -lz

fastq_filter_n: fastq_filter_n.o fastq.o hash.o
	gcc  $(CFLAGS) $^ -o $@  -lz

fastq_num_reads: fastq_num_reads.o hash.o fastq.o
	gcc  $(CFLAGS) $^  -o $@  -lz

fastq_truncate:  fastq_truncate.o  hash.o fastq.o
	gcc  $(CFLAGS) $^  -o $@  -lz


# deprecated
fastq_validator:  hash.o fastq_validator.o
	gcc  $(CFLAGS) $^ -o $@


fastq.o: fastq.c fastq.h hash.h
	gcc $(CFLAGS) -c $< -lz

hash.o: hash.c hash.h
	gcc $(CFLAGS) -c $<

%.o: %.c
	gcc $(CFLAGS) -c $<


fastq_randsample:  fastq_randsample.o 
	gcc  $(CFLAGS) -o $@ $^

clean:
	rm -f *.o $(TARGETS) *~

###########################################################################

tests:
	./run_tests.sh

test6:
	@./fastq_validator tests/test_e1.fastq || \
	./fastq_validator tests/test_e2.fastq || \
	./fastq_validator tests/test_e3.fastq || \
	./fastq_validator tests/test_e4.fastq || \
	./fastq_validator tests/test_e5.fastq || \
	./fastq_validator tests/test_e6.fastq || \
	./fastq_validator tests/test_e7.fastq || \
	./fastq_validator tests/test_e8.fastq || \
	./fastq_validator tests/test_e9.fastq || \
	./fastq_validator tests/test_e10.fastq || \
	./fastq_validator tests/test_e14.fastq || \
	./fastq_validator tests/test_e15.fastq || \
	./fastq_validator tests/test_e16.fastq || \
	./fastq_validator tests/casava.1.8i_err.Fastq || \
	time -p ./fastq_validator tests/c18_1M_1.fastq tests/c18_1M_2.fastq || \
	echo "OK! all errors detected"
	./fastq_validator tests/test_1.fastq && \
	./fastq_validator tests/test_2.fastq && \
	./fastq_validator tests/test_13.fastq && \
	./fastq_validator tests/test_17.fastq && \
	time -p ./fastq_validator tests/pe_bug14.fastq tests/pe_bug14.fastq && \
	time -p ./fastq_validator tests/c18_1M_1.fastq && \
	time -p ./fastq_validator tests/c18_1M_2.fastq && \
	time -p ./fastq_validator tests/casava.1.8i.fastq pe && \
	time -p ./fastq_validator tests/solexa_1.fastq tests/solexa_2.fastq && \
	echo "OK! No errors detected"


