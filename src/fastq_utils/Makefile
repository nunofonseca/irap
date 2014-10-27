# avoid having boost in  -I : it may result on including the wrong regexp.h
CFLAGS=-O3 -L$(IRAP_DIR)/lib
#-g

ifdef DEBUG
CFLAGS+= -g -DDEBUG=1
endif

TARGETS=fastq_truncate fastq_filterpair fastq_randsample fastq_trim fastq_filter_n fastq_validator fastq_num_reads fastq_info

############################################################################
all: $(TARGETS)


fastq_validator:  hash.o fastq_validator.o
	gcc  $(CFLAGS) $^ -o $@

fastq_info:  hash.o fastq_info.o
	gcc  $(CFLAGS) $^ -o $@ -lz

fastq_filterpair: hash.o fastq_filterpair.o
	gcc  $(CFLAGS) $^ -o $@

fastq_num_reads: fastq_num_reads.c
	gcc  $(CFLAGS)  -o $@  $<

hash.o: hash.c hash.h
	gcc $(CFLAGS) -c $<

%.o: %.c
	gcc $(CFLAGS) -c $<

fastq_truncate:  fastq_truncate.c
	gcc  $(CFLAGS) -o $@ $<

fastq_trim:  fastq_trim.c
	gcc  $(CFLAGS) -o $@ $<

fastq_filter_n: fastq_filter_n.o
	gcc  $(CFLAGS) -o $@ $<

fastq_randsample:  fastq_randsample.o 
	gcc  $(CFLAGS) -o $@ $^

clean:
	rm -f *.o $(TARGETS) *~

###########################################################################

test1:
	./fastq_filterpair tests/test_2.fastq tests/test_2.fastq  f1.fastq f2.fastq up.fastq
	diff f1.fastq tests/test_2.fastq

test2:
	./fastq_filterpair tests/a_1.fastq tests/a_2.fastq  f1.fastq f2.fastq up.fastq
	diff f2.fastq tests/a_2.fastq
	diff f1.fastq tests/a_1.fastq

test3:
	time -p ./fastq_filterpair tests/casava.1.8_2.fastq tests/casava.1.8_2.fastq  f1.fastq f2.fastq up.fastq
	time -p ./fastq_filterpair tests/casava.1.8_1.fastq tests/casava.1.8_1.fastq  f1.fastq f2.fastq up.fastq
#	diff f2.fastq tests/casava.1.8_2.fastq
#	diff f1.fastq tests/casava.1.8_1.fastq
# the order can be different
#	diff up.fastq test.fq

test4:
	echo "It should not take more than 30 seconds...."
	time -p ./fastq_filterpair tests/c18_1M_1.fastq tests/c18_1M_2.fastq  f1.fastq f2.fastq up.fastq
	time -p ./fastq_filterpair tests/c18_1M_2.fastq tests/c18_1M_2.fastq  f1.fastq f2.fastq up.fastq
	time -p ./fastq_filterpair tests/c18_1M_1.fastq tests/c18_1M_1.fastq  f1.fastq f2.fastq up.fastq

test5:
	echo "It should not take more than 30 seconds...."
	time -p ./fastq_filterpair tests/c18_10000_1.fastq tests/c18_10000_2.fastq  f1.fastq f2.fastq up.fastq

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
	echo "OK! No error detected"

test7:
	@./fastq_info tests/test_e1.fastq || \
	./fastq_info tests/test_e2.fastq || \
	./fastq_info tests/test_e3.fastq || \
	./fastq_info tests/test_e4.fastq || \
	./fastq_info tests/test_e5.fastq || \
	./fastq_info tests/test_e6.fastq || \
	./fastq_info tests/test_e7.fastq || \
	./fastq_info tests/test_e8.fastq || \
	./fastq_info tests/test_e9.fastq || \
	./fastq_info tests/test_e10.fastq || \
	./fastq_info tests/test_e14.fastq || \
	./fastq_info tests/test_e15.fastq || \
	./fastq_info tests/test_e16.fastq || \
	time -p ./fastq_info tests/c18_1M_1.fastq tests/c18_1M_2.fastq || \
	./fastq_info tests/casava.1.8i_err.Fastq || \
	echo "OK! all errors detected"
	./fastq_info -f tests/test_dot.fastq && \
	./fastq_info tests/test_dot.fastq && \
	./fastq_info tests/test_1.fastq && \
	./fastq_info tests/test_2.fastq && \
	./fastq_info tests/test_13.fastq && \
	./fastq_info tests/test_17.fastq && \
	time -p ./fastq_info tests/pe_bug14.fastq tests/pe_bug14.fastq && \
	time -p ./fastq_info tests/c18_1M_1.fastq && \
	time -p ./fastq_info tests/c18_1M_2.fastq && \
	time -p ./fastq_info tests/c18_1M_1.fastq.gz && \
	time -p ./fastq_info tests/c18_1M_2.fastq.gz && \
	time -p ./fastq_info tests/c18_1M_1.fastq.gz tests/c18_1M_1.fastq.gz && \
	time -p ./fastq_info tests/casava.1.8i.fastq pe && \
	time -p ./fastq_info tests/solexa_1.fastq tests/solexa_2.fastq && \
	 echo "OK! No error detected"


tests: test1 test2 test3 test4 test5 test6 test7
