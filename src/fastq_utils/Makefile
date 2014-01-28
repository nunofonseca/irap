CFLAGS=-O3
#-g

ifdef DEBUG
CFLAGS+= -g -DDEBUG=1
endif

TARGETS=fastq_truncate fastq_filterpair fastq_randsample fastq_trim fastq_filter_n fastq_validator

############################################################################
all: $(TARGETS)


fastq_validator:  hash.o fastq_validator.o
	gcc  $(CFLAGS) $^ -o $@

fastq_filterpair: hash.o fastq_filterpair.o
	gcc  $^ -o $@


hash.o: hash.c hash.h
	gcc $(CFLAGS) -c $<

%.o: %.c
	gcc $(CFLAGS) -c $<

fastq_qual_enc:  fastq_qual_enc.c
	gcc  $(CFLAGS) -o $@ $<

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
	echo "OK! all errors detected"
	./fastq_validator tests/test_1.fastq && \
	./fastq_validator tests/test_2.fastq && \
	./fastq_validator tests/test_13.fastq && \
	time -p ./fastq_validator tests/pe_bug14.fastq tests/pe_bug14.fastq && \
	time -p ./fastq_validator tests/c18_1M_1.fastq && \
	time -p ./fastq_validator tests/c18_1M_2.fastq
	time -p ./fastq_validator tests/c18_1M_1.fastq tests/c18_1M_2.fastq || echo "OK! Error detected"


tests: test1 test2 test3 test4 test5 test6
