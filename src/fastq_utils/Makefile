CFLAGS=-O2
#-g

ifdef DEBUG
CFLAGS+= -g -DDEBUG=1
endif

TARGETS=fastq_truncate fastq_filterpair fastq_randsample fastq_trim fastq_filter_n 

############################################################################
all: $(TARGETS)

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




tests: test1 test2 test3 test4 test5
