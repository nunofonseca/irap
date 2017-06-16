#!/bin/bash

function must_fail {
    cmd=$*
    STATUS=PASSED
    bash -c "$cmd" 2> /dev/null
    if [ 0 -eq $? ];  then	
	STATUS=FAILED
    fi
    echo $STATUS $cmd
}

function must_succeed {
    cmd=$*
    STATUS=PASSED
    bash -c "$cmd" 2> /dev/null
    if [ 0 -ne $? ];  then	
	STATUS=FAILED
    fi
    echo $STATUS $cmd
}
#############################################
##
echo "*** fastq_filter_n"
must_succeed "./fastq_filter_n tests/test_21_2.fastq > tmp && diff /dev/null tmp"
must_fail "./fastq_filter_n -n 100 tests/test_21_2.fastq > tmp && diff /dev/null tmp"
must_fail "./fastq_filter_n tests/test_21_2.fastq > tmp && diff  tests/test_21_2.fastq tmp"
must_succeed "./fastq_filter_n tests/test_1.fastq > tmp && diff tests/test_1.fastq tmp"
##
echo "*** fastq_num_reads"
must_succeed "[ `./fastq_num_reads tests/test_21_2.fastq` -eq 2 ]"
must_fail "[ `./fastq_num_reads tests/test_21_2.fastq` -ne 2 ]"
must_succeed "[ `./fastq_num_reads tests/c18_10000_1.fastq` -eq 10000 ]"
must_fail "[ `./fastq_num_reads tests/c18_10000_1.fastq` -ne 10000 ]"
must_succeed "[ `./fastq_num_reads tests/one.fastq` -eq 1 ]"
##
echo "*** fastq_truncate"
must_succeed "[ `./fastq_truncate tests/test_21_2.fastq 1|wc -l` -eq 4 ]"
must_succeed "[ `./fastq_truncate tests/test_21_2.fastq 0|wc -l` -eq 0 ]"
must_succeed "[ `./fastq_truncate tests/test_21_2.fastq 2|wc -l` -eq 8 ]"
must_fail "./fastq_truncate tests/test_21_2.fastq"


echo "*** fastq_info"
must_fail ./fastq_info tests/test_e1.fastq 
must_fail ./fastq_info tests/test_e2.fastq
must_fail ./fastq_info tests/test_e3.fastq 
must_fail ./fastq_info tests/test_e4.fastq 
must_fail ./fastq_info tests/test_e5.fastq 
must_fail ./fastq_info tests/test_e6.fastq 
must_fail ./fastq_info tests/test_e7.fastq 
must_fail ./fastq_info tests/test_e8.fastq 
must_fail ./fastq_info tests/test_e9.fastq 
must_fail ./fastq_info tests/test_e10.fastq 
must_fail ./fastq_info tests/test_e14.fastq 
must_fail ./fastq_info tests/test_e15.fastq 
must_fail ./fastq_info tests/test_e16.fastq 
must_fail ./fastq_info tests/test_e17.fastq 
must_fail ./fastq_info -f tests/test_dot.fastq 
must_fail ./fastq_info tests/c18_1M_1.fastq tests/c18_1M_2.fastq 
##
## just checks the exit status
must_succeed ./fastq_info tests/test_dot.fastq 
must_succeed 	./fastq_info tests/test_1.fastq 
must_succeed 	./fastq_info tests/test_2.fastq 
must_succeed 	./fastq_info tests/test_13.fastq
must_succeed 	./fastq_info tests/test_17.fastq
must_succeed 	./fastq_info tests/test_21_1.fastq
must_succeed 	./fastq_info tests/test_21_1.fastq tests/test_21_2.fastq 
must_succeed 	time -p ./fastq_info tests/pe_bug14.fastq tests/pe_bug14.fastq 
must_succeed 	time -p ./fastq_info tests/c18_1M_1.fastq 
must_succeed 	time -p ./fastq_info tests/c18_1M_2.fastq 
must_succeed 	time -p ./fastq_info tests/c18_1M_1.fastq 
must_succeed 	time -p ./fastq_info tests/c18_1M_2.fastq 
must_succeed 	time -p ./fastq_info tests/c18_1M_1.fastq tests/c18_1M_1.fastq 
must_succeed 	time -p ./fastq_info tests/casava.1.8i.fastq pe 
must_succeed 	time -p ./fastq_info tests/solexa_1.fastq tests/solexa_2.fastq 
##
echo "*** fastq_filterpair"
must_succeed "./fastq_filterpair tests/test_2.fastq tests/test_2.fastq  f1.fastq.gz f2.fastq.gz up.fastq.gz && diff <(zcat f1.fastq.gz) tests/test_2.fastq"
must_succeed "./fastq_filterpair tests/a_1.fastq tests/a_2.fastq  f1.fastq.gz f2.fastq.gz up.fastq.gz &&	diff <(zcat f2.fastq.gz) tests/a_2.fastq && diff <(zcat f1.fastq.gz) tests/a_1.fastq"
must_succeed "./fastq_filterpair tests/casava.1.8_2.fastq tests/casava.1.8_2.fastq  f1.fastq.gz f2.fastq.gz up.fastq.gz"
must_succeed "./fastq_filterpair tests/casava.1.8_1.fastq tests/casava.1.8_1.fastq  f1.fastq.gz f2.fastq.gz up.fastq.gz"
must_succeed ./fastq_filterpair tests/c18_1M_1.fastq tests/c18_1M_2.fastq  f1.fastq f2.fastq up.fastq
must_succeed ./fastq_filterpair tests/c18_1M_2.fastq tests/c18_1M_2.fastq  f1.fastq f2.fastq up.fastq 
must_succeed ./fastq_filterpair tests/c18_1M_1.fastq tests/c18_1M_1.fastq  f1.fastq f2.fastq up.fastq 
must_succeed ./fastq_filterpair tests/c18_10000_1.fastq tests/c18_10000_2.fastq  f1.fastq.gz f2.fastq.gz up.fastq.gz


exit 0


