#!/bin/bash
set -e
GTF_FILE=$1
if [ "$GTF_FILE-" = "-" ]; then
    echo "Error! Usage: gtf2geneclass.sh gtf_file" >&2
    exit 1
fi
echo class gene_id 	
cut -f 2,9  $GTF_FILE | cut -f 1,3 -d\ | tr -d "\";\t"  | sort   -u 
exit 0
