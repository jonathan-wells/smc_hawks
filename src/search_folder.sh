#!/bin/bash

for i in `ls *.fasta`; do
    base=$(echo $i | awk -F'.' '{print $1}')
    hhblits -i "$base".fasta -d /usr/local/hh-suite/db/uniprot20_2016_02/uniprot20_2016_02 -oa3m "$base".a3m -n $1 -cpu 3 
    addss.pl "$base".a3m "$base".a3m
    hhsearch -i "$base".a3m -d $2 -cpu 3 $3
done

