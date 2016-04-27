#!/bin/bash

hhblits -i "$1".fasta -d /usr/local/hh-suite/db/uniprot20_2016_02/uniprot20_2016_02 -oa3m "$1".a3m -n $3 -cpu 3 
addss.pl "$1".a3m "$1".a3m
hhsearch -i "$1".a3m -d $2 -cpu 3 $4

