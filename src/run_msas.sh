#!/bin/bash

for i in `ls *.fasta`; do
    ID=$( cut -d '.' -f 1 <<< $i )
    FID=$ID'.a3m'
    hhblits -i $i -d $HHLIB/db/uniprot20_2016_02/uniprot20_2016_02 -n 3 -cpu 4 -oa3m $FID
done

