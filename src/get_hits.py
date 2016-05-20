#!/usr/bin/env python3

import sys

"Returns list of uniprot ids from fasta labels"
def get_hits(filename):
    with open(filename) as infile:
        data = infile.read().split('\n\n')
    hits = data[1].split('\n')[1:]
    uprids = [hit.split('|')[1] for hit in hits]
    return uprids

def write_download_list(hits):
    baseurl = 'http://www.uniprot.org/uniprot/'
    suffix = '.fasta'
    with open('tmp_download.dat', 'w') as outfile:
        for uprid in hits:
            line = baseurl + uprid + suffix + '\n'
            outfile.write(line)

if __name__ == '__main__':
    hits = get_hits(sys.argv[1])
    write_download_list(hits)

