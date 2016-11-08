#!/usr/bin/env python3

import re
import os

def parse_hhrepid_file(filename, pval = 1e-3):
    basedir = '../data/phylogeny/repeats/'
    seqpat = re.compile(r'(\w+)\s+(sp[\w\|]+)+\s+([\d-]+)\s+\+\d+\s(.+)')
    scorepat = re.compile(r'(\w+)\s+[\d\.]+\s+([\d\.\w-]+)')
    repeat_pvals = {}
    repeat_locs = {}
    with open(basedir + filename) as infile:
        for line in infile:
            scorehit = scorepat.match(line)
            if scorehit:
                repeat_pvals[scorehit.groups()[0]] = float(scorehit.groups()[1])
            seqhit = seqpat.match(line)
            if seqhit == None:
                continue
            repeat = seqhit.groups()[0]
            if repeat_pvals[repeat] > pval:
                continue
            uniprot = seqhit.groups()[1].split('|')[1]
            position = [int(i) for i in seqhit.groups()[2].split('-')]
            repeat_locs[repeat] = (position[0] - 1, position[1])
    return uniprot, repeat_locs

def get_repeats_from_fasta(uniprot, repeats):
    basedir = '../data/phylogeny/repeats/done/'
    with open(basedir + uniprot + '.a3m') as infile:
        label = infile.readline().strip().split()[0]
        line = infile.readline().strip()
        seq = line.strip()
    for rep in repeats:
        start, end = repeats[rep]
        repeat = seq[start:end]
        print(label + '|' + rep)
        print(repeat)

def get_all_repeats():
    pat = re.compile(r'.+hhrepid')
    for f in os.listdir('../data/phylogeny/repeats/'):
        if pat.match(f):
            uniprot, repeats = parse_hhrepid_file(f)
            get_repeats_from_fasta(uniprot, repeats)

if __name__ == '__main__':
    get_all_repeats()
