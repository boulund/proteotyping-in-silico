#!/usr/bin/env python
# Fredrik Boulund 2015
# Convert output from EMBOSS digest in SRS format to FASTA.

from sys import argv, exit
from textwrap import TextWrapper

if len(argv)<2:
    print "usage: convert_srs_to_fasta.py FILE"
    exit()

with open(argv[1]) as srs:
    tw = TextWrapper()
    for line in srs:
        if line.startswith("# Sequence:"):
            _, _, header, _, _, _, _ = line.split()
        elif line.startswith("Feature:"):
            _, feature = line.split()
            print ">{h}{f}".format(h=header, f=feature)
        elif line.startswith("Sequence:"):
            _, seq = line.split()
            for seqline in tw.wrap(seq):
                print seqline
