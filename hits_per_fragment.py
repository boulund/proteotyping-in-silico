#!/usr/bin/env python
# Fredrik Boulund 2015
# Compute the number of BLAT hits for each fragment in the source file

from sys import argv, exit
from read_fasta import read_fasta

def parse_blat_output(filename):
    """Parses blat output in blast8 format. 
    Returns all hits for each fragment.
    """

    hit_counter = 0
    hits = {}

    with open(filename) as blast8:
        for line in blast8:
            blast8_line = line.split()
            fragment_id = blast8_line[0]
            target_id = blast8_line[1]
            identity = float(blast8_line[2])
            matches = int(blast8_line[3])
            mismatches = int(blast8_line[4])
            startpos = int(blast8_line[8])
            endpos = int(blast8_line[9])
            #fragment_length = int(fragment_id.split("_")[1])
            #fragment_coverage = matches/fragment_length
            hit = (target_id, identity, matches, mismatches, 
                    0, startpos, endpos )
            try:
                hits[fragment_id].append(hit)
            except KeyError:
                hits[fragment_id] = [hit]
            hit_counter += 1
    print "Parsed {} hits for {} fragments from {}.".format(hit_counter, len(hits), filename)
    return hits


def count_hits_per_sequence(fasta, blast8):
    """Compare source FASTA with blast8 output to find sequences without hits.
    """
    seqs = read_fasta(fasta)
    seqs = [(h, s) for h,s in seqs]
    maxlen = max([len(s) for h,s in seqs])

    hits = parse_blat_output(blast8)

    hits_per_length = {}
    for n in xrange(6,maxlen+1):
        hits_per_length[n] = [0,0]
    nohits = []

    for header, seq in seqs:
        if header in hits.keys():
            #print "{} hits for {}".format(len(hits[header]), header)
            hits_per_length[len(seq)][0] += 1
        else:
            #print "! NO HITS FOR {}".format(header)
            hits_per_length[len(seq)][1] += 1
            nohits.append((header,seq))

    return nohits, hits_per_length


if __name__ == "__main__":
    if len(argv) <2:
        print "usage: script.py FASTA BLAST8"
        exit()

    nohits, hpl = count_hits_per_sequence(argv[1], argv[2])

    print "Len Percentage Hits  Misses"
    for key, value in hpl.iteritems():
        print "{:<3} {:>10.2f} {:>4} {:>4}".format(key, float(value[1]*100)/float(sum(value)), value[0], value[1])

    if len(nohits) > 0:
        print "Statistics for fragments with no hits\n-----------------------------"
        seqlengths = [len(seq) for _, seq in nohits]
        print "Min length: {}".format(min(seqlengths))
        print "Max length: {}".format(max(seqlengths))
        print "Avg length: {}".format(float(sum(seqlengths)) / len(nohits))
    
