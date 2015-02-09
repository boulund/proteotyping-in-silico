#!/usr/bin/env python
# Fredrik Boulund 2015
# Sample sequences from a FASTA file 

from read_fasta import read_fasta
from sys import argv, exit, maxint
import argparse
from random import choice


def parse_args(argv):
    """Parse commandline arguments.
    """

    desc = """Sample sequences from FASTA files with replacement. Fredrik Boulund 2015"""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("FASTA", 
            help="FASTA file to sample from.")
    parser.add_argument("-n", metavar="N", required=True, type=int,
            help="Number of sequences to sample from FASTA file [%(default)s].")
    parser.add_argument("--maxlength", metavar="M", type=int,
            default=0,
            help="Maximum length of sequences to sample from, 0 means no limit [%(default)s], cant be bigger than {}.".format(maxint))
    parser.add_argument("--minlength", metavar="m", type=int,
            default=0,
            help="Minimum length of sequences to sample from, 0 means no limit [%(default)s].")
    parser.add_argument("-o", "--outfile", metavar="FILE",
            default="",
            help="Write output to FILE instead of STDOUT.")

    if len(argv)<2:
        parser.print_help()
        exit()
    
    options = parser.parse_args()
    return options


def sample_fasta(fastafile, outfile, options):
    """Sample sequences from FASTA.
    """

    seqs = []
    for header, seq in read_fasta(fastafile):
        seqlen = len(seq)
        if not options.maxlength:
            options.maxlength = maxint
        if seqlen >= options.minlength and seqlen <= options.maxlength:
            seqs.append((header,seq))

    for n in xrange(0,options.n):
        header, seq = choice(seqs)
        print ">"+header
        print seq


if __name__ == "__main__":
    options = parse_args(argv)
    sample_fasta(options.FASTA, options.outfile, options)
