#!/usr/bin/env python
# Fredrik Boulund 2015
# Sample sequences from a FASTA file 

from read_fasta import read_fasta
from sys import argv, exit, maxint
import argparse
from random import sample, choice as pychoice
from numpy.random import binomial, choice


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
    parser.add_argument("-p", metavar="P", type=float,
            default=0.0,
            help="Probability of mutation (per amino acid) [%(default)s].")
    parser.add_argument("--matrix", metavar="M", 
            default="/local/blast-2.2.26/data/PAM30",
            help="Location of BLAST reference PAM or BLOSUM matrix to use for point mutations [%(default)s].")
    parser.add_argument("-o", "--outfile", metavar="FILE", dest="outfile",
            default="",
            help="Write output to FILE instead of STDOUT.")

    if len(argv)<2:
        parser.print_help()
        exit()
    
    options = parser.parse_args()
    return options


def sample_fasta(fastafile, outfile, matrix, options):
    """Sample sequences from FASTA.
    """

    seqs = []
    for header, seq in read_fasta(fastafile, keep_formatting=False):
        seqlen = len(seq)
        if not options.maxlength:
            options.maxlength = maxint
        if seqlen >= options.minlength and seqlen <= options.maxlength:
            seqs.append((header,seq))

    if options.outfile:
        with open(outfile, 'w') as f:
            for n in xrange(0,options.n):
                header, seq = pychoice(seqs)
                if options.p > 0.0:
                    seq, mutations = mutate_seq(seq, options.p, matrix)
                    f.write(">"+header+" mutations=", str(mutations), "\n")
                else:
                    f.write(">"+header+"\n")
                f.write(seq+"\n")
    else:
        for n in xrange(0,options.n):
            header, seq = pychoice(seqs)
            if options.p > 0.0:
                seq, mutations = mutate_seq(seq, options.p, matrix)
                print ">"+header+" mutations="+str(mutations)
            else:
                print ">"+header
            print seq


def mutate(aa, m):
    """Mutate a single amino acid.
    """
    likelihoods = matrix[aa].items()
    minlikelihood = min(likelihoods, key=lambda v: v[1])[1]
    adjusted_likelihoods = [l+abs(minlikelihood) for a, l in likelihoods]
    normalizer = sum(adjusted_likelihoods)
    probs = [float(l)/normalizer for l in adjusted_likelihoods]
    aas = [amino_acid for amino_acid, likelihood in likelihoods]
    new_aa = choice(aas, p=probs)
    if new_aa == aa:
        return mutate(aa, m)
    else:
        return new_aa

    
def mutate_seq(seq, p, m):
    """Mutate sequence as positions chosen by sampling binomial with probability p,
    using the substitution matrix m.
    """
    mutations = binomial(len(seq), p)
    positions = sample(xrange(len(seq)), mutations)
    seq = list(seq.upper())
    for pos in positions:
        seq[pos] = mutate(seq[pos], m)
    seq = ''.join(seq)
    return seq, mutations


def read_substitution_matrix(filename, normalize=False, remove=[]):
    """Read substitution matrix from filename into a nested dictionary.
    The likelihoods can be normalized to "probabilities".
    """
    with open(filename) as f:
        line = f.readline()
        if line.startswith("# Entries"):
            matrix = {}
            aas = f.readline().split()
            for aa1 in aas:
                for aa2 in aas:
                    matrix[aa1] = {aa2: 0}
            for line in f:
                likelihoods = line.split()
                cur_aa = likelihoods[0]
                likelihoods = likelihoods[1:]
                for aa in aas:
                    matrix[cur_aa][aa] = int(likelihoods.pop(0))
        else:
            raise IOError("{} doesn't appear to be a BLAST substitution matrix.".format(filename))
    for code in remove:
        matrix.pop(code, None)
        for subdict in matrix.itervalues():
            subdict.pop(code, None)
    return matrix


if __name__ == "__main__":
    options = parse_args(argv)
    matrix = read_substitution_matrix(options.matrix, remove=["X", "*"])
    sample_fasta(options.FASTA, options.outfile, matrix, options)
