"""
test whether two codon annotations may be the same nucleotide position.
"""
import sys, argparse, re
from utils import *

def main(args):

    codon1, codon2 = args.codons
    name2gene, thash = parse_annotation(args)

    gene_name, codon_pos = codon1.split('.p')
    print gene_name, codon_pos
    codon_pos = int(codon_pos)
    gene = name2gene[gene_name]
    locs1 = set()
    for i, t in enumerate(gene.tpts):
        codon = t.cpos2codon(codon_pos)
        if codon.index > 0:
            print "transcript [%s] %d \tcodon: %s" % (t.name, i, "-".join(map(str, codon.locs)))
            locs1.add(codon.locs)

    gene_name, codon_pos = codon2.split('.p')
    print gene_name, codon_pos
    codon_pos = int(codon_pos)
    gene = name2gene[gene_name]
    locs2 = set()
    for i, t in enumerate(gene.tpts):
        codon = t.cpos2codon(codon_pos)
        if codon.index > 0:
            print "transcript [%s] %d \tcodon: %s" % (t.name, i, ",".join(map(str, codon.locs)))
            locs2.add(codon.locs)

    if locs1 & locs2:
        print "Genomic location might be the same."
    else:
        print "Genomic location shouldn't be the same."


def add_parser_codoneq(subparsers):

    parser = subparsers.add_parser("codoneq", help=__doc__)
    parser.add_argument("-c", nargs=2, help="two codons to test. Format: [gene_name].p[codon_location], e.g., MET.p1010", dest='codons')
    parser_add_annotation(parser)
    parser.set_defaults(func=main)
