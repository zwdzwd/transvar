"""
search alternative codonpositions due to different transcript usage
"""
import sys, re, argparse
from mutation import parser_add_mutation, parse_tok_mutation_str, list_parse_mutation
from transcripts import *
from utils import *
from revanno_snv import __core_annotate_codon_snv
from anno import pos2codon
from record import Query

outformat="{altid}\t{chrm}\t{codon1}\t{codon2}\t{tptstr}"

def _main_core_(args, q, thash):

    k2transcripts = {}
    if q.is_codon:
        for t1, c1 in __core_annotate_codon_snv(args, q):
            for cind in xrange(3):
                for t2, c2 in pos2codon(thash, t1.chrm, c1.locs[cind]):
                    if t1 == t2: continue
                    if c2.region != 'coding': continue
                    if c1.index == c2.index: continue
                    if q.ref and q.ref != standard_codon_table[c2.seq]: continue
                    altid = t1.gene.name+'.p.'
                    if q.ref: altid += q.ref
                    altid += str(c2.index)
                    k = (altid, c1.chrm, tuple(c1.locs), tuple(c2.locs))
                    tpair = '%s[%s]/%s[%s]' % (t1.name, t1.source, t2.name, t2.source)
                    if k in k2transcripts:
                        if tpair not in k2transcripts[k]:
                            k2transcripts[k].append(tpair)
                    else:
                        k2transcripts[k] = [tpair]

    for k, tpairs in k2transcripts.iteritems():
        altid, chrm, c1, c2 = k
        if q.op: s = q.op+'\t'
        else: s = ''
        s += outformat.format(altid=altid, tptstr=','.join(tpairs), chrm=chrm,
                              codon1='-'.join(map(str,c1)), codon2='-'.join(map(str,c2)))
        print s

def main(args):

    name2gene, thash = parse_annotation(args)

    if args.l:
        main_list(args, name2gene, thash)
    if args.i:
        main_one(args, name2gene, thash)

def main_list(args, name2gene, thash):

    for q, line in list_parse_mutation(args):

        if q.tok not in name2gene:
            sys.stderr.write("Gene %s is not recognized.\n" % q.tok)
            continue
        q.gene = name2gene[q.tok]
        try:
            _main_core_(args, q, thash)
        except UnImplementedError as e:
            err_print(line)
            raise e

def main_one(args, name2gene, thash):

    q = parse_tok_mutation_str(args.i)
    q.op = args.i
    if q.tok not in name2gene:
        sys.stderr.write("Gene %s not recognized.\n" % q.tok)
        return
        return

    q.gene = name2gene[q.tok]
    q.op = args.i

    _main_core_(args, q, thash)

def add_parser_codonsearch(subparsers, d):

    parser = subparsers.add_parser('codonsearch', help=__doc__)
    parser_add_mutation(parser)
    parser_add_annotation(parser, d)
    parser.set_defaults(func=main)
