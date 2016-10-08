"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Wanding Zhou, Tenghui Chen, Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
import sys, re, argparse
from .mutation import parser_add_mutation, parse_tok_mutation_str, list_parse_mutation
from .parser import parser_add_annotation
from .annodb import AnnoDB
from .transcripts import *
from .err import *
from .utils import *
from .config import read_config
from .snv import __core_annotate_codon_snv
from .record import Query, QueryREG

outformat="{altid}\t{chrm}\t{codon1}\t{codon2}\t{tptstr}"

def _main_core_(args, q, db):

    k2transcripts = {}
    if isinstance(q, QueryREG):
        q.pos = q.beg
        q.ref = q.refseq
        q.alt = ''

    for t1, c1 in __core_annotate_codon_snv(args, q, db):
        # search any of the 3 positions
        for cind in range(3):
            gpos = c1.locs[cind]
            for t2 in db.get_transcripts(t1.chrm, gpos):
                c2, p = t2.gpos2codon(gpos)
                if t1 == t2: continue
                if p.tpos != 0: continue
                # if c2.region != 'coding': continue
                if c1.index == c2.index: continue
                if len(c2.seq) != 3: continue # often due to last incomplete codon
                if q.ref and q.ref != codon2aa(c2.seq): continue
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

    for k, tpairs in k2transcripts.items():
        altid, chrm, c1, c2 = k
        if q.op: s = q.op+'\t'
        else: s = ''
        s += outformat.format(altid=altid, tptstr=','.join(tpairs), chrm=chrm,
                              codon1='-'.join(map(str,c1)), codon2='-'.join(map(str,c2)))
        print(s)

def main_list(args, db): #name2gene, thash):

    if not args.noheader:
        print('origin_id\talt_id\tchrm\tcodon1\tcodon2\ttranscripts_choice')
    for q, line in list_parse_mutation(args, 'p'):

        genefound = False
        for q.gene in db.get_gene(q.tok):
            genefound = True
            try:
                _main_core_(args, q, db)
            except UnImplementedError as e:
                wrap_exception(e, q.op, args)
            except SequenceRetrievalError as e:
                wrap_exception(e, q.op, args)
        if not genefound:
            err_warn('gene %s is not recognized.' % q.tok)


def main_one(args, db): #name2gene, thash):

    if not args.noheader:
        print('origin_id\talt_id\tchrm\tcodon1\tcodon2\ttranscripts_choice')
    q = parse_tok_mutation_str(args.i, 'p')
    q.op = args.i
    genefound = False
    for q.gene in db.get_gene(q.tok):
        genefound = True
        q.op = args.i
        _main_core_(args, q, db)
        
    if not genefound:
        err_die('gene %s is not recognized.' % q.tok)
        return

def main(args):

    config = read_config()
    db = AnnoDB(args, config)
    # name2gene, thash = parse_annotation(args)

    if args.l:
        main_list(args, db) #name2gene, thash)
    if args.i:
        main_one(args, db) #name2gene, thash)

def add_parser_codonsearch(subparsers, config):

    parser = subparsers.add_parser('codonsearch', help="search equivalent codon representations")
    parser_add_mutation(parser)
    parser_add_annotation(parser)
    parser.set_defaults(func=main)
