"""
annotate nucleotide position(s) or mutations
"""
import sys, argparse, re
from transcripts import *
from utils import *
from config import read_config
from mutation import parse_tok_mutation_str, list_parse_mutation, parser_add_mutation
from record import *
from anno_reg import _annotate_reg
from anno_snv import _annotate_snv
from anno_del import _annotate_del
from anno_ins import _annotate_ins

def _main_core_(args, db, q):

    if isinstance(q, QuerySNV):
        return _annotate_snv(args, q, db)
    elif isinstance(q, QueryDEL):
        return _annotate_del(args, q, db)
    elif isinstance(q, QueryINS):
        return _annotate_ins(args, q, db)
    elif isinstance(q, QueryMNV):
        return _annotate_mnv(args, q, db)
    else:
        return _annotate_reg(args, q, db)

def main_list(args, db):

    for q, line in list_parse_mutation(args, muttype='g'):
        q.tok = normalize_chrm(q.tok)
        _main_core_(args, db, q)
        # _main_core_(args, thash, q)

def main_one(args, db):

    q = parse_tok_mutation_str(args.i, muttype='g')
    q.op = args.i
    q.tok = normalize_chrm(q.tok)
    # _main_core_(args, thash, q)
    _main_core_(args, db, q)

def main(args):

    config = read_config()
    replace_defaults(args, config)
    db = AnnoDB(args)
    # name2gene, thash = parse_annotation(args)

    if args.l:
        main_list(args, db) #thash)

    if args.i:
        main_one(args, db) #thash)

# def main(args):

#     name2gene, thash = parse_annotation(args)
#     for line in args.l:
#         fields = line.strip().split('\t')
#         name = fields[int(args.i)-1]
#         tn = fields[int(args.t)-1]
#         if name in name2gene:
#             gene = name2gene[name]
#             ts = [t for t in gene.tpts if tn == t.name]
#             if ts:
#                 o = len([t for t in gene.tpts if len(t) > len(ts[0])])
#                 print o
    
def add_parser_anno(subparsers, config):

    parser = subparsers.add_parser('anno', help=__doc__)
    parser_add_annotation(parser)
    parser_add_mutation(parser)
    parser.add_argument('--longest', action="store_true",
                        help='use longest transcript instead of reporting all transcripts')

    parser.set_defaults(func=main)
