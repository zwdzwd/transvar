"""
annotate nucleotide position(s) or mutations
"""
import sys, argparse, re
from transcripts import *
from utils import *
from record import *
from config import read_config
from mutation import parse_tok_mutation_str, list_parse_mutation, parser_add_mutation

from mnv import annotate_mnv_gdna
from snv import annotate_snv_gdna
from insertion import annotate_insertion_gdna
from deletion import annotate_deletion_gdna
from region import annotate_region_gdna

def _main_core_(args, db, q):

    if isinstance(q, QuerySNV):
        return annotate_snv_gdna(args, q, db)
    elif isinstance(q, QueryDEL):
        return annotate_deletion_gdna(args, q, db)
    elif isinstance(q, QueryINS):
        return annotate_insertion_gdna(args, q, db)
    elif isinstance(q, QueryMNV):
        return annotate_mnv_gdna(args, q, db)
    else:
        return annotate_region_gdna(args, q, db)

def main_list(args, db):

    for q, line in list_parse_mutation(args, muttype='g'):
        q.tok = normalize_chrm(q.tok)
        try:
            _main_core_(args, db, q)
        except IncompatibleTranscriptError:
            err_print(line)
            raise Exception()

def main_one(args, db):

    q = parse_tok_mutation_str(args.i, muttype='g')
    q.op = args.i
    q.tok = normalize_chrm(q.tok)
    _main_core_(args, db, q)

def main(args):

    config = read_config()
    db = AnnoDB(args, config)

    if not args.noheader:
        print_header()
    
    if args.l:
        main_list(args, db)

    if args.i:
        main_one(args, db)

def add_parser_anno(subparsers, config):

    parser = subparsers.add_parser('anno', help=__doc__)
    parser_add_annotation(parser)
    parser_add_mutation(parser)
    parser.add_argument('--longest', action="store_true",
                        help='use longest transcript instead of reporting all transcripts')

    parser.set_defaults(func=main)
