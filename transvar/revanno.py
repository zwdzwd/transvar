""" 
annotate a codon position or an amino acid change
"""
import sys, argparse, re
from transcripts import *
from utils import *
from record import *
from config import read_config
from mutation import parser_add_mutation, parse_tok_mutation_str, list_parse_mutation

from snv import 
from deletion import 
from insertion import
from mnv import
from region import 
from frameshift import annotate_frameshift

def _core_annotate_codon(args, q, db):

    if args.longest: tpts = [q.gene.longest_tpt()]
    else: tpts = q.gene.tpts
    

def _core_annotate_nuc(args, q, db):

    

def _main_core_(args, q, db):


        
    elif q.is_codon:              # in codon coordinate
        _core_annotate_codon(args, q, db)
    else:                       # in nucleotide coordinate
        _core_annotate_nuc(args, q, db)

def main_list(args, db):

    for q, line in list_parse_mutation(args):
        q.tok = q.tok.upper()
        q.gene = db.get_gene(q.tok)
        if not q.gene:
            r = Record()
            r.info = "GeneNotRecognized"
            r.format(q.op)
            continue
        try:
            _main_core_(args, q, db)
        except UnImplementedError as e:
            err_print(line)
            raise e
        except Exception as e:
            err_print('line:'+line)
            raise e

def main_one(args, db):

    try:
        q = parse_tok_mutation_str(args.i)
    except InvalidInputError:
        r = Record()
        r.append_info('invalid_mutation_string_%s' % args.i)
        r.format(args.i)
        return

    if not q:
        return
    
    q.tok = q.tok.upper()
    q.gene = db.get_gene(q.tok)
    if not q.gene:
        r = Record()
        r.append_info('gene_not_recognized_(%s)' % q.tok)
        err_warn('gene %s not recognized. make sure the right (if any) transcript database is used.' % q.tok)
        r.format(q.op)
        return

    q.op = args.i
    _main_core_(args, q, db)

def main(args):

    config = read_config()
    db = AnnoDB(args, config)

    if not args.noheader:
        print_header()

    if args.l:
        main_list(args, db)

    if args.i:
        main_one(args, db)

def add_parser_revanno(subparsers, config):

    parser = subparsers.add_parser("revanno", help=__doc__)
    parser_add_annotation(parser)
    parser_add_mutation(parser)
    parser.add_argument("--longest", action="store_true",
                        help="consider only longest transcript")
    parser.set_defaults(func=main)


