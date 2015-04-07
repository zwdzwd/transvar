""" 
annotate a codon position or an amino acid change
"""
import sys, argparse, re
from transcripts import *
from utils import *
from record import *
from config import read_config
from mutation import parser_add_mutation, parse_tok_mutation_str, list_parse_mutation

from snv import annotate_snv_protein, annotate_snv_cdna
from deletion import annotate_deletion_protein, annotate_deletion_cdna
from insertion import annotate_insertion_protein, annotate_insertion_cdna, annotate_duplication_cdna
from mnv import annotate_mnv_protein, annotate_mnv_cdna
from region import annotate_region_protein, annotate_region_cdna
from frameshift import annotate_frameshift

def _core_annotate_codon(args, q, db):

    if args.longest: tpts = [q.gene.longest_tpt()]
    else: tpts = q.gene.tpts
    if isinstance(q, QuerySNV):
        return annotate_snv_protein(args, q, tpts, db)
    elif isinstance(q, QueryDEL):
        return annotate_deletion_protein(args, q, tpts, db)
    elif isinstance(q, QueryINS):
        return annotate_insertion_protein(args, q, tpts, db)
    elif isinstance(q, QueryMNV):
        return annotate_mnv_protein(args, q, tpts, db)
    elif isinstance(q, QueryFrameShift):
        return annotate_frameshift(args, q, tpts, db)
    else:
        return annotate_region_protein(args, q, tpts, db)

def _core_annotate_nuc(args, q, db):

    if args.longest: tpts = [q.gene.longest_tpt()]
    else: tpts = q.gene.tpts

    if isinstance(q, QuerySNV):
        return annotate_snv_cdna(args, q, tpts, db)
    elif isinstance(q, QueryDEL):
        return annotate_deletion_cdna(args, q, tpts, db)
    elif isinstance(q, QueryINS):
        return annotate_insertion_cdna(args, q, tpts, db)
    elif isinstance(q, QueryMNV):
        return annotate_mnv_cdna(args, q, tpts, db)
    elif isinstance(q, QueryDUP):
        return annotate_duplication_cdna(args, q, tpts, db)
    else:
        return annotate_region_cdna(args, q, tpts, db)

def _main_core_(args, q, db):

    if q.is_codon:              # in codon coordinate
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

    q = parse_tok_mutation_str(args.i)
    q.tok = q.tok.upper()
    if not q: return
    q.gene = db.get_gene(q.tok)
    if not q.gene:
        r = Record()
        r.info = 'GeneNotRecognized'
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


