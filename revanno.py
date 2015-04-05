""" 
annotate a codon position or an amino acid change
"""
import sys, argparse, re
from transcripts import *
from utils import *
from record import *
from config import read_config
from mutation import parser_add_mutation, parse_tok_mutation_str, list_parse_mutation
from revanno_snv import _core_annotate_nuc_snv, _core_annotate_codon_snv
from revanno_del import _core_annotate_nuc_del, _core_annotate_codon_del
# from revanno_ins import _core_annotate_nuc_ins, _core_annotate_codon_ins
from insertion import annotate_insertion_protein, annotate_insertion_cdna, annotate_duplication_cdna
from revanno_mnv import _core_annotate_nuc_mnv, _core_annotate_codon_mnv
from revanno_fs import _core_annotate_codon_fs
# from revanno_dup import _core_annotate_nuc_dup
from revanno_reg import _core_annotate_nuc_reg, _core_annotate_codon_reg

def _core_annotate_codon(args, q, db):

    if args.longest: tpts = [q.gene.longest_tpt()]
    else: tpts = q.gene.tpts
    if isinstance(q, QuerySNV):
        return _core_annotate_codon_snv(args, q, tpts, db)
    elif isinstance(q, QueryDEL):
        return _core_annotate_codon_del(args, q, tpts, db)
    elif isinstance(q, QueryINS):
        return annotate_insertion_protein(args, q, tpts, db)
    elif isinstance(q, QueryMNV):
        return _core_annotate_codon_mnv(args, q, tpts, db)
    elif isinstance(q, QueryFrameShift):
        return _core_annotate_codon_fs(args, q, tpts, db)
    else:
        return _core_annotate_codon_reg(args, q, tpts, db)

def _core_annotate_nuc(args, q, db):

    if args.longest: tpts = [q.gene.longest_tpt()]
    else: tpts = q.gene.tpts

    if isinstance(q, QuerySNV):
        return _core_annotate_nuc_snv(args, q, tpts, db)
    elif isinstance(q, QueryDEL):
        return _core_annotate_nuc_del(args, q, tpts, db)
    elif isinstance(q, QueryINS):
        return annotate_insertion_cdna(args, q, tpts, db)
    elif isinstance(q, QueryMNV):
        return _core_annotate_nuc_mnv(args, q, tpts, db)
    elif isinstance(q, QueryDUP):
        return annotate_duplication_cdna(args, q, tpts, db)
    else:
        return _core_annotate_nuc_reg(args, q, tpts, db)

def _main_core_(args, q, db):

    if q.is_codon:                # in codon coordinate
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
    replace_defaults(args, config)
    db = AnnoDB(args)

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


