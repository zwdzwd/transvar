""" 
annotate a codon position or an amino acid change
"""
import sys, argparse, re
from transcripts import *
from utils import *
from record import *
from mutation import parser_add_mutation, parse_tok_mutation_str, list_parse_mutation
from revanno_snv import _core_annotate_nuc_snv, _core_annotate_codon_snv
from revanno_del import _core_annotate_nuc_del, _core_annotate_codon_del
from revanno_ins import _core_annotate_nuc_ins, _core_annotate_codon_ins
from revanno_mnv import _core_annotate_nuc_mnv, _core_annotate_codon_mnv
from revanno_fs import _core_annotate_codon_fs
from revanno_dup import _core_annotate_nuc_dup
from revanno_reg import _core_annotate_nuc_reg

def _core_annotate_codon(args, q):

    if args.longest: tpts = [q.gene.longest_tpt()]
    else: tpts = q.gene.tpts

    if isinstance(q, QuerySNV):
        return _core_annotate_codon_snv(args, q, tpts)
    elif isinstance(q, QueryDEL):
        return _core_annotate_codon_del(args, q, tpts)
    elif isinstance(q, QueryINS):
        return _core_annotate_codon_ins(args, q, tpts)
    elif isinstance(q, QueryMNV):
        return _core_annotate_codon_mnv(args, q, tpts)
    elif isinstance(q, QueryFrameShift):
        return _core_annotate_codon_fs(args, q, tpts)

def _core_annotate_nuc(args, q):

    if args.longest: tpts = [q.gene.longest_tpt()]
    else: tpts = q.gene.tpts

    if isinstance(q, QuerySNV):
        return _core_annotate_nuc_snv(args, q, tpts)
    elif isinstance(q, QueryDEL):
        return _core_annotate_nuc_del(args, q, tpts)
    elif isinstance(q, QueryINS):
        return _core_annotate_nuc_ins(args, q, tpts)
    elif isinstance(q, QueryMNV):
        return _core_annotate_nuc_mnv(args, q, tpts)
    elif isinstance(q, QueryDUP):
        return _core_annotate_nuc_dup(args, q, tpts)
    else:
        return _core_annotate_nuc_reg(args, q, tpts)

def _main_core_(args, q):

    if q.is_codon:                # in codon coordinate
        _core_annotate_codon(args, q)
    else:                       # in nucleotide coordinate
        _core_annotate_nuc(args, q)

def main_list(args, name2gene):

    for q, line in list_parse_mutation(args):
        q.tok = q.tok.upper()
        if q.tok not in name2gene:
            r = Record()
            r.info = "GeneNotRecognized"
            r.format(q.op)
            continue
        q.gene = name2gene[q.tok]
        try:
            _main_core_(args, q)
        except UnImplementedError as e:
            err_print(line)
            raise e
        except Exception as e:
            err_print('line:'+line)
            raise e

def main_one(args, name2gene):

    q = parse_tok_mutation_str(args.i)
    q.tok = q.tok.upper()
    if not q: return

    if q.tok not in name2gene:
        r = Record()
        r.info = 'GeneNotRecognized'
        r.format(q.op)
        return

    q.gene = name2gene[q.tok]
    q.op = args.i

    _main_core_(args, q)

def main(args):

    name2gene, thash = parse_annotation(args)

    if args.l:
        main_list(args, name2gene)

    if args.i:
        main_one(args, name2gene)

def add_parser_revanno(subparsers, d):

    parser = subparsers.add_parser("revanno", help=__doc__)
    parser_add_annotation(parser, d)
    parser_add_mutation(parser)
    parser.add_argument("--longest", action="store_true",
                        help="consider only longest transcript")
    parser.set_defaults(func=main)


