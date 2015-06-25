"""
The MIT License

Copyright (c) 2015 by The University of Texas MD Anderson Cancer Center (kchen3@mdanderson.org)

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

import sys, argparse, re
from annodb import AnnoDB
# from transcripts import *
from parser import parser_add_annotation
# import parser
from record import *
from err import *
from config import read_config
from mutation import parse_tok_mutation_str, list_parse_mutation, vcf_parse_mutation, parser_add_mutation

from mnv import annotate_mnv_gdna, annotate_mnv_protein, annotate_mnv_cdna
from snv import annotate_snv_gdna, annotate_snv_protein, annotate_snv_cdna
from insertion import annotate_insertion_gdna, annotate_insertion_protein, annotate_insertion_cdna, annotate_duplication_cdna
from deletion import annotate_deletion_gdna, annotate_deletion_protein, annotate_deletion_cdna
from region import annotate_region_gdna, annotate_region_protein, annotate_region_cdna, annotate_gene
from frameshift import annotate_frameshift
from functools import partial

def _main_core_(args, q, db, at):

    if at == 'c':

        if args.longest:
            tpts = [q.gene.longest_tpt()]
        elif args.longestcoding:
            tpts = [q.gene.longest_coding_tpt()]
        else: tpts = q.gene.tpts
        
        if isinstance(q, QueryGENE):
            return annotate_gene(args, q, tpts, db)
        elif isinstance(q, QuerySNV):
            return annotate_snv_cdna(args, q, tpts, db)
        elif isinstance(q, QueryDEL):
            return annotate_deletion_cdna(args, q, tpts, db)
        elif isinstance(q, QueryINS):
            return annotate_insertion_cdna(args, q, tpts, db)
        elif isinstance(q, QueryMNV):
            return annotate_mnv_cdna(args, q, tpts, db)
        elif isinstance(q, QueryDUP):
            return annotate_duplication_cdna(args, q, tpts, db)
        elif isinstance(q, QueryREG):
            return annotate_region_cdna(args, q, tpts, db)
        else:
            raise Exception('mutation type inference error for %s' % q.op)
        
    elif at == 'p':

        if args.longest or args.longestcoding:
            tpts = [q.gene.longest_coding_tpt()]
        else:
            tpts = q.gene.coding_tpts()
        
        if isinstance(q, QueryGENE):
            return annotate_gene(args, q, tpts, db)
        elif isinstance(q, QuerySNV):
            return annotate_snv_protein(args, q, tpts, db)
        elif isinstance(q, QueryDEL):
            return annotate_deletion_protein(args, q, tpts, db)
        elif isinstance(q, QueryINS):
            return annotate_insertion_protein(args, q, tpts, db)
        elif isinstance(q, QueryMNV):
            return annotate_mnv_protein(args, q, tpts, db)
        elif isinstance(q, QueryFrameShift):
            return annotate_frameshift(args, q, tpts, db)
        elif isinstance(q, QueryREG):
            return annotate_region_protein(args, q, tpts, db)
        else:
            raise Exception('mutation type inference error for %s' % q.op)

    elif at == 'g':

        if isinstance(q, QuerySNV):
            return annotate_snv_gdna(args, q, db)
        elif isinstance(q, QueryDEL):
            return annotate_deletion_gdna(args, q, db)
        elif isinstance(q, QueryINS):
            return annotate_insertion_gdna(args, q, db)
        elif isinstance(q, QueryMNV):
            return annotate_mnv_gdna(args, q, db)
        elif isinstance(q, QueryREG):
            return annotate_region_gdna(args, q, db)
        else:
            raise Exception('mutation type inference error for %s' % q.op)

    
def main_list(args, db, at, mutation_parser):

    for q, line in mutation_parser:

        if at == 'g':
            q.tok = normalize_chrm(q.tok)
        else:
            q.tok = q.tok.upper()
            q.gene = db.get_gene(q.tok)
            if not q.gene:
                r = Record()
                r.append_info('gene_not_recognized_(%s)' % q.tok)
                err_warn('gene %s not recognized. make sure the right (if any) transcript database is used.' % q.tok)
                r.format(q.op)
                continue
            
        # try:
        _main_core_(args, q, db, at)
        # except:
        # err_print('exception %s' % line)
        # raise Exception()

def main_one(args, db, at):

    try:
        q = parse_tok_mutation_str(args.i, at)
    except InvalidInputError:
        err_die('invalid mutation string %s (type:%s)' % (args.i, at))

    if not q:
        return

    if at == 'g':
        q.tok = normalize_chrm(q.tok)
    else:
        q.tok = q.tok.upper()
        q.gene = db.get_gene(q.tok)
        if not q.gene:
            err_die('gene %s not recognized. make sure the right (if any) transcript database is used.' % q.tok)

    q.op = args.i
    _main_core_(args, q, db, at)

def main(args, at):

    config = read_config()
    db = AnnoDB(args, config)

    if (not args.vcf) and (not args.noheader):
        print print_header()

    if args.l:
        main_list(args, db, at, list_parse_mutation(args, at))

    if args.vcf:
        if at != 'g':
            err_raise("can apply on ganno to VCF input")
        main_list(args, db, at, vcf_parse_mutation(args, 'g'))

    if args.i:
        main_one(args, db, at)

def add_parser_anno(subparsers, config):

    parser = subparsers.add_parser('ganno', help='annotate gDNA element')
    parser_add_annotation(parser)
    parser_add_mutation(parser)
    parser.set_defaults(func=partial(main, at='g'))

    parser = subparsers.add_parser("canno", help='annotate cDNA elements')
    parser_add_annotation(parser)
    parser_add_mutation(parser)
    parser.set_defaults(func=partial(main, at='c'))

    parser = subparsers.add_parser("panno", help='annotate protein element')
    parser_add_annotation(parser)
    parser_add_mutation(parser)
    parser.set_defaults(func=partial(main, at='p'))

