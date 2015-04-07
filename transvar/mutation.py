import re, sys, argparse
from utils import *
from record import *
from err import *

def parse_mutation_str(mut_str, muttype=None):
    mut_str = mut_str.strip()
    mp = re.match(r'(p)?(\.)?([A-Z*?]*)(\d+)(_([A-Z*?]*)(\d+))?(del([A-Z*?\d]*))?(ins([A-Z*?]+))?>?([A-Z*?]+)?(fs\*(\d+))?(ref([A-Zx*]*))?$', mut_str)
    mn = re.match(r'(c)?(\.)?([\d+-]+)(_([\d+-]+))?(\.)?(del([atgcATGC\d]*))?(ins([atgcATGC]*))?(([atgcATGC?]*)>([atgcATGC?]*))?(dup([atgcATGC\d]*))?$', mut_str)
    mg = re.match(r'(g)?(\.)?(\d+)(_(\d+))?(\.)?(del([atgcATGC\d]*))?(ins([atgcATGC]*))?(([atgcATGC?]*)>([atgcATGC?]*))?(dup([atgcATGC\d]*))?$', mut_str)
    
    if (mp and not muttype) or muttype=='p':
        #print mp.groups()
        if not mp:
            raise InvalidInputError("Invalid protein level identifier: %s" % mut_str)
        (_, _, _beg_aa, _beg_i, _end_s, _end_aa, _end_i, 
         _is_del, _d, _is_ins, _i, _alt, _is_fs, _stop_i, _has_ref, _ref) = mp.groups()
        if _is_fs:
            #print 'fs'
            q = QueryFrameShift()
            q.pos = int(_beg_i)
            q.ref = _beg_aa
            q.alt = _alt if _alt else ''
            q.stop_index = int(_stop_i)
        elif _is_del and not _is_ins:
            #print 'del'
            q = QueryDEL()
            q.beg = int(_beg_i)
            q.end = int(_end_i) if _end_i else q.beg
            q.beg_aa = _beg_aa
            q.end_aa = _end_aa
            if _d: q.delseq = '' if _d.isdigit() else _d.upper()
        elif (_is_del and _is_ins and 
              (_d == '1' or (_d and len(_d) == 1) or not _end_s)
              and len(_i) == 1):
            #print 'snv'
            q = QuerySNV()
            q.pos = int(_beg_i)
            q.ref = _beg_aa
            q.alt = _alt
        elif _is_del and _is_ins:
            #print 'mnv'
            q = QueryMNV()
            q.beg = int(_beg_i)
            q.end = int(_end_i) if _end_i else q.beg
            q.beg_aa = _beg_aa.upper() if _beg_aa else ''
            q.end_aa = _end_aa.upper() if _end_aa else ''
            q.refseq = _d.upper() if _d else ''
            q.altseq = _i.upper() if _i else ''
        elif _is_ins and not _is_del:
            # print 'ins'
            q = QueryINS()
            q.beg = int(_beg_i)
            q.beg_aa = _beg_aa
            q.end = int(_end_i)
            q.end_aa = _end_aa
            q.insseq = _i.upper() if _i else ''
        elif not _end_s and (not _beg_aa or len(_beg_aa) == 1) and _alt:
            # print 'snv'
            q = QuerySNV()
            q.pos = int(_beg_i)
            q.ref = _beg_aa
            q.alt = _alt.upper()
        elif _i:
            # print 'mnv'
            q = QueryMNV()
            q.beg = int(_beg_i)
            q.end = int(_end_i) if _end_i else q.beg
            q.beg_aa = _beg_aa.upper() if _beg_aa else ''
            q.end_aa = _end_aa.upper() if _end_aa else ''
            q.refseq = _d.upper() if _d else ''
            q.altseq = _i.upper() if _i else ''
        else:
            # print 'reg'
            q = QueryREG()
            q.beg = int(_beg_i)
            q.end = int(_end_i) if _end_i else q.beg
            q.beg_aa = _beg_aa.upper() if _beg_aa else ''
            q.end_aa = _end_aa.upper() if _end_aa else ''
            q.refseq = _ref if _ref else ''
            if (not q.refseq) and q.beg == q.end:
                q.refseq = q.beg_aa

        q.is_codon = True
    elif (mn and not muttype) or muttype=='n':
        if not mn:
            raise InvalidInputError("Invalid cDNA level input: %s" % mut_str)
        (_, _, _beg, _end_s, _end, _, _is_del, _d,
         _is_ins, _i, _is_sub, _ref, _alt, _is_dup, _dupseq) = mn.groups()
        if _is_sub and len(_ref) <= 1 and len(_alt) <= 1:
            q = QuerySNV()
            q.pos = parse_pos(_beg)
            q.ref = _ref if _ref and _ref != '?' else ''
            q.alt = _alt if _alt and _alt != '?' else ''
        elif (_is_del and _is_ins and
              (_d == '1' or len(_d) == 1) and len(_i) == 1):
            q = QuerySNV()
            q.pos = parse_pos(_beg)
            q.ref = '' if _d.isdigit() else _d.upper()
            q.alt = _i.upper() if _i else ''
        elif _is_del and not _is_ins:
            q = QueryDEL()
            q.beg = parse_pos(_beg)
            q.end = parse_pos(_end) if _end else q.beg
            q.delseq = '' if _d.isdigit() else _d.upper()
        elif _is_ins and not _is_del:
            q = QueryINS()
            q.pos = parse_pos(_beg)
            if _i: q.insseq = _i.upper()
            else: err_die('insertion without inserted sequence: %s.' % mut_str)
        elif _is_ins and _is_del:
            q = QueryMNV()
            q.beg = parse_pos(_beg)
            q.end = parse_pos(_end) if _end else q.beg
            if _d and not _d.isdigit(): q.refseq = _d.upper()
            q.altseq = _i.upper() if _i else ''
        elif _is_sub:
            q = QueryMNV()
            q.beg = parse_pos(_beg)
            q.end = parse_pos(_end) if _end else q.beg
            q.refseq = _ref.upper() if _ref else ''
            q.altseq = _alt.upper() if _alt else ''
        elif _is_dup:
            q = QueryDUP()
            q.beg = parse_pos(_beg)
            q.end = parse_pos(_end) if _end else q.beg
            q.dupseq = _dupseq.upper() if _dupseq else ''
        else:
            # use only the beg and end
            # treat input as a region
            q = Query()
            q.beg = parse_pos(_beg)
            q.end = parse_pos(_end) if _end else q.beg
            q.ref = _ref.upper() if _ref else ''
            # err_raise(InvalidInputError,
            #           'Invalid nucleotide mutation: "%s".' % mut_str, __name__)
        q.is_codon = False
    elif (mg and not muttype) or muttype == 'g':
        if not mg:
            raise InvalidInputError("Invalid genomic level input: %s" % mut_str)
        (_, _, _beg, _end_s, _end, _, _is_del, _d,
         _is_ins, _i, _is_sub, _ref, _alt, _is_dup, _dupseq) = mg.groups()
        if _is_sub and len(_ref) <= 1 and len(_alt) <= 1:
            q = QuerySNV()
            q.pos = int(_beg)
            q.ref = _ref if _ref and _ref != '?' else ''
            q.alt = _alt if _alt and _alt != '?' else ''
        elif (_is_del and _is_ins and
              (_d == '1' or len(_d) == 1) and len(_i) == 1):
            q = QuerySNV()
            q.pos = int(_beg)
            q.ref = '' if _d.isdigit() else _d.upper()
            q.alt = _i.upper() if _i else ''
        elif _is_del and not _is_ins:
            q = QueryDEL()
            q.beg = int(_beg)
            q.end = int(_end) if _end else q.beg
            q.delseq = '' if _d.isdigit() else _d.upper()
        elif _is_ins and not _is_del:
            q = QueryINS()
            q.pos = int(_beg)
            if _i:
                q.insseq = _i.upper()
            else:
                err_raise(InvalidInputError, 'insertion without inserted sequence: %s.' % mut_str)

        elif _is_ins and _is_del:
            q = QueryMNV()
            q.beg = int(_beg)
            q.end = int(_end) if _end else q.beg
            if _d and not _d.isdigit(): q.refseq = _d.upper()
            q.altseq = _i.upper() if _i else ''
        elif _is_sub:
            q = QueryMNV()
            q.beg = int(_beg)
            q.end = int(_end) if _end else q.beg
            q.refseq = _ref.upper() if _ref else ''
            q.altseq = _alt.upper() if _alt else ''
        elif _is_dup:
            q = QueryDUP()
            q.beg = int(_beg)
            q.end = int(_end) if _end else q.beg
            q.dupseq = _dupseq.upper() if _dupseq else ''
        else:
            # use only the beg and end
            # treat input as a region
            q = Query()
            q.beg = int(_beg)
            q.end = int(_end) if _end else q.beg
            q.ref = _ref.upper() if _ref else ''
            # err_raise(InvalidInputError,
            #           'Invalid nucleotide mutation: "%s".' % mut_str, __name__)
        q.is_codon = False
    else:
        err_raise(InvalidInputError, 'invalid mutation: "%s".' % mut_str)

    return q

def parse_tok_mutation_str(s, muttype=None):

    if s.find(":") >= 0:
        tok, mut_str = s.split(':', 1)
    elif s.find(".") >= 0:
        tok, mut_str = s.split('.', 1)
    else:
        try:
            tok, mut_str = s.split(None, 1)
        except:
            err_raise(InvalidInputError, 'invalid mutation string: "%s".' % s)

    q = parse_mutation_str(mut_str, muttype)
    q.tok = tok

    return q

def _list_parse_mutation(args, fields, indices, muttype=None):

    # protein
    if args.g > 0 and args.p > 0: # gene, position, ref, alt in separate columns for SNV

        q = QuerySNV()
        q.op = '|'.join(indices.extract(fields))
        q.tok = fields[args.g-1].strip()
        q.pos = parse_pos(fields[args.p-1].strip())
        if args.r > 0: q.ref = fields[args.r-1].strip()
        if args.v > 0: q.alt = fields[args.v-1].strip()
        if args.t > 0: q.tpt = fields[args.t-1].strip()
        q.is_codon = True

    elif args.g > 0 and args.n > 0:
        if muttype == 'g':      # gDNA
            q = QuerySNV()
            q.op = '|'.join(indices.extract(fields))
            q.tok = fields[args.g-1].strip()
            q.pos = int(fields[args.n-1].strip())
            if args.r > 0: q.ref = fields[args.r-1].strip()
            if args.v > 0: q.alt = fields[args.v-1].strip()
            if args.t > 0: q.tpt = fields[args.t-1].strip()
            q.is_codon = False
        else:                   # cDNA
            q = QuerySNV()
            q.op = '|'.join(indices.extract(fields))
            q.tok = fields[args.g-1].strip()
            q.pos = parse_pos(fields[args.n-1].strip())
            if args.r > 0: q.ref = fields[args.r-1].strip()
            if args.v > 0: q.alt = fields[args.v-1].strip()
            if args.t > 0: q.tpt = fields[args.t-1].strip()
            q.is_codon = False

    elif args.g > 0 and args.m > 0: # gene and mutation string

        q = parse_mutation_str(fields[args.m-1].strip(), muttype)
        q.op = '|'.join(indices.extract(fields))
        q.tok = fields[args.g-1].strip()
        if args.t > 0: q.tpt = fields[args.t-1].strip()

    elif args.m > 0:

        q = parse_tok_mutation_str(fields[args.m-1].strip(), muttype)
        q.op = '|'.join(indices.extract(fields))
        if args.t > 0: q.tpt = fields[args.t-1].strip()

    else:
        err_raise(InvalidInputError, "invalid line: %s" % line)

    return q

def list_parse_mutation(args, muttype=None):

    indices = parse_indices(args.o)
    if args.skipheader:
        args.l.readline()

    for line in args.l:
        # print line.strip()
        if args.d == 's':
            fields = line.strip().split()
        else:
            fields = line.strip().split(args.d)
        try:
            q = _list_parse_mutation(args, fields, indices, muttype)
        except IndexError as e:
            err_print(str(e))
            err_print("may be different delimiter in the input file? (try -d s)")
            sys.exit(1)
        except InvalidInputError as e:
            err_print(str(e))
            continue

        yield q, line

def parser_add_mutation(parser):

    parser.add_argument('--noheader', action='store_true', help='repress header print')
    parser.add_argument('-i', default=None,
                        help='<gene/chrm>:<mutation>, E.g., MET:1010, PIK3CA:E545K, PIK3CA:c.1633G>A, chr12:25398285')
    parser.add_argument('-l', default=None,
                        type = argparse.FileType('r'), 
                        help='mutation list file')
    parser.add_argument('-d', default="\t",
                        help="table delimiter [\\t], use 's' for space.")
    parser.add_argument('-g', type=int,
                        default=-1, help='column for gene/chromosome (1-based)')
    parser.add_argument('-p', type=int, default=-1,
                        help='column for amino acid position (1-based)')
    parser.add_argument('-n', type=int, default=-1,
                        help='column for nucleotide position (1-based)')
    parser.add_argument('-r', type=int, default=-1,
                        help='column for reference base/amino acid (1-based)')
    parser.add_argument('-v', type=int, default=-1,
                        help='column for variant base/amino acid (1-based)')
    parser.add_argument('-t', type=int, default=-1,
                        help='columns for preferred transcript (1-based)')
    parser.add_argument('-m', type=int, default=1,
                        help='column for <gene/chrm>:<mutation> (1-based)')
    parser.add_argument('-o', default='-', 
                        help='columns to be printed in output (1-based), e.g., 3,4,5-10')
    parser.add_argument('--skipheader', action='store_true',
                        help='skip header')

