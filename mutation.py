import re, sys, argparse
from utils import *
from record import *

def parse_mutation_str(mut_str):
    mut_str = mut_str.strip()
    mp = re.match(r'(p)?(\.)?([A-Z*?]?)(\d+)([A-Z*?]?)$', mut_str)
    # mn = re.match(r'(c)?(\.)?(\d+)(\.)?([ATGC?]?)>([ATGC?]?)$', mut_str)
    mn = re.match(r'(c)?(\.)?([\d+-]+)(_([\d+-]+))?(\.)?(del([atgcATGC\d]*))?(ins([atgcATGC]*))?(([atgcATGC?]*)>([atgcATGC?]*))?$', mut_str)

    if mp:
        is_codon = True
        ref = mp.group(3) if mp.group(3) and mp.group(3) != '?' else ''
        pos = int(mp.group(4))
        alt = mp.group(5) if mp.group(5) and mp.group(5) != '?' else ''
    elif mn:
        (_, _, _beg, _end_s, _end, _, _is_del, _d,
         _is_ins, _i, _is_sub, _ref, _alt) = mn.groups()
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
            else: err_die('Insertion without inserted sequence: %s.' % mut_str, __name__)
        elif _is_ins and _is_del:
            q = QueryMNV()
            q.beg = parse_pos(_beg)
            q.end = parse_pos(_end) if _end else q.beg
            q.refseq = _ref.upper() if _ref else ''
            q.altseq = _alt.upper() if _alt else ''
        elif _is_sub:
            q = QueryMNV()
            q.beg = parse_pos(_beg)
            q.end = parse_pos(_end) if _end else q.beg
            if _d and not _d.isdigit(): q.refseq = _d.upper()
            q.altseq = _i.upper() if _i else ''
        else:
            err_raise(InvalidInputError,
                      'Invalid nucleotide mutation: "%s".' % mut_str, __name__)
    else:
        err_raise(InvalidInputError,
                  'Invalid mutation: "%s".' % mut_str, __name__)

    q.is_codon = False
    return q

def parse_tok_mutation_str(s):

    if s.find(":") >= 0:
        tok, mut_str = s.split(':', 1)
    elif s.find(".") >= 0:
        tok, mut_str = s.split('.', 1)
    else:
        try:
            tok, mut_str = s.split(None, 1)
        except:
            err_raise(InvalidInputError,
                      'Invalid mutation string: "%s".' % s, __name__)

    q = parse_mutation_str(mut_str)
    q.tok = tok

    return q

def _list_parse_mutation(args, fields, indices):

    if args.g > 0 and args.p > 0: # gene, position, ref, alt in separate columns for SNV

        q = QuerySNV()
        q.op = '\t'.join(indices.extract(fields))
        q.tok = fields[args.g-1].strip()
        q.pos = parse_pos(fields[args.p-1].strip())
        if args.r > 0: q.ref = fields[args.r-1].strip()
        if args.v > 0: q.alt = fields[args.v-1].strip()
        if args.t > 0: q.tpt = fields[args.t-1].strip()
        q.is_codon = True

    elif args.g > 0 and args.n > 0:

        q = QuerySNV()
        q.op = '\t'.join(indices.extract(fields))
        q.tok = fields[args.g-1].strip()
        q.pos = parse_pos(fields[args.n-1].strip())
        if args.r > 0: q.ref = fields[args.r-1].strip()
        if args.v > 0: q.alt = fields[args.v-1].strip()
        if args.t > 0: q.tpt = fields[args.t-1].strip()
        q.is_codon = False

    elif args.g > 0 and args.m > 0:

        q = parse_mutation_str(fields[args.m-1].strip())
        q.op = '\t'.join(indices.extract(fields))
        q.tok = fields[args.g-1].strip()
        if args.t > 0: q.tpt = fields[args.t-1].strip()

    elif args.m > 0:

        q = parse_tok_mutation_str(fields[args.m-1].strip())
        q.op = '\t'.join(indices.extract(fields))
        if args.t > 0: q.tpt = fields[args.t-1].strip()

    else:
        err_raise(InvalidInputError, "Invalid line: %s" % line, __name__)

    return q

def list_parse_mutation(args):

    indices = parse_indices(args.o)
    if args.skipheader:
        args.l.readline()

    for line in args.l:
        fields = line.strip().split(args.d)
        try:
            q = _list_parse_mutation(args, fields, indices)
        except InvalidInputError as e:
            err_print(str(e))
            continue

        print line
        print type(q), q.tok
        yield q

def parser_add_mutation(parser):

    parser.add_argument('-i', default=None,
                        help='<gene>:<mutation>, E.g., MET:1010, PIK3CA:E545K, PIK3CA:c.1633G>A')
    parser.add_argument('-l', default=None,
                        type = argparse.FileType('r'), 
                        help='mutation list file')
    parser.add_argument('-d', default="\t",
                        help="table delimiter [\\t]")
    parser.add_argument('-g', type=int,
                        default=-1, help='column for gene (1-based)')
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
                        help='column for <gene>:<mutation> (1-based)')
    parser.add_argument('-o', default='-', 
                        help='columns to be printed in output (1-based), e.g., 3,4,5-10')
    parser.add_argument('--skipheader', action='store_true',
                        help='skip header')

