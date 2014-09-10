import re, sys, argparse
from utils import *
from record import Query

def parse_mutation_str(mut_str):

    mp = re.match(r'(p)?(\.)?([A-Z*]?)(\d+)([A-Z*]?)$', mut_str)
    mn = re.match(r'(c)?(\.)?(\d+)(\.)?([ATGC]?)>([ATGC]?)$', mut_str)

    if mn and not mp:
        is_codon = False
        pos = int(mn.group(3))
        ref = mn.group(5) if mn.group(5) else ''
        alt = mn.group(6) if mn.group(6) else ''
    elif mp:
        is_codon = True
        ref = mp.group(3) if mp.group(3) else ''
        pos = int(mp.group(4))
        alt = mp.group(5) if mp.group(5) else ''
    else:
        sys.stderr.write('Cannot infer mutation type "%s", skip.\n' % mut_str)
        return None

    return (is_codon, pos, ref, alt)

def parse_tok_mutation_str(s):

    if s.find(":") >= 0:
        tok, mut_str = s.split(':', 1)
    elif s.find(".") >= 0:
        tok, mut_str = s.split('.', 1)
    else:
        try:
            tok, mut_str = s.split(None, 1)
        except:
            sys.stderr.write('[%s] unacceptable mutation string: "%s", skip.\n' % 
                             (__name__, s))
            return None

    ret = parse_mutation_str(mut_str)
    if ret:
        is_codon, pos, ref, alt = ret
        return (tok, is_codon, pos, ref, alt)
    else:
        return None


def list_parse_mutation(args):

    indices = parse_indices(args.o)
    if args.skipheader:
        args.l.readline()

    for line in args.l:
        fields = line.strip().split(args.d)
        q = Query()
        q.op = '\t'.join(indices.extract(fields))
        if args.g > 0 and args.p > 0: # gene, position, ref, alt in separate columns
            q.gn_name = fields[args.g-1].strip()
            if not q.set_pos(fields[args.p-1].strip()): continue
            if args.r > 0: q.ref = fields[args.r-1].strip()
            if args.v > 0: q.alt = fields[args.v-1].strip()
            if args.t > 0: q.tpt = fields[args.t-1].strip()
            q.is_codon = True
            yield q

        elif args.g > 0 and args.n > 0:
            q.gn_name = fields[args.g-1].strip()
            if not q.set_pos(fields[args.n-1].strip()): continue
            if args.r > 0: q.ref = fields[args.r-1].strip()
            if args.v > 0: q.alt = fields[args.v-1].strip()
            if args.t > 0: q.tpt = fields[args.t-1].strip()
            q.is_codon = False
            yield q

        elif args.g > 0 and args.m > 0:
            q.gn_name = fields[args.g-1].strip()
            if args.t > 0: q.tpt = fields[args.t-1].strip()
            ret = parse_mutation_str(fields[args.m-1].strip())
            if ret:
                q.is_codon, q.pos, q.ref, q.alt = ret
                yield q

        elif args.m > 0:

            ret = parse_tok_mutation_str(fields[args.m-1].strip())
            if ret:
                q.gn_name, q.is_codon, q.pos, q.ref, q.alt = ret
                if args.t > 0: q.tpt = fields[args.t-1].strip()
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

