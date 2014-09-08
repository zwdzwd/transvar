"""
annotate nucleotide position(s) or mutations
"""
import sys, argparse, re
from transcripts import *
from utils import *
from mutation import parse_mutation_str, mtemplate

def pos2codon(thash, chrm, pos):

    for t in thash.get_transcripts(chrm, pos):
        c = t.npos2codon(chrm, pos)
        yield t, c

def _main_core_(args, thash, op, chrm, pos, ref, alt):

    for t, c in pos2codon(thash, chrm, pos):

        if isinstance(c, Codon):

            ref_aa = standard_codon_table[c.seq]
            if alt != '.':
                if c.strand == "+":
                    alt_seq = set_seq(c.seq, pos-c.locs[0], alt)
                else:
                    alt_seq = set_seq(c.seq, 2-(pos-c.locs[0]), complement(alt))
                alt_aa = standard_codon_table[alt_seq]
            else:
                alt_aa = '.'

            s = op+'\t' if op else ''
            s += mtemplate.format(t=t, c=c.format(), ref=ref_aa, alt=alt_aa,
                                  mutloc='%d\t%s\t%s' % (pos, ref, alt))
        elif isinstance(c, NonCoding):
            s += mtemplate.format(t=t, c=c.format(), ref='.', alt='.',
                                  mutloc='%d\t%s\t%s' % (pos, ref, alt))
        
        print s

def list_parse_genomic_mutation(args):

    indices = parse_indices(args.o)
    if args.skipheader:
        args.l.readline()

    for line in args.l:
        fields = line.strip().split(args.d)
        op = '\t'.join(indices.extract(fields))
        if args.c > 0 and args.p > 0:
            chrm = fields[args.c-1].strip()
            pos = int(fields[args.p-1].strip())
            yield op, chrm, pos, ref, alt
        elif args.m > 0:
            ret = parse_mutation_str(fields[args.m-1])
            if ret: 
                chrm, is_codon, pos, ref, alt = ret
                yield op, chrm, pos, ref, alt

def main_list(args, thash):

    for op, chrm, pos, ref, alt in list_parse_genomic_mutation(args):
        _main_core_(args, thash, op, chrm, pos, ref, alt)

def oldmain():
    if args.skipheader:
        args.npos_list.readline()

    outindices = parse_indices(args.outcol)
    for line in args.npos_list:

        fields = line.strip().split(args.d)
        if args.col_c > 0 and args.col_p > 0: # separate columns
            chrm = fields[args.col_c-1]
            pos = int(fields[args.col_p-1])
            ref = fields[args.col_r-1].strip() if args.col_r > 0 else None
            alt = fields[args.col_v-1].strip() if args.col_v > 0 else None
        else:                   # <chrm>:<pos> format
            m = re.match(r'([^:]*):([ATGC]?)(\d+)([ATGC]?)',
                         fields[args.col_cp-1])
            chrm = m.group(1)
            ref = m.group(2)
            pos = int(m.group(3))
            alt = m.group(4)

        tpts = thash.get_transcripts(chrm, pos, args.standard)
        if not tpts:
            continue

        if not args.alltrans:
            tpt = tpts[0]
            if len(tpts) > 1:
                for _tpt in tpts:
                    if _tpt.is_standard():
                        tpt = _tpt
                        break
            tpts = [tpt]

        for tpt in tpts:
            prncol = outindices.extract(fields)
            c = tpt.npos2codon(chrm, pos)
            if isinstance(c, Codon):
                if alt:
                    prncol.append(nuc_mutation(c, pos, ref, alt))
                else:
                    prncol.append(c.format())
            elif isinstance(c, NonCoding):
                prncol.append(c.format())

            print '\t'.join(prncol)


def main_one(args, thash):

    ret = parse_mutation_str(args.i)
    if ret:
        chrm, is_codon, pos, ref, alt = ret
        _main_core_(args, thash, args.i, chrm, pos, ref, alt)

def oldmainone():
    m = re.match(r'([^:]*):([ATGC]?)(\d+)([ATGC]?)', args.npos)
    chrm = m.group(1)
    ref = m.group(2)
    pos = int(m.group(3))
    alt = m.group(4)

    tpts = thash.get_transcripts(chrm, pos, args.standard)
    if not tpts:
        return

    if not args.alltrans:
        tpt = tpts[0]
        if len(tpts) > 1:
            for _tpt in tpts:
                if _tpt.is_standard():
                    tpt = _tpt
                    break
        tpts = [tpt]

    for tpt in tpts:
        prnstr = args.npos
        c = tpt.npos2codon(chrm, pos)
        if isinstance(c, Codon):
            if alt:
                prnstr += '\t'
                prnstr += nuc_mutation(c, pos, ref, alt)
            else:
                prnstr += '\t'
                prnstr += c.format()
        elif isinstance(c, NonCoding):
            prnstr += '\t'+c.format()
            if alt:
                prnstr += "\t%s\t%s" % (ref, alt)

        print prnstr


def main(args):

    name2gene, thash = parse_annotation(args)

    if args.l:
        main_list(args, thash)

    if args.i:
        main_one(args, thash)
    
def add_parser_anno(subparsers):

    parser = subparsers.add_parser('anno', help=__doc__)
    parser_add_annotation(parser)
    parser.add_argument('-i', default=None,
                        help="<chrm>:[<ref>]<pos>[<alt>], E.g., chr12:25398285")
    parser.add_argument('-l', default=None,
                        type = argparse.FileType('r'),
                        help='mutation list file')
    parser.add_argument('-d', default="\t", 
                        help="table delimiter of mutation list [\\t]")
    parser.add_argument('-c', type=int, default=-1,
                        help='column for chromosome (1-based)')
    parser.add_argument("-p", type=int, default=-1,
                        help='column for position (1-based)')
    parser.add_argument('-r', type=int, default=-1,
                        help='column for reference base (1-based)')
    parser.add_argument('-v', type=int, default=-1,
                        help='column for variant base (1-based)')
    parser.add_argument("-m", type=int, default=-1,
                        help="column for mutation in format <chrm>:[<ref>]<pos>[<alt>] (1-based)")
    parser.add_argument('-o', default='-',
                        help='columns in the original table to be output (1-based)')
    parser.add_argument('--skipheader', action='store_true',
                        help='skip header')
    parser.add_argument('--longest', action="store_true", help='use longest transcript instead of reporting all transcripts')

    parser.set_defaults(func=main)
