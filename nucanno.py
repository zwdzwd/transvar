"""
annotate nucleotide position(s) or mutations
"""
import sys, argparse, re
from transcripts import *
from utils import parse_indices

def nuc_mutation(codon, pos, ref, alt):

    ref_aa = standard_codon_table[codon.seq]
    if codon.strand == "+":
        alt_seq = set_seq(codon.seq, pos-codon.locs[0], alt)
    else:
        alt_seq = set_seq(codon.seq, 2-(pos-codon.locs[0]), complement(alt))
    alt_aa = standard_codon_table[alt_seq]

    return "%s\t%s=>%s\t%s" % (codon.format(), ref_aa, alt_aa, alt_seq)


def main_list(args, thash):

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

    name2gene, thash = parse_annotation(args.annotation_file)

    if args.npos_list:
        main_list(args, thash)

    if args.npos:
        main_one(args, thash)
    
def add_parser_nucanno(subparsers):

    parser = subparsers.add_parser('nucanno', help=__doc__)
    parser.add_argument('-a',
                        metavar='annotation',
                        required = True,
                        dest='annotation_file', 
                        help='protein annotation file')
    parser.add_argument('-n',
                        dest="npos",
                        default=None,
                        help="<chrm>:[<ref>]<pos>[<alt>], E.g., chr12:25398285")
    parser.add_argument('-l',
                        dest="npos_list",
                        type = argparse.FileType('r'),
                        default=None,
                        help='nucleotide position table')
    parser.add_argument('-d',
                        default="\t", 
                        help="table delimiter [\\t]")
    parser.add_argument('-c',
                        dest='col_c',
                        type=int,
                        default=-1,
                        help='column for chromosome (1-based)')
    parser.add_argument("-p",
                        dest='col_p',
                        type=int,
                        default=-1,
                        help='column for position (1-based)')
    parser.add_argument('-r',
                        dest='col_r',
                        type=int,
                        default=-1,
                        help='column for reference base (1-based)')
    parser.add_argument('-v',
                        dest='col_v',
                        type=int,
                        default=-1,
                        help='column for variant base (1-based)')
    parser.add_argument("-cp",
                        dest="col_cp",
                        type=int,
                        default=-1,
                        help="column for <chrm>:[<ref>]<pos>[<alt>] (1-based)")
    parser.add_argument('-oc',
                        dest='outcol',
                        default='-',
                        help='columns in the original table to be output (1-based) [everything]')
    parser.add_argument('--standard', action="store_true")
    parser.add_argument('--skipheader',
                        action='store_true',
                        help='skip header')
    parser.add_argument('--alltrans', action='store_true', help='report results from all transcripts')

    parser.set_defaults(func=main)
