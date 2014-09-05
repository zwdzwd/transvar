""" 
annotate a codon position or an amino acid change
"""
import sys, argparse, re
from transcripts import *
from utils import *

def codon_mutation(args, gene, pos, ref, alt):

    if alt not in reverse_codon_table:
        return "NA\tunknown alternative: %s" % alt

    if args.alltrans:
        tpts = gene.tpts
    else:
        tpts = [gene.longest_tpt()]

    for tpt in tpts:

        codon = tpt.cpos2codon(pos)

        # skip if reference codon sequence does not
        # generate reference amino acid
        if codon.seq not in reverse_codon_table[ref]: # codon.seq is natural sequence
            continue

        tgt_codon_seqs = reverse_codon_table[alt]
        diffs = [codondiff(x, codon.seq) for x in tgt_codon_seqs]
        baseloc_list = []
        refbase_list = []
        varbase_list = []
        for i, diff in enumerate(diffs):
            if len(diff) == 1:
                if codon.strand == "+":
                    baseloc_list.append(str(codon.locs[diff[0]]))
                    refbase_list.append(codon.seq[diff[0]])
                    varbase_list.append(tgt_codon_seqs[i][diff[0]])
                else:
                    baseloc_list.append(str(codon.locs[2-diff[0]]))
                    refbase_list.append(complement(codon.seq[diff[0]]))
                    varbase_list.append(complement(tgt_codon_seqs[i][diff[0]]))
        
        mutloc = "%s\t%s\t%s" % (','.join(baseloc_list), ','.join(refbase_list), ','.join(varbase_list))
        candidates = ','.join(tgt_codon_seqs)

        return "%s\t%s\t%s\t%s=>%s\t%s\t%s" % (tpt.source, tpt.name, codon.format(), ref, alt, mutloc, candidates)

    return 'no transcripts matching reference amino acid'

def main_list(args, name2gene):

    if args.skipheader:
        args.codon_list.readline()

    for line in args.codon_list:

        fields = line.strip().split(args.d)
        if args.col_g > 0 and args.col_p > 0: # separate columns
            gn_name = fields[args.col_g-1].strip()
            posstr = fields[args.col_p-1].strip()
            if posstr.isdigit() and int(posstr) > 0:
                pos = int(posstr)
            else:
                sys.stderr.write("[Warning] abnormal position %s. skip.\n" % posstr)
                continue
            ref = fields[args.col_r-1].strip() if args.col_r > 0 else None
            alt = fields[args.col_v-1].strip() if args.col_v > 0 else None
            if gn_name not in name2gene:
                sys.stderr.write("Gene: %s not recognized.\n" % gn_name)
                continue
            gene = name2gene[gn_name]

        else:                   # <gene>:<pos> format
            m = re.match(r'([^:]*):([A-Z*]?)(\d+)([A-Z*]?)',
                         fields[args.col_gp-1])
            if not m:
                sys.stderr.write("[Warning] abnormal input %s. skip.\n" % fields[args.col_gp-1])
                continue
            gn_name = m.group(1)
            ref = m.group(2)
            pos = int(m.group(3))
            if pos <= 0:
                sys.stderr.write("[Warning] abnormal position %d. skip.\n" % pos)
                continue
            alt = m.group(4)
            if gn_name not in name2gene:
                sys.stderr.write("Gene: %s not recognized. skip.\n" % gn_name)
                continue
            gene = name2gene[gn_name]

        codon = gene.cpos2codon(pos)
        prnstr = line.strip()
        if alt:                 # with mutation
            prnstr += '\t'
            prnstr += codon_mutation(args, gene, pos, ref, alt)
        else:                   # without mutation
            prnstr += '\t'
            prnstr += codon.format()

        print prnstr

    return

def main_one(args, name2gene):

    m = re.match(r'([^:]*):([A-Z*]?)(\d+)([A-Z*]?)', args.codon)
    gn_name = m.group(1)
    ref = m.group(2)
    pos = int(m.group(3))
    alt = m.group(4)
    if gn_name not in name2gene:
        sys.stderr.write("Gene: %s not recognized.\n" % gn_name)
        return
    gene = name2gene[gn_name]
    if alt:                 # with mutation
        prnstr = args.codon
        prnstr += '\t'
        prnstr += codon_mutation(args, gene, pos, ref, alt)
    else:                   # without mutation
        prnstr = ''
        if args.alltrans:
            for tpt in gene.tpts:
                prnstr += '%s\t' % args.codon
                prnstr += tpt.cpos2codon(pos).format()
                prnstr += '\n'
        else:
            prnstr = '%s\t' % args.codon
            prnstr += gene.cpos2codon(pos).format()
            prnstr += '\n'

    sys.stdout.write(prnstr)


def main(args):

    name2gene, thash = parse_annotation(args)

    if args.codon_list:
        main_list(args, name2gene)

    if args.codon:
        main_one(args, name2gene)

def add_parser_codonanno(subparsers):

    parser = subparsers.add_parser("codonanno", help=__doc__)
    parser_add_annotation(parser)
    parser.add_argument('-c',
                        dest="codon",
                        default=None,
                        help='<gene>:[<ref>]<pos>[<alt>], E.g., MET:1010, PIK3CA:E545K')
    parser.add_argument('-l',
                        dest="codon_list",
                        default=None,
                        type = argparse.FileType('r'), 
                        help='codon list file')
    parser.add_argument('-d',
                        default="\t",
                        help="table delimiter [\\t]")
    parser.add_argument('-g', 
                        dest='col_g',
                        type=int,
                        default=-1,
                        help='column for gene (1-based)')
    parser.add_argument('-p',
                        dest='col_p',
                        type=int,
                        default=-1,
                        help='column for position (1-based)')
    parser.add_argument('-r',
                        dest='col_r',
                        type=int,
                        default=-1,
                        help='column for reference amino acid (1-based)')
    parser.add_argument('-v',
                        dest='col_v',
                        type=int,
                        default=-1,
                        help='column for variant amino acid (1-based)')
    parser.add_argument('-gp',
                        dest='col_gp',
                        type=int,
                        default=1,
                        help='column for <gene>:<ref><pos><alt> (1-based)')
    parser.add_argument("--alltrans",
                        action="store_true",
                        help="consider all transcripts")
    parser.add_argument('--skipheader',
                        action='store_true',
                        help='skip header')
    parser.set_defaults(func=main)

