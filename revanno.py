""" 
annotate a codon position or an amino acid change
"""
import sys, argparse, re
from transcripts import *
from utils import *
from mutation import parser_add_mutation, parse_mutation_str, list_parse_mutation, mtemplate

def codon_mutation(args, gene, pos, ref, alt, longest_only=False):

    if alt != '.' and alt not in reverse_codon_table:
        sys.stderr.write("Unknown alternative: %s, ignore alternative.\n" % alt)
        alt = '.'

    if longest_only: tpts = [gene.longest_tpt()]
    else: tpts = gene.tpts

    for tpt in tpts:

        tpt.ensure_seq()
        codon = tpt.cpos2codon(pos)

        # skip if reference amino acid is given
        # and codon sequence does not generate reference aa
        if ref != '.' and codon.seq not in reverse_codon_table[ref]: # codon.seq is natural sequence
            continue

        # if alternative amino acid is given
        # filter the target mutation set to those give the
        # alternative aa
        if alt == '.':          # no alternative amino acid
            mutloc = '.\t.\t.\t.'
        else:
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

            mutloc = "%s\t%s\t%s\t%s" % (','.join(baseloc_list) if baseloc_list else '.',
                                         ','.join(refbase_list) if refbase_list else '.',
                                         ','.join(varbase_list) if varbase_list else '.',
                                         ','.join(tgt_codon_seqs))

        yield tpt, codon, mutloc

def nuc_mutation(args, gene, pos, ref, alt):

    if args.alltrans:
        tpts = gene.tpts
    else:
        tpts = [gene.longest_tpt()]

    for tpt in tpts:

        tpt.ensure_seq()
        # skip if reference base is given
        if ref != tpt.seq[pos-1]:
            continue

        codon = tpt.cpos2codon((pos-1)/3+1)
        if codon.strand == '+':
            loc = codon.locs[0]+(pos-1)%3
        else:
            loc = codon.locs[-1]-(pos-1)%3
        refaa = standard_codon_table[codon.seq]
        if alt == '.':
            mutloc = '.\t.\t.\t.'
            altaa = '.'
        else:
            mutloc = '%d\t%s\t%s' % (loc, ref, alt)
            mut_seq = list(codon.seq[:])
            mut_seq[(pos-1) % 3] = alt
            altaa = standard_codon_table[''.join(mut_seq)]

        yield tpt, codon, refaa, altaa, mutloc


def main_list(args, name2gene):

        for op, is_codon, gn_name, pos, ref, alt in list_parse_mutation(args):
            if gn_name not in name2gene:
                sys.stderr.write("Gene %s is not recognized.\n" % gn_name)
                continue
            gene = name2gene[gn_name]

            _main_core_(args, op, is_codon, gene, pos, ref, alt)

    #     if args.col_g > 0 and args.col_p > 0: # separate columns
    #         gn_name = fields[args.col_g-1].strip()
    #         posstr = fields[args.col_p-1].strip()
    #         if posstr.isdigit() and int(posstr) > 0:
    #             pos = int(posstr)
    #         else:
    #             sys.stderr.write("[Warning] abnormal position %s. skip.\n" % posstr)
    #             continue
    #         ref = fields[args.col_r-1].strip() if args.col_r > 0 else None
    #         alt = fields[args.col_v-1].strip() if args.col_v > 0 else None
    #         if gn_name not in name2gene:
    #             sys.stderr.write("Gene: %s not recognized.\n" % gn_name)
    #             continue
    #         gene = name2gene[gn_name]

    #     else:                   # <gene>:<pos> format
    #         m = re.match(r'([^:]*):([A-Z*]?)(\d+)([A-Z*]?)',
    #                      fields[args.col_gp-1])
    #         if not m:
    #             sys.stderr.write("[Warning] abnormal input %s. skip.\n" % fields[args.col_gp-1])
    #             continue
    #         gn_name = m.group(1)
    #         ref = m.group(2)
    #         pos = int(m.group(3))
    #         if pos <= 0:
    #             sys.stderr.write("[Warning] abnormal position %d. skip.\n" % pos)
    #             continue
    #         alt = m.group(4)
    #         if gn_name not in name2gene:
    #             sys.stderr.write("Gene: %s not recognized. skip.\n" % gn_name)
    #             continue
    #         gene = name2gene[gn_name]

    #     codon = gene.cpos2codon(pos)
    #     prnstr = line.strip()
    #     if alt:                 # with mutation
    #         prnstr += '\t'
    #         prnstr += codon_mutation(args, gene, pos, ref, alt)
    #     else:                   # without mutation
    #         prnstr += '\t'
    #         prnstr += codon.format()

    #     print prnstr

    # return

def _main_core_(args, op, is_codon, gene, pos, ref, alt):

    found = False
    if is_codon:                # in codon coordinate
        for t, codon, mutloc in codon_mutation(args, gene, pos, ref, alt, args.longest_only):
            found = True
            s = op+'\t' if op else ''
            s += mtemplate.format(t=t, c=codon.format(), mutloc=mutloc, alt=alt,
                                  ref=ref if ref!='.' else standard_codon_table[codon.seq])
            print s
    else:                       # in nucleotide coordinate
        for t, codon, refaa, altaa, mutloc in nuc_mutation(args, gene, pos, ref, alt):
            found = True
            s = op+'\t' if op else ''
            s += mtemplate.format(t=t, c=codon.format(), mutloc=mutloc, ref=refaa, alt=altaa)
            print s

    if not found:
        print notemplate.format(gene=gene, pos=pos, ref=ref, alt=alt)

def main_one(args, name2gene):

    gn_name, mut_str = args.i.split(':', 1)
    is_codon, pos, ref, alt = parse_mutation_str(mut_str)
    if not pos: return

    if gn_name not in name2gene:
        sys.stderr.write("Gene %s not recognized.\n" % gn_name)
        return
    gene = name2gene[gn_name]

    _main_core_(args, args.i, is_codon, gene, pos, ref, alt)

def main(args):

    name2gene, thash = parse_annotation(args)

    if args.l:
        main_list(args, name2gene)

    if args.i:
        main_one(args, name2gene)

def add_parser_revanno(subparsers):

    parser = subparsers.add_parser("revanno", help=__doc__)
    parser_add_annotation(parser)
    parser_add_mutation(parser)
    parser.add_argument("--longest_only", action="store_true",
                        help="consider only longest transcript")
    parser.set_defaults(func=main)

