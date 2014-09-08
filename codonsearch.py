"""
search alternative codonpositions due to different transcript usage
"""
import sys, re, argparse
from mutation import parser_add_mutation, parse_mutation_str, list_parse_mutation
from transcripts import *
from utils import *
from revanno import codon_mutation
from anno import pos2codon

outformat="{altid}\t{tptstr}"

def main(args):

    name2gene, thash = parse_annotation(args)

    for op, is_codon, gn_name, pos, ref, alt in list_parse_mutation(args):
        if gn_name not in name2gene:
            sys.stderr.write("Gene %s is not recognized.\n" % gn_name)
            continue
        gene = name2gene[gn_name]

        altid2transcripts = {}
        if is_codon:
            for t1, c1, mutloc in codon_mutation(args, gene, pos, ref, alt):
                for t2, c2 in pos2codon(thash, t1.chrm, c1.locs[0]):
                    if t1 == t2: continue
                    if c2.region == 'coding' and c2.index != c1.index:
                        if c2.index in altid2transcripts:
                            altid2transcripts[c2.index] += ',%s[%s]/%s[%s]' % (t1.name, t2.source, t2.name, t2.source)
                        else:
                            altid2transcripts[c2.index] = '%s[%s]/%s[%s]' % (t1.name, t2.source, t2.name, t2.source)
        for altid, tptstr in altid2transcripts.iteritems():
            if op: s = op+'\t'
            else: s = ''
            s += outformat.format(altid=altid, tptstr=tptstr)
            print s

def oldmain():

    for line in args.l:
        fields = line.strip().split(args.d)
        op = '\t'.join(indices.extract(fields))

        for line in f:

            # if line.strip() != "TP53:281":
                # continue

            gene = name2gene[gene_name]
            codon_pos = int(codon_pos)
            loc2trans = {}
            locs = set()
            for i, trans in enumerate(gene.transcripts):
                nucpos = trans.aa_pos2nuc_pos(codon_pos)
                if nucpos:
                    locs.add(tuple(nucpos))

                    if tuple(nucpos) in loc2trans:
                        loc2trans[tuple(nucpos)].append(trans)
                    else:
                        loc2trans[tuple(nucpos)] = [trans]

                    # print trans, trans.cds_beg, trans.cds_end, nucpos, trans.exons

            aa_poses = set()
            aa_pos2trans = {}
            for loc in locs:
                for i, trans in enumerate(gene.transcripts):
                    result = trans.nuc_pos2aa_pos(loc[0])
                    if not result:
                        continue
                    aa_pos, codon_loc = result
                    # print i, trans, "nucleotide loc: ", loc, "aa loc: ", result

                    # if aa_pos <=0:
                    #     print trans, trans.cds_beg, trans.cds_end, loc, aa_pos
                    #     print trans.exons
                    #     import sys
                    #     sys.exit(1)

                    if aa_pos != codon_pos:
                        # aa_poses.add("%s:%d:%d" % (gene_name, aa_pos, codon_loc))
                        aa_poses.add("%s.p%d" % (gene_name, aa_pos))

                        aa_pos_str = "%s.p%d" % (gene_name, aa_pos)
                        if aa_pos_str in aa_pos2trans:
                            aa_pos2trans[aa_pos_str].append((trans))
                        else:
                            aa_pos2trans[aa_pos_str] = [(aa_pos2trans, trans)]


            print "%s.p%s\t%s\t%d\t%s" % (gene_name, codon_pos, gene.strand(), len(gene.transcripts), '\t'.join(aa_poses))


def add_parser_codonsearch(subparsers):

    parser = subparsers.add_parser('codonsearch', help=__doc__)
    parser_add_mutation(parser)
    parser_add_annotation(parser)
    parser.set_defaults(func=main)
