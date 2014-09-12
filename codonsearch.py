"""
search alternative codonpositions due to different transcript usage
"""
import sys, re, argparse
from mutation import parser_add_mutation, parse_tok_mutation_str, list_parse_mutation
from transcripts import *
from utils import *
from revanno import codon_mutation
from anno import pos2codon
from record import Query

outformat="{altid}\t{chrm}\t{codon1}\t{codon2}\t{tptstr}"

def _main_core_(args, q, thash):

    k2transcripts = {}
    if q.is_codon:
        for t1, c1, mutloc, tnuc, gnuc in codon_mutation(args, q):
            for cind in xrange(3):
                for t2, c2 in pos2codon(thash, t1.chrm, c1.locs[cind]):
                    if t1 == t2: continue
                    if c2.region != 'coding': continue
                    if c1.index == c2.index: continue
                    if q.ref and q.ref != standard_codon_table[c2.seq]: continue
                    altid = 'p.'
                    if q.ref: altid += q.ref
                    altid += str(c2.index)
                    k = (altid, c1.chrm, tuple(c1.locs), tuple(c2.locs))
                    tpair = '%s[%s]/%s[%s]' % (t1.name, t1.source, t2.name, t2.source)
                    if k in k2transcripts:
                        if tpair not in k2transcripts[k]:
                            k2transcripts[k].append(tpair)
                    else:
                        k2transcripts[k] = [tpair]

    for k, tpairs in k2transcripts.iteritems():
        altid, chrm, c1, c2 = k
        if q.op: s = q.op+'\t'
        else: s = ''
        s += outformat.format(altid=altid, tptstr=','.join(tpairs), chrm=chrm,
                              codon1='-'.join(map(str,c1)), codon2='-'.join(map(str,c2)))
        print s

def main(args):

    name2gene, thash = parse_annotation(args)

    if args.l:
        main_list(args, name2gene, thash)
    if args.i:
        main_one(args, name2gene, thash)

def main_list(args, name2gene, thash):

    for q in list_parse_mutation(args):

        if q.gn_name not in name2gene:
            sys.stderr.write("Gene %s is not recognized.\n" % q.gn_name)
            continue
        q.gene = name2gene[q.gn_name]

        _main_core_(args, q, thash)

def main_one(args, name2gene, thash):

    q = Query()
    ret = parse_tok_mutation_str(args.i)
    if not ret: return
    q.op = args.i
    q.gn_name, q.is_codon, q.pos, q.ref, q.alt = ret
    if q.gn_name not in name2gene:
        sys.stderr.write("Gene %s not recognized.\n" % q.gn_name)
        return
    q.gene = name2gene[q.gn_name]
    q.op = args.i

    _main_core_(args, q, thash)

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
