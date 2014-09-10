""" 
annotate a codon position or an amino acid change
"""
import sys, argparse, re
from transcripts import *
from utils import *
from record import *
from mutation import parser_add_mutation, parse_tok_mutation_str, list_parse_mutation

def codon_mutation(args, q):

    """ find all the mutations given a codon position, yield records """

    if q.alt and q.alt not in reverse_codon_table:
        sys.stderr.write("Unknown alternative: %s, ignore alternative.\n" % q.alt)
        q.alt = ''

    if args.longest: tpts = [q.gene.longest_tpt()]
    else: tpts = q.gene.tpts

    for tpt in tpts:

        # when there's a transcript specification
        if q.tpt and tpt.name != q.tpt:
            continue

        tpt.ensure_seq()
        codon = tpt.cpos2codon(q.pos)
        if not codon: continue

        # skip if reference amino acid is given
        # and codon sequence does not generate reference aa
        if q.ref and codon.seq not in reverse_codon_table[q.ref]: # codon.seq is natural sequence
            continue

        mutloc = ''
        tnuc_pos = ''
        tnuc_ref = ''
        tnuc_alt = ''
        gnuc_ref = ''
        gnuc_alt = ''
        gnuc_pos = ''

        # if alternative amino acid is given
        # filter the target mutation set to those give the
        # alternative aa
        if q.alt:
            tgt_codon_seqs = reverse_codon_table[q.alt]
            diffs = [codondiff(x, codon.seq) for x in tgt_codon_seqs]
            baseloc_list = []
            refbase_list = []
            varbase_list = []
            cdd_muts = []
            for i, diff in enumerate(diffs):
                if len(diff) == 1:
                    tnuc_pos = (codon.index-1)*3 + 1 + diff[0]
                    tnuc_ref = codon.seq[diff[0]]
                    tnuc_alt = tgt_codon_seqs[i][diff[0]]

                    if codon.strand == "+":
                        gnuc_ref = codon.seq[diff[0]]
                        gnuc_alt = tgt_codon_seqs[i][diff[0]]
                        gnuc_pos = codon.locs[diff[0]]
                        cdd_muts.append('%s:%s%d%s' % (\
                                tpt.chrm, gnuc_ref, gnuc_pos, gnuc_alt))
                    else:
                        gnuc_ref = complement(codon.seq[diff[0]])
                        gnuc_alt = complement(tgt_codon_seqs[i][diff[0]])
                        gnuc_pos = codon.locs[2-diff[0]]
                        cdd_muts.append('%s:%s%d%s' % (\
                                tpt.chrm, gnuc_ref, gnuc_pos, gnuc_alt))

            mutloc = "CddMuts=%s;NCodonSeq=%s;NCddSeqs=%s" % (\
                ','.join(cdd_muts), codon.seq, ','.join(tgt_codon_seqs))

        yield (tpt, codon, mutloc, 
               (tnuc_pos, tnuc_ref, tnuc_alt), 
               (gnuc_pos, gnuc_ref, gnuc_alt))

def nuc_mutation(args, q):

    if args.longest: tpts = [q.gene.longest_tpt()]
    else: tpts = q.gene.tpts

    for tpt in tpts:

        tpt.ensure_seq()
        # skip if reference base is given
        if q.ref != tpt.seq[q.pos-1]:
            continue

        codon = tpt.cpos2codon((q.pos-1)/3+1)
        if not codon: continue
        if codon.strand == '+':
            gnuc_pos = codon.locs[0]+(q.pos-1)%3
        else:
            gnuc_pos = codon.locs[-1]-(q.pos-1)%3

        refaa = standard_codon_table[codon.seq]
        if not q.alt:
            mutloc = ''
            altaa = ''
        else:
            mut_seq = list(codon.seq[:])
            mut_seq[(q.pos-1) % 3] = q.alt
            altaa = standard_codon_table[''.join(mut_seq)]
            mutloc = 'CddLocs=%d;CddRefs=%s;CddAlts=%s;NCodonSeq=%s;NCddSeqs=%s' % (loc, q.ref, q.alt, codon.seq, mut_seq)

        yield tpt, codon, refaa, altaa, mutloc


def _main_core_(args, q):

    if q.is_codon:                # in codon coordinate
        found = False
        for t, codon, mutloc, tnuc, gnuc in codon_mutation(args, q):
            found = True

            r = Record()
            r.chrm = t.chrm
            r.tname = t.name
            r.reg = '%s (%s, coding)' % (t.gene.name, t.strand)
            r.pos = '-'.join(map(str, codon.locs))
            r.taa_ref = q.ref if q.ref else standard_codon_table[codon.seq]
            r.taa_alt = q.alt
            r.taa_pos = q.pos
            r.tnuc_pos, r.tnuc_ref, r.tnuc_alt = tnuc
            r.gnuc_pos, r.gnuc_ref, r.gnuc_alt = gnuc
            r.info = mutloc
            r.format(q.op)

        if not found:
            r = Record()
            r.taa_ref = q.ref
            r.taa_alt = q.alt
            r.taa_pos = q.pos
            r.info = 'status=NoValidTranscriptFound'
            r.format(q.op)

    else:                       # in nucleotide coordinate
        found = False
        for t, codon, taa_ref, taa_alt, mutloc in nuc_mutation(args, q):
            found = True

            r = Record()
            r.chrm = t.chrm
            r.tname = t.name
            r.reg = '%s (coding)' % t.gene_name
            r.pos = '-'.join(map(str, codon.locs))
            r.taa_ref = taa_ref if q.ref else standard_codon_table[codon.seq]
            r.taa_alt = taa_alt
            r.taa_pos = codon.index
            r.tnuc_pos = q.pos
            r.tnuc_ref = q.ref
            r.tnuc_alt = q.alt
            r.info = mutloc
            r.format(q.op)

        if not found:
            r = Record()
            r.tnuc_pos = q.pos
            r.tnuc_ref = q.ref
            r.tnuc_alt = q.alt
            r.info = 'status=NoValidTranscriptFound'
            r.format(q.op)

def main_list(args, name2gene):

        for q in list_parse_mutation(args):
            if q.gn_name not in name2gene:
                sys.stderr.write("Gene %s is not recognized.\n" % q.gn_name)
                continue
            q.gene = name2gene[q.gn_name]
            _main_core_(args, q)

def main_one(args, name2gene):

    q = Query()
    ret = parse_tok_mutation_str(args.i)
    if not ret: return
    q.gn_name, q.is_codon, q.pos, q.ref, q.alt = ret

    if q.gn_name not in name2gene:
        sys.stderr.write("Gene %s not recognized.\n" % q.gn_name)
        return
    q.gene = name2gene[q.gn_name]
    q.op = args.i

    _main_core_(args, q)

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
    parser.add_argument("--longest", action="store_true",
                        help="consider only longest transcript")
    parser.set_defaults(func=main)

