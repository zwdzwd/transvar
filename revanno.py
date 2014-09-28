""" 
annotate a codon position or an amino acid change
"""
import sys, argparse, re
from transcripts import *
from utils import *
from record import *
from mutation import parser_add_mutation, parse_tok_mutation_str, list_parse_mutation
from revanno_snv import _core_annotate_nuc_snv
from revanno_del import _core_annotate_nuc_del
from revanno_ins import _core_annotate_nuc_ins
from revanno_mnv import _core_annotate_nuc_mnv

def codon_mutation(args, q):

    """ find all the mutations given a codon position, yield records """

    if q.alt and q.alt not in reverse_codon_table:
        sys.stderr.write("Unknown alternative: %s, ignore alternative.\n" % q.alt)
        q.alt = ''

    if hasattr(args, 'longest') and args.longest:
        tpts = [q.gene.longest_tpt()]
    else: tpts = q.gene.tpts

    for tpt in tpts:

        # when there's a transcript specification
        if q.tpt and tpt.name != q.tpt: continue

        if not tpt.ensure_seq(): continue
        codon = tpt.cpos2codon(q.pos)
        if not codon: continue

        # skip if reference amino acid is given
        # and codon sequence does not generate reference aa
        # codon.seq is natural sequence
        if q.ref and codon.seq not in reverse_codon_table[q.ref]: continue

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

def _core_annotate_codon(args, q):

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
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)
        
    return

def _core_annotate_nuc(args, q):

    if args.longest: tpts = [q.gene.longest_tpt()]
    else: tpts = q.gene.tpts

    if isinstance(q, QuerySNV):
        return _core_annotate_nuc_snv(args, q, tpts)
    elif isinstance(q, QueryDEL):
        return _core_annotate_nuc_del(args, q, tpts)
    elif isinstance(q, QueryINS):
        return _core_annotate_nuc_ins(args, q, tpts)
    elif isinstance(q, QueryMNV):
        return _core_annotate_nuc_mnv(args, q, tpts)

def _main_core_(args, q):

    if q.is_codon:                # in codon coordinate
        _core_annotate_codon(args, q)
    else:                       # in nucleotide coordinate
        _core_annotate_nuc(args, q)

def main_list(args, name2gene):

    for q, line in list_parse_mutation(args):
        if q.tok not in name2gene:
            r = Record()
            r.info = "GeneNotRecognized"
            r.format(q.op)
            continue
        q.gene = name2gene[q.tok]
        try:
            _main_core_(args, q)
        except UnImplementedError as e:
            err_print(line)
            raise e
        except Exception as e:
            err_print(line)
            raise e

def main_one(args, name2gene):

    q = parse_tok_mutation_str(args.i)
    if not q: return

    if q.tok not in name2gene:
        r = Record()
        r.info = 'GeneNotRecognized'
        r.format(q.op)
        return

    q.gene = name2gene[q.tok]
    q.op = args.i

    _main_core_(args, q)

def main(args):

    name2gene, thash = parse_annotation(args)

    if args.l:
        main_list(args, name2gene)

    if args.i:
        main_one(args, name2gene)

def add_parser_revanno(subparsers, d):

    parser = subparsers.add_parser("revanno", help=__doc__)
    parser_add_annotation(parser, d)
    parser_add_mutation(parser)
    parser.add_argument("--longest", action="store_true",
                        help="consider only longest transcript")
    parser.set_defaults(func=main)

