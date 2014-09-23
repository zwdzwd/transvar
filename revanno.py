""" 
annotate a codon position or an amino acid change
"""
import sys, argparse, re
from transcripts import *
from utils import *
from record import *
from mutation import parser_add_mutation, parse_tok_mutation_str, list_parse_mutation
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

def nuc_mutation_snv_coding(r, tpt, codon, q):

    r.reg = '%s (%s coding)' % (tpt.gene.name, tpt.strand)
    r.pos = '-'.join(map(str, codon.locs))
    r.tnuc_pos = q.cpos()
    r.tnuc_ref = q.ref
    r.tnuc_alt = q.alt

    if (q.ref and q.ref != tpt.seq[q.cpos()-1]): return False
    if codon.strand == '+':
        r.gnuc_pos = codon.locs[0]+(q.cpos()-1)%3
        r.gnuc_ref = q.ref if q.ref else codon.seq[(q.cpos()-1)%3]
        if q.alt: r.gnuc_alt = q.alt
    else:
        r.gnuc_pos = codon.locs[-1]-(q.cpos()-1)%3
        r.gnuc_ref = complement(q.ref if q.ref else codon.seq[(q.cpos()-1)%3])
        if q.alt: r.gnuc_alt = complement(q.alt)

    r.taa_ref = standard_codon_table[codon.seq]
    r.taa_pos = codon.index
    if not q.alt:
        r.info = ''
        r.taa_alt = ''
    else:
        mut_seq = list(codon.seq[:])
        mut_seq[(q.cpos()-1) % 3] = q.alt
        r.taa_alt = standard_codon_table[''.join(mut_seq)]
        r.info = 'NCodonSeq=%s;NAltCodonSeq=%s' % (codon.seq, ''.join(mut_seq))
        
    return True

def nuc_mutation_snv_intronic(r, tpt, codon, q):

    r.reg = '%s (%s intronic)' % (tpt.gene.name, tpt.strand)
    i = q.cpos() - (codon.index-1)*3 - 1
    if q.pos.tpos > 0:
        if tpt.strand == '+':
            if codon.locs[i]+1 == codon.locs[i+1]: return False
            r.gnuc_pos = codon.locs[i] + q.pos.tpos
            r.pos = '%s-(%d)-%s' % ('-'.join(map(str, codon.locs[:i+1])), r.gnuc_pos,
                                    '-'.join(map(str, codon.locs[i+1:])))
        elif tpt.strand == '-':
            ir = 2-i
            if codon.locs[ir-1]+1 == codon.locs[ir]: return False
            r.gnuc_pos = codon.locs[ir] - q.pos.tpos
            r.pos = '%s-(%d)-%s' % ('-'.join(map(str, codon.locs[:ir])), r.gnuc_pos,
                                    '-'.join(map(str, codon.locs[ir:])))
    elif q.pos.tpos < 0:
        if tpt.strand == '+':
            if codon.locs[i-1]+1 == codon.locs[i]: return False
            r.gnuc_pos = codon.locs[i] + q.pos.tpos
            r.pos = '%s-(%d)-%s' % ('-'.join(map(str, codon.locs[:i])), r.gnuc_pos,
                                    '-'.join(map(str, codon.locs[i:])))
        elif tpt.strand == '-':
            ir = 2-i
            if codon.locs[ir]+1 == codon.locs[ir+1]: return False
            r.gnuc_pos = codon.locs[ir] - q.pos.tpos
            r.pos = '%s-(%d)-%s' % ('-'.join(map(str, codon.locs[:ir+1])), r.gnuc_pos,
                                    '-'.join(map(str, codon.locs[ir+1:])))
    else:
        raise Exception('Conflicting region: coding vs intronic')

    r.gnuc_ref = faidx.refgenome.fetch_sequence(tpt.chrm, r.gnuc_pos, r.gnuc_pos)
    if tpt.strand == '+':
        if q.ref and r.gnuc_ref != q.ref: return False
        r.gnuc_alt = q.alt if q.alt else ''
    else:
        if q.ref and r.gnuc_ref != complement(q.ref): return False
        r.gnuc_alt = complement(q.alt) if q.alt else ''
    r.tnuc_pos = q.pos
    r.tnuc_ref = r.gnuc_ref if tpt.strand == '+' else complement(r.gnuc_ref)
    r.tnuc_alt = q.alt

    return True

def nuc_mutation_snv(args, q, tpt):

    if q.tpt and tpt.name != q.tpt: return None
    tpt.ensure_seq()

    if (q.cpos() <= 0 or q.cpos() > len(tpt)): return None
    codon = tpt.cpos2codon((q.cpos()+2)/3)
    if not codon: return None

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name

    if q.pos.tpos == 0:                # coding region
        if not nuc_mutation_snv_coding(r, tpt, codon, q): return None
    else:          # coordinates are with respect to the exon boundary
        if not nuc_mutation_snv_intronic(r, tpt, codon, q): return None

    return r

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
        r.info = 'status=NoValidTranscriptFound'
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

def _core_annotate_nuc_snv(args, q, tpts):

    found = False
    for tpt in tpts:
        r = nuc_mutation_snv(args, q, tpt)
        if r:
            found = True
            r.format(q.op)

    if not found:
        r = Record()
        r.tnuc_pos = q.pos
        r.tnuc_ref = q.ref
        r.tnuc_alt = q.alt
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

    return


def _main_core_(args, q):

    if q.is_codon:                # in codon coordinate
        _core_annotate_codon(args, q)
    else:                       # in nucleotide coordinate
        _core_annotate_nuc(args, q)

def main_list(args, name2gene):

    for q in list_parse_mutation(args):
        if q.tok not in name2gene:
            sys.stderr.write("Gene %s is not recognized.\n" % q.tok)
            continue
        q.gene = name2gene[q.tok]
        _main_core_(args, q)

def main_one(args, name2gene):

    q = parse_tok_mutation_str(args.i)
    if not q: return

    if q.tok not in name2gene:
        sys.stderr.write("Gene %s not recognized.\n" % q.tok)
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

