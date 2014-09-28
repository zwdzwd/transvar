from transcripts import *
from utils import *
from record import *

def nuc_mutation_snv_coding(r, tpt, codon, q):

    r.reg = '%s (%s coding)' % (tpt.gene.name, tpt.strand)
    r.pos = '-'.join(map(str, codon.locs))
    r.tnuc_pos = q.cpos()
    r.tnuc_ref = q.ref
    r.tnuc_alt = q.alt

    if (q.ref and q.ref != tpt.seq[q.cpos()-1]):
        raise IncompatibleTranscriptError('SNV ref not matched')
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

def nuc_mutation_snv_intronic(r, tpt, codon, q):

    r.reg = '%s (%s intronic)' % (tpt.gene.name, tpt.strand)
    i = q.cpos() - (codon.index-1)*3 - 1
    np = tpt.position_array()
    check_exon_boundary(np, q.pos)
    r.gnuc_pos = tnuc2gnuc2(np, q.pos, tpt)
    r.gnuc_ref = faidx.refgenome.fetch_sequence(tpt.chrm, r.gnuc_pos, r.gnuc_pos)
    r.pos = r.gnuc_pos
    r.gnuc_ref = faidx.refgenome.fetch_sequence(tpt.chrm, r.gnuc_pos, r.gnuc_pos)
    if tpt.strand == '+':
        if q.ref and r.gnuc_ref != q.ref: raise IncompatibleTranscriptError()
        r.gnuc_alt = q.alt if q.alt else ''
    else:
        if q.ref and r.gnuc_ref != complement(q.ref): raise IncompatibleTranscriptError()
        r.gnuc_alt = complement(q.alt) if q.alt else ''
    r.tnuc_pos = q.pos
    r.tnuc_ref = r.gnuc_ref if tpt.strand == '+' else complement(r.gnuc_ref)
    r.tnuc_alt = q.alt

def nuc_mutation_snv(args, q, tpt):

    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError()
    tpt.ensure_seq()

    if (q.cpos() <= 0 or q.cpos() > len(tpt)):
        raise IncompatibleTranscriptError()
    codon = tpt.cpos2codon((q.cpos()+2)/3)
    if not codon:
        raise IncompatibleTranscriptError()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name

    if q.pos.tpos == 0:                # coding region
        nuc_mutation_snv_coding(r, tpt, codon, q)
    else:          # coordinates are with respect to the exon boundary
        nuc_mutation_snv_intronic(r, tpt, codon, q)

    return r

def _core_annotate_nuc_snv(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = nuc_mutation_snv(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
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
