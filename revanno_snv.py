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

    r.taa_ref = codon2aa(codon.seq)
    r.taa_pos = codon.index
    if not q.alt:
        r.info = ''
        r.taa_alt = ''
    else:
        mut_seq = list(codon.seq[:])
        mut_seq[(q.cpos()-1) % 3] = q.alt
        r.taa_alt = codon2aa(''.join(mut_seq))
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
        raise IncompatibleTranscriptError('transcript id unmatched')
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
        except UnknownChromosomeError:
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

def codon_mutation_snv(args, q, tpt):

    """ find all the mutations given a codon position, yield records """

    if q.alt and q.alt not in reverse_codon_table:
        sys.stderr.write("Unknown alternative: %s, ignore alternative.\n" % q.alt)
        q.alt = ''

    # when there's a transcript specification
    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError('transcript id unmatched')

    tpt.ensure_seq()

    if (q.pos <= 0 or q.pos > len(tpt)):
        raise IncompatibleTranscriptError('codon nonexistent')
    codon = tpt.cpos2codon(q.pos)
    if not codon:
        raise IncompatibleTranscriptError('codon nonexistent')

    # skip if reference amino acid is given
    # and codon sequence does not generate reference aa
    # codon.seq is natural sequence
    if q.ref and codon.seq not in reverse_codon_table[q.ref]:
        raise IncompatibleTranscriptError('reference amino acid unmatched')

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name
    r.pos = '-'.join(map(str, codon.locs))
    # if alternative amino acid is given
    # filter the target mutation set to those give
    # the alternative aa
    if q.alt:
        tgt_codon_seqs = reverse_codon_table[q.alt]
        diffs = [codondiff(x, codon.seq) for x in tgt_codon_seqs]
        baseloc_list = []
        refbase_list = []
        varbase_list = []
        cdd_muts = []
        for i, diff in enumerate(diffs):
            if len(diff) == 1:
                r.tnuc_pos = (codon.index-1)*3 + 1 + diff[0]
                r.tnuc_ref = codon.seq[diff[0]]
                r.tnuc_alt = tgt_codon_seqs[i][diff[0]]

                if codon.strand == "+":
                    r.gnuc_ref = codon.seq[diff[0]]
                    r.gnuc_alt = tgt_codon_seqs[i][diff[0]]
                    r.gnuc_pos = codon.locs[diff[0]]
                    cdd_muts.append('%s:%s%d%s' % (
                        tpt.chrm, r.gnuc_ref, r.gnuc_pos, r.gnuc_alt))
                else:
                    r.gnuc_ref = complement(codon.seq[diff[0]])
                    r.gnuc_alt = complement(tgt_codon_seqs[i][diff[0]])
                    r.gnuc_pos = codon.locs[2-diff[0]]
                    cdd_muts.append('%s:%s%d%s' % (
                        tpt.chrm, r.gnuc_ref, r.gnuc_pos, r.gnuc_alt))

        r.info = "CddMuts=%s;NCodonSeq=%s;NCddSeqs=%s" % (\
            ','.join(cdd_muts), codon.seq, ','.join(tgt_codon_seqs))
    else:
        r.gnuc_range = '%d-%d' % (codon.locs[0], codon.locs[2])
        r.tnuc_range = '%d-%d' % ((codon.index-1)*3+1, (codon.index-1)*3+3)

    return r, codon

def __core_annotate_codon_snv(args, q):
    for tpt in q.gene.tpts:
        try:
            r, c = codon_mutation_snv(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
        except UnknownChromosomeError:
            continue
        yield tpt, c

def _core_annotate_codon_snv(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r, c = codon_mutation_snv(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
        except UnknownChromosomeError:
            continue
        r.taa_pos = q.pos
        r.taa_ref = q.ref
        r.taa_alt = q.alt
        r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
        r.format(q.op)
        found = True

    if not found:
        r = Record()
        r.taa_pos = q.pos
        r.taa_ref = q.ref
        r.taa_alt = q.alt
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

