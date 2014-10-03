from transcripts import *
from utils import *
from record import *
from revanno_ins import nuc_mutation_ins_coding
from copy import copy

def nuc_mutation_dup(args, q, tpt):
    
    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    tpt.ensure_seq()

    r = Record('dup')
    r.chrm = tpt.chrm
    r.tname = tpt.name

    if q.beg.pos > len(tpt) or q.end.pos > len(tpt):
        raise IncompatibleTranscriptError('codon nonexistent')

    np = tpt.position_array()
    check_exon_boundary(np, q.beg)
    check_exon_boundary(np, q.end)

    # identify duplicate seq
    if q.beg.tpos == 0 and q.end.tpos == 0:
        natdupseq = tpt.seq[q.beg.pos-1:q.end.pos]
        refdupseq = reverse_complement(natdupseq) if tpt.strand == '-' else natdupseq
        r.gnuc_beg, r.gnuc_end = tnuc_range2gnuc_range_(np, q.beg.pos, q.end.pos)
    else:
        _gnuc_beg = tnuc2gnuc2(np, q.beg, tpt)
        _gnuc_end = tnuc2gnuc2(np, q.end, tpt)
        r.gnuc_beg = min(_gnuc_beg, _gnuc_end)
        r.gnuc_end = max(_gnuc_beg, _gnuc_end)
        refdupseq = faidx.refgenome.fetch_sequence(tpt.chrm, r.gnuc_beg, r.gnuc_end)
        natdupseq = reverse_complement(refdupseq) if tpt.strand == '-' else refdupseq
    if q.dupseq and natdupseq != q.dupseq:
        raise IncompatibleTranscriptError('unmatched reference')

    r.gnuc_range = '%d-%ddup%s' % (r.gnuc_beg, r.gnuc_end, refdupseq)
    r.pos = '%d-%d (dup)' % (r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s, %s)' % (tpt.gene.name, tpt.strand, tpt.region(r.gnuc_beg, r.gnuc_end))
    r.info = 'RefDupSeq=%s;NatDupSeq=%s' % (refdupseq, natdupseq)

    # assume no change in splice site
    if q.end.tpos == 0:
        # insert to exonic region, change to coding
        if (q.end.pos % 3 == 0 and len(natdupseq) % 3 == 0 and
            q.beg.tpos == 0 and q.end.tpos == 0):
            taa_dupseq = translate_seq(natdupseq)
            beg_index = (q.beg.pos + 2)/3
            end_index = (q.end.pos + 2)/3
            if beg_index == end_index:
                r.taa_range = '%s%ddup%s' % (taa_dupseq[0], beg_index, taa_dupseq)
            else:
                r.taa_range = '%s%d_%s%ddup%s' % (taa_dupseq[0], beg_index,
                                                  taa_dupseq[-1], end_index, taa_dupseq)
        else:
            q2 = copy(q)
            q2.insseq = natdupseq
            q2.pos = q.end
            r2 = Record()
            nuc_mutation_ins_coding(args, q2, tpt, r2)
            r.taa_range = r2.taa_range

    return r

def _core_annotate_nuc_dup(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = nuc_mutation_dup(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
        except UnknownChromosomeError:
            continue
        found = True
        r.tnuc_range = '%s_%sdup%s' % (q.beg, q.end, q.dupseq)
        r.format(q.op)

    if not found:
        r = Record()
        r.tnuc_range = '%s_%sdup%s' % (q.beg, q.end, q.dupseq)
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)


