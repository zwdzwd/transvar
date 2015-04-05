from transcripts import *
from utils import *
from record import *
from revanno_ins import nuc_mutation_ins_coding
from copy import copy

def nuc_mutation_dup(args, q, t):
    
    if q.tpt and t.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    t.ensure_seq()

    r = Record()
    r.chrm = t.chrm
    r.tname = t.format()
    r.gene = t.gene.name
    r.strand = t.strand

    if q.beg.pos > len(t) or q.end.pos > len(t):
        raise IncompatibleTranscriptError('codon nonexistent')

    t.ensure_position_array()
    check_exon_boundary(self.np, q.beg)
    check_exon_boundary(self.np, q.end)

    _gnuc_beg = t.tnuc2gnuc(q.beg)
    _gnuc_end = t.tnuc2gnuc(q.end)
    gnuc_beg = min(_gnuc_beg, _gnuc_end)
    gnuc_end = max(_gnuc_beg, _gnuc_end)
    gnuc_dupseq = faidx.getseq(t.chrm, gnuc_beg, gnuc_end)
    tnuc_dupseq = gnuc_dupseq if t.strand == '+' else reverse_complement(gnuc_dupseq)
    if q.dupseq and tnuc_dupseq != q.dupseq:
        raise IncompatibleTranscriptError('unmatched reference')

    if t.strand == '+':
        gnuc_ins = gnuc_set_ins(t.chrm, gnuc_end, gnuc_dupseq, r)
    else:
        gnuc_ins = gnuc_set_ins(t.chrm, gnuc_beg-1, gnuc_dupseq, r)

    tnuc_set_ins(gnuc_ins, t, r)
    r.reg = describe_genic(args, t.chrm, gnuc_ins.beg_r, gnuc_ins.end_r, t, db)
    tnuc_coding_ins(args, q, t, r, db)
    
    # # identify duplicate seq
    # if q.beg.tpos == 0 and q.end.tpos == 0:
    #     tnuc_dupseq = t.seq[q.beg.pos-1:q.end.pos]
    #     gnuc_dupseq = reverse_complement(tnuc_dupseq) if t.strand == '-' else tnuc_dupseq
    #     r.gnuc_beg, r.gnuc_end = tnuc_range2gnuc_range_(np, q.beg.pos, q.end.pos)
    # else:
    #     _gnuc_beg = tnuc2gnuc2(np, q.beg, t)
    #     _gnuc_end = tnuc2gnuc2(np, q.end, t)
    #     r.gnuc_beg = min(_gnuc_beg, _gnuc_end)
    #     r.gnuc_end = max(_gnuc_beg, _gnuc_end)
    #     gnuc_dupseq = faidx.getseq(t.chrm, r.gnuc_beg, r.gnuc_end)
    #     tnuc_dupseq = reverse_complement(gnuc_dupseq) if t.strand == '-' else gnuc_dupseq

    r.gnuc_range = '%d-%ddup%s' % (r.gnuc_beg, r.gnuc_end, gnuc_dupseq)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s, %s)' % (t.gene.name, t.strand, t.region(r.gnuc_beg, r.gnuc_end))
    r.append_info('duplication_gDNA=%s;duplication_cDNA=%s' % (gnuc_dupseq, tnuc_dupseq))

    # assume no change in splice site
    if q.end.tpos == 0:
        # insert to exonic region, change to coding
        if (q.end.pos % 3 == 0 and len(tnuc_dupseq) % 3 == 0 and
            q.beg.tpos == 0 and q.end.tpos == 0):
            taa_dupseq = translate_seq(tnuc_dupseq)
            beg_index = (q.beg.pos + 2)/3
            end_index = (q.end.pos + 2)/3
            if beg_index == end_index:
                r.taa_range = '%s%ddup%s' % (taa_dupseq[0], beg_index, taa_dupseq)
            else:
                r.taa_range = '%s%d_%s%ddup%s' % (taa_dupseq[0], beg_index,
                                                  taa_dupseq[-1], end_index, taa_dupseq)
        else:
            q2 = copy(q)
            q2.insseq = tnuc_dupseq
            q2.pos = q.end
            r2 = Record()
            nuc_mutation_ins_coding(args, q2, t, r2)
            r.taa_range = r2.taa_range

    return r

def _core_annotate_nuc_dup(args, q, tpts):

    found = False
    for t in tpts:
        try:
            r = nuc_mutation_dup(args, q, t)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue
        except UnknownChromosomeError:
            continue
        found = True
        r.tnuc_range = '%s_%sdup%s' % (q.beg, q.end, q.dupseq)
        r.format(q.op)

    if not found:
        r = Record()
        r.tnuc_range = '%s_%sdup%s' % (q.beg, q.end, q.dupseq)
        r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))
        r.format(q.op)


