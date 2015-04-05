""" annotate block substitution,
assume this does not affect splice site """
from err import *
from record import *
from transcripts import *
from describe import *

def nuc_mutation_mnv(args, q, t, db):

    if q.tpt and t.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    t.ensure_seq()

    r = Record()
    r.chrm = t.chrm
    r.tname = t.format()
    r.gene = t.gene.name
    r.strand = t.strand

    t.ensure_position_array()
    check_exon_boundary(t.np, q.beg)
    check_exon_boundary(t.np, q.end)

    _gnuc_beg = t.tnuc2gnuc(q.beg)
    _gnuc_end = t.tnuc2gnuc(q.end)
    gnuc_beg = min(_gnuc_beg, _gnuc_end)
    gnuc_end = max(_gnuc_beg, _gnuc_end)
    r.pos = '%d-%d' % (gnuc_beg, gnuc_end)
    
    gnuc_refseq = faidx.getseq(t.chrm, gnuc_beg, gnuc_end)
    tnuc_refseq = reverse_complement(gnuc_refseq) if t.strand == '-' else gnuc_refseq
    gnuc_altseq = reverse_complement(q.altseq) if t.strand == '-' else q.altseq
    if q.refseq and tnuc_refseq != q.refseq:
        raise IncompatibleTranscriptError()

    r.gnuc_range = '%d_%d%s>%s' % (gnuc_beg, gnuc_end, gnuc_refseq, gnuc_altseq)
    r.tnuc_range = '%s_%s%s>%s' % (q.beg, q.end, tnuc_refseq, q.altseq)

    r.reg = describe_genic(args, t.chrm, gnuc_beg, gnuc_end, t, db)
    expt = r.set_splice()
    if (not expt) and r.reg.entirely_in_cds():
        try:
            t.tnuc_mnv_coding(q.beg.pos, q.end.pos, q.altseq, r)
        except IncompatibleTranscriptError as inst:
            _beg, _end, _seqlen = inst
            r.append_info('mnv_(%s-%s)_at_truncated_refseq_of_length_%d' % (_beg, _end, _seqlen))
        
    return r

def _core_annotate_nuc_mnv(args, q, tpts, db):

    found = False
    for t in tpts:
        try:
            r = nuc_mutation_mnv(args, q, t, db)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue
        except UnknownChromosomeError:
            continue
        found = True
        r.format(q.op)

    if not found:
        r = Record()
        r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))
        r.format(q.op)

    return

def codon_mutation_mnv(args, q, t):
    
    if q.tpt and t.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    t.ensure_seq()

    r = Record()
    r.chrm = t.chrm
    r.tname = t.format()
    r.gene = t.gene.name
    r.strand = t.strand

    if q.beg*3 > len(t) or q.end*3 > len(t):
        raise IncompatibleTranscriptError('codon nonexistent')

    tnuc_beg = q.beg*3-2
    tnuc_end = q.end*3
    gnuc_beg, gnuc_end = t.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
    tnuc_refseq = t.seq[tnuc_beg-1:tnuc_end]
    gnuc_refseq = reverse_complement(tnuc_refseq) if t.strand == '-' else tnuc_refseq
    taa_refseq = translate_seq(tnuc_refseq)
    if q.beg_aa and q.beg_aa != taa_refseq[0]:
        raise IncompatibleTranscriptError('beginning reference amino acid unmatched')
    if q.end_aa and q.end_aa != taa_refseq[-1]:
        raise IncompatibleTranscriptError('ending reference amino acid unmatched')
    if q.refseq and taa_refseq != q.refseq:
        raise IncompatibleTranscriptError('reference sequence unmatched')
    # reverse translate
    tnuc_altseq = []
    cdd_altseq = []
    for aa in q.altseq:
        tnuc_altseq.append(aa2codon(aa)[0])
        cdd_altseq.append('/'.join(aa2codon(aa)))
    tnuc_altseq = ''.join(tnuc_altseq)
    gnuc_altseq = reverse_complement(tnuc_altseq) if t.strand == '-' else tnuc_altseq
    r.tnuc_range = '%d_%d%s>%s' % (tnuc_beg, tnuc_end, tnuc_refseq, tnuc_altseq)
    r.gnuc_range = '%d_%d%s>%s' % (gnuc_beg, gnuc_end, gnuc_refseq, gnuc_altseq)
    r.pos = '%d-%d' % (gnuc_beg, gnuc_end)
    if len(cdd_altseq) <= 2:
        r.append_info('candidate_alternative_sequence=%s' % ('+'.join(cdd_altseq), ))

    return r

def _core_annotate_codon_mnv(args, q, tpts, db):

    found = False
    for t in tpts:
        try:
            r = codon_mutation_mnv(args, q, t)
        except IncompatibleTranscriptError:
            continue
        except UnknownChromosomeError:
            continue
        r.taa_range = '%s%s_%s%sdel%sins%s' % (
            q.beg_aa, str(q.beg), q.end_aa, str(q.end), q.refseq, q.altseq)
        r.reg = RegCDSAnno(t)
        r.reg.from_taa_range(q.beg, q.end)
        r.append_info('imprecise')
        r.format(q.op)
        found = True

    if not found:
        r = Record()
        r.taa_range = '%s%s_%s%sdel%sins%s' % (
            q.beg_aa, str(q.beg), q.end_aa, str(q.end), q.refseq, q.altseq)
        r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))

        r.format(q.op)

