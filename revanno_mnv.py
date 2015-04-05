""" annotate block substitution,
assume this does not affect splice site """
from err import *
from record import *
from transcripts import *

# def nuc_mutation_mnv_coding_inframe(beg, end, altseq, tpt, r):


# def nuc_mutation_mnv_coding_frameshift(beg, end, altseq, tpt, r):

#     beg_codon_index = (beg + 2) / 3
#     beg_codon_beg = beg_codon_index * 3 - 2
#     old_seq = tpt.seq[beg_codon_beg-1:]
#     new_seq = tpt.seq[beg_codon_beg-1:beg-1]+altseq+tpt.seq[end:]

#     ret = tpt.extend_taa_seq(beg_codon_beg, old_seq, new_seq)
#     if ret:
#         taa_pos, taa_ref, taa_alt, termlen = ret
#         r.taa_range = '%s%d%sfs*%s' % (taa_ref, taa_pos, taa_alt, termlen)
#     else:
#         r.taa_range = '(=)'

#     return

# def nuc_mutation_mnv_coding(beg, end, altseq, tpt, r):

#     if (len(q.altseq) - (q.end.pos-q.beg.pos+1))  % 3 == 0:
#         nuc_mutation_mnv_coding_inframe(beg, end, altseq, tpt, r)
#     else:
#         nuc_mutation_mnv_coding_frameshift(beg, end, altseq, tpt, r)

def nuc_mutation_mnv(args, q, tpt):

    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    tpt.ensure_seq()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name
    r.muttype = 'mnv'

    np = tpt.position_array()
    check_exon_boundary(np, q.beg)
    check_exon_boundary(np, q.end)

    gnuc_beg = tnuc2gnuc2(np, q.beg, tpt)
    gnuc_end = tnuc2gnuc2(np, q.end, tpt)
    r.gnuc_beg = min(gnuc_beg, gnuc_end)
    r.gnuc_end = max(gnuc_beg, gnuc_end)
    tnuc_coding_beg = q.beg.pos if q.beg.tpos <= 0 else q.beg.pos+1 # base after intron
    tnuc_coding_end = q.end.pos if q.end.tpos >= 0 else q.end.pos-1 # base before intron
    gnuc_coding_beg = tnuc2gnuc(np, tnuc_coding_beg)
    gnuc_coding_end = tnuc2gnuc(np, tnuc_coding_end)

    r.refrefseq = faidx.refgenome.fetch_sequence(tpt.chrm, r.gnuc_beg, r.gnuc_end)
    r.natrefseq = reverse_complement(r.refrefseq) if tpt.strand == '-' else r.refrefseq
    r.refaltseq = reverse_complement(q.altseq) if tpt.strand == '-' else q.altseq
    if q.refseq and r.natrefseq != q.refseq:
        raise IncompatibleTranscriptError()

    r.gnuc_range = '%d_%d%s>%s' % (r.gnuc_beg, r.gnuc_end, r.refrefseq, r.refaltseq)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.tnuc_range = '%s_%s%s>%s' % (q.beg, q.end, r.natrefseq, q.altseq)

    r._reg_ = tpt.describe_span(r.gnuc_beg, r.gnuc_end, args)
    r.reg = '%s (%s, %s)' % (tpt.gene.name, tpt.strand, r._reg_.format())
    if r._reg_.in_exon():
        r.taa_range = tpt.tnuc_mnv_coding(q.beg.pos, q.end.pos, q.altseq, r)

    # if tnuc_region_in_exon(np, q.beg, q.end):
    #     .tnunuc_mutation_mnv_coding(q.beg.pos, q.end.pos, q.altseq, tpt, r)
    #     r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
    # elif not tnuc_region_in_intron(np, q.beg, q.end): # if block mutation occurs across splice site
    #     r.info = 'CrossSpliceSite'
    #     r.reg = '%s (%s, coding;intronic)' % (tpt.gene.name, tpt.strand)
    #     # raise UnImplementedError('Block mutation occurs across splice site')
    # else:
    #     r.reg = '%s (%s, intronic)' % (tpt.gene.name, tpt.strand)

    return r

def _core_annotate_nuc_mnv(args, q, tpts, db):

    found = False
    for tpt in tpts:
        try:
            r = nuc_mutation_mnv(args, q, tpt)
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
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

    return

def codon_mutation_mnv(args, q, tpt):
    
    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    tpt.ensure_seq()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name

    if q.beg*3 > len(tpt) or q.end*3 > len(tpt):
        raise IncompatibleTranscriptError('codon nonexistent')

    tnuc_beg = q.beg*3-2
    tnuc_end = q.end*3
    r.gnuc_beg, r.gnuc_end = tpt.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
    r.natrefseq = tpt.seq[tnuc_beg-1:tnuc_end]
    r.refrefseq = reverse_complement(r.natrefseq) if tpt.strand == '-' else r.natrefseq
    taa_natrefseq = translate_seq(r.natrefseq)
    if q.beg_aa and q.beg_aa != taa_natrefseq[0]:
        raise IncompatibleTranscriptError('beginning reference amino acid unmatched')
    if q.end_aa and q.end_aa != taa_natrefseq[-1]:
        raise IncompatibleTranscriptError('ending reference amino acid unmatched')
    if q.refseq and taa_natrefseq != q.refseq:
        raise IncompatibleTranscriptError('reference sequence unmatched')
    # reverse translate
    tnuc_altseq = []
    cdd_altseq = []
    for aa in q.altseq:
        tnuc_altseq.append(aa2codon(aa)[0])
        cdd_altseq.append('/'.join(aa2codon(aa)))
    r.nataltseq = ''.join(tnuc_altseq)
    r.refaltseq = reverse_complement(r.nataltseq) if tpt.strand == '-' else r.nataltseq
    r.tnuc_range = '%d_%d%s>%s' % (tnuc_beg, tnuc_end, r.natrefseq, r.nataltseq)
    r.gnuc_range = '%d_%d%s>%s' % (r.gnuc_beg, r.gnuc_end, r.refrefseq, r.refaltseq)
    r.pos = '%d-%d (block substitution)' % (r.gnuc_beg, r.gnuc_end)
    r.info = 'CddNatAlt=%s' % ('+'.join(cdd_altseq), )

    return r

def _core_annotate_codon_mnv(args, q, tpts, db):

    found = False
    for tpt in tpts:
        try:
            r = codon_mutation_mnv(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
        except UnknownChromosomeError:
            continue
        r.muttype = 'delins'
        r.taa_range = '%s%s_%s%sdel%sins%s' % (
            q.beg_aa, str(q.beg), q.end_aa, str(q.end), q.refseq, q.altseq)
        r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
        r.info = 'Uncertain' if r.info == '.' else (r.info+';Uncertain')
        r.format(q.op)
        found = True

    if not found:
        r = Record()
        r.taa_range = '%s%s_%s%sdel%sins%s' % (
            q.beg_aa, str(q.beg), q.end_aa, str(q.end), q.refseq, q.altseq)
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

