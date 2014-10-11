""" annotate region """
from err import *
from record import *
from transcripts import *


def nuc_revanno_reg(args, q, tpt):

    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError('Transcript name unmatched')
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
    if q.ref and r.natrefseq != q.ref:
        raise IncompatibleTranscriptError()

    r.gnuc_range = '%d_%d%s' % (r.gnuc_beg, r.gnuc_end, r.refrefseq)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.tnuc_range = '%s_%s%s' % (q.beg, q.end, r.natrefseq)

    reg = ''
    if region_in_exon(np, q.beg, q.end):
        beg_codon = tpt.cpos2codon(q.beg.pos/3)
        end_codon = tpt.cpos2codon(q.end.pos/3)
        print q.beg, q.end, len(tpt)
        if beg_codon.index == end_codon.index:
            r.taa_pos = beg_codon.index
            r.taa_ref = codon2aa(beg_codon.seq)
        else:
            r.taa_ref = translate_seq(tpt.seq[beg_codon.index*3-3:end_codon.index*3])
            r.taa_range = 'p.%d_%d%s' % (beg_codon.index, end_codon.index, r.taa_ref)
    elif not region_in_intron(np, q.beg, q.end): # if block mutation occurs across splice site
        r.info = 'CrossSplitSite'
        r.reg = '%s (%s, coding;intronic)' % (tpt.gene.name, tpt.strand)
    else:
        r.reg = '%s (%s, intronic)' % (tpt.gene.name, tpt.strand)

    return r

def _core_annotate_nuc_reg(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = nuc_revanno_reg(args, q, tpt)
        except IncompatibleTranscriptError:
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



def codon_revanno_reg(args, q, tpt):

    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError('Transcript name unmatched')
    tpt.ensure_seq()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name

    if q.beg*3 > len(tpt) or q.end*3 > len(tpt):
        raise IncompatibleTranscriptError('codon nonexistent')

    tnuc_beg = q.beg*3 - 2
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
    r.tnuc_range = '%d_%d' % (tnuc_beg, tnuc_end)
    r.gnuc_range = '%d_%d' % (r.gnuc_beg, r.gnuc_end)
    r.taa_range = '%d_%d' % (q.beg, q.end)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
    r.info = 'PRefSeq=%s;NRefSeq=%s;RefSeq=%s' % (printseq(taa_natrefseq),
                                                  printseq(r.natrefseq),
                                                  printseq(r.refrefseq))

    return r


def _core_annotate_codon_reg(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = codon_revanno_reg(args, q, tpt)
        except IncompatibleTranscriptError:
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


