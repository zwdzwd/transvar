""" annotate region """
from err import *
from record import *
from transcripts import *
from describe import *

def nuc_revanno_reg(args, q, tpt, db):

    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError('Transcript name unmatched')
    tpt.ensure_seq()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name
    # r.muttype = 'mnv'
    r.gene = tpt.gene.name
    r.strand = tpt.strand

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

    r.gnuc_range = '%d_%d%s' % (r.gnuc_beg, r.gnuc_end, r.refrefseq) if r.gnuc_beg != r.gnuc_end else '%d%s' % (r.gnuc_beg, r.refrefseq)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.tnuc_range = '%s_%s%s' % (q.beg, q.end, r.natrefseq) if q.beg != q.end else '%s%s' % (q.beg, r.natrefseq)

    r.reg = describe_genic(args, tpt.chrm, r.gnuc_beg, r.gnuc_end, tpt, db)
    expt = r.set_splice('included')
    if hasattr(r.reg, 'cover_cds') and r.reg.cover_cds: #region_in_exon(np, q.beg, q.end):
        c1, p1 = tpt.intronic_lean(q.beg, 'c_greater')
        c2, p2 = tpt.intronic_lean(q.end, 'c_smaller')
        
        # beg_codon = tpt.cpos2codon(q.beg.pos/3)
        # end_codon = tpt.cpos2codon(q.end.pos/3)
        if c1.index == c2.index:
            r.taa_pos = c1.index
            r.taa_ref = codon2aa(c1.seq)
        else:
            r.taa_ref = translate_seq(tpt.seq[c1.index*3-3:c2.index*3])
            r.taa_range = '%d_%d%s' % (c1.index, c2.index, r.taa_ref)
    # elif not region_in_intron(np, q.beg, q.end): # if block mutation occurs across splice site
    #     r.append_info('CrossSplitSite')
    #     r.reg = '%s (%s, coding;intronic)' % (tpt.gene.name, tpt.strand)
    # else:
    #     r.reg = '%s (%s, intronic)' % (tpt.gene.name, tpt.strand)

    if tpt.gene.dbxref:
        r.append_info('DBXref=%s' % tpt.gene.dbxref)

    return r

def _core_annotate_nuc_reg(args, q, tpts, db):

    found = False
    for tpt in tpts:
        try:
            r = nuc_revanno_reg(args, q, tpt, db)
        except IncompatibleTranscriptError:
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
    if q.refseq and not re.match(q.refseq.replace('x','[A-Z]'), taa_natrefseq): # != q.refseq:
        raise IncompatibleTranscriptError('reference sequence unmatched')
    r.tnuc_range = '%d_%d' % (tnuc_beg, tnuc_end)
    r.gnuc_range = '%d_%d' % (r.gnuc_beg, r.gnuc_end)
    r.taa_range = '%s%d_%s%d' % (taa_natrefseq[0], q.beg, taa_natrefseq[-1], q.end) if q.beg != q.end else '%d%s' % (q.beg, taa_natrefseq[0])
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.gene = tpt.gene.name
    r.strand = tpt.strand

    r.reg = RegCDSAnno(tpt)
    r.reg.from_taa_range(q.beg, q.end)
    r.append_info('protein_sequence=%s;cDNA_sequence=%s;gDNA_sequence=%s' % (printseq(taa_natrefseq), printseq(r.natrefseq), printseq(r.refrefseq)))
    if tpt.gene.dbxref:
        r.append_info('DBXref=%s' % tpt.gene.dbxref)

    return r


def _core_annotate_codon_reg(args, q, tpts, db):

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
        r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))

        r.format(q.op)

    return


