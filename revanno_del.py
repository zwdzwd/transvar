from transcripts import *
from utils import *
from record import *
from copy import copy

def nuc_mutation_del_coding_inframe_inphase(args, q, tpt, r):

    """  deletion starts at the 1st base of the codon """

    beg_codon_index = (q.beg.pos + 2) / 3
    end_codon_index = (q.end.pos + 2) / 3
    if beg_codon_index == end_codon_index:
        beg_codon = tpt.cpos2codon(beg_codon_index)
        if not beg_codon: raise IncompatibleTranscriptError()
        end_codon = beg_codon
        r.taa_range = '%s%ddel' % (codon2aa(beg_codon.seq), beg_codon.index)
    else:
        beg_codon = tpt.cpos2codon(beg_codon_index)
        end_codon = tpt.cpos2codon(end_codon_index)
        if not beg_codon or not end_codon: raise IncompatibleTranscriptError()
        r.taa_range = '%s%d_%s%ddel' % (codon2aa(beg_codon.seq), beg_codon.index,
                                        codon2aa(end_codon.seq), end_codon.index)
    r.natdelseq = tpt.seq[q.beg.pos-1:q.end.pos]
    r.refdelseq = r.natdelseq if tpt.strand == '+' else reverse_complement(r.natdelseq)
    r.info = 'RefDelSeq=%s;NatDelSeq=%s' % (r.refdelseq, r.natdelseq)
    r.tnuc_range = '%d_%ddel' % (q.beg.pos, q.end.pos)
    if tpt.strand == '+':
        r.gnuc_beg = beg_codon.locs[0]
        r.gnuc_end = end_codon.locs[2]
    else:
        r.gnuc_beg = end_codon.locs[0]
        r.gnuc_end = beg_codon.locs[2]
    r.gnuc_range = '%d_%ddel' % (r.gnuc_beg, r.gnuc_end)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)

def nuc_mutation_del_coding_inframe_outphase(args, q, tpt, r):

    """ deletion starts at the 2nd/3rd base of the codon """

    beg_codon_index = (q.beg.pos + 2) / 3
    end_codon_index = (q.end.pos + 2) / 3
    beg_codon_beg = beg_codon_index*3 - 2
    end_codon_end = end_codon_index*3 # 1 past the last codon
    # print q.beg.pos, q.end.pos
    # print beg_codon_index, end_codon_index
    # print beg_codon_beg, end_codon_end
    newcodonseq = tpt.seq[beg_codon_beg-1:q.beg.pos-1]+tpt.seq[q.end.pos:end_codon_end]
    if len(newcodonseq) != 3: raise IncompatibleTranscriptError()
    r.taa_alt = codon2aa(newcodonseq)
    # beg_codon_seq = tpt.seq[beg_codon_beg-1:beg_codon_beg+2]
    # end_codon_seq = tpt.seq[end_codon_end-3:end_codon_end]
    tnuc_delseq = tpt.seq[beg_codon_beg-1:end_codon_end]
    taa_delseq = ''
    for i in xrange(len(tnuc_delseq)/3):
        aa = codon2aa(tnuc_delseq[i*3:i*3+3])
        taa_delseq += aa
    # print taa_delseq, end_codon_seq, beg_codon_seq, r.taa_alt
    if r.taa_alt == taa_delseq[-1]:
        # G100_S200delinsS becomes a pure deletion G100_D199del
        if len(taa_delseq) == 2:
            r.taa_range = '%s%ddel' % (taa_delseq[0], beg_codon_index)
        else:
            r.taa_range = '%s%d_%s%ddel' % (taa_delseq[0], beg_codon_index,
                                         taa_delseq[-2], end_codon_index-1)
    elif r.taa_alt == taa_delseq[0]:
        # S100_G200delinsS becomes a pure deletion D101_G200del
        if len(taa_delseq) == 2:
            r.taa_range = '%s%ddel' % (taa_delseq[1], beg_codon_index+1)
        else:
            r.taa_range = '%s%d_%s%ddel' % (taa_delseq[1], beg_codon_index+1,
                                         taa_delseq[-1], end_codon_index)
    else:
        r.taa_range = '%s%d_%s%ddelins%s' % (taa_delseq[0], beg_codon_index, 
                                             taa_delseq[-1], end_codon_index, r.taa_alt)

    r.natdelseq = tpt.seq[q.beg.pos-1:q.end.pos]
    r.refdelseq = r.natdelseq if tpt.strand == '+' else reverse_complement(r.natdelseq)
    beg_codon = tpt.cpos2codon(beg_codon_index)
    end_codon = tpt.cpos2codon(end_codon_index)
    r.tnuc_range = '%d_%ddel' % (q.beg.pos, q.end.pos)
    r.info = 'RefDelSeq=%s;NatDelSeq=%s' % (r.refdelseq, r.natdelseq)
    gnuc_del_beg = reverse_tnuc_pos(beg_codon, q.beg.pos)
    gnuc_del_end = reverse_tnuc_pos(end_codon, q.end.pos)
    if tpt.strand == '+':
        r.gnuc_beg = gnuc_del_beg
        r.gnuc_end = gnuc_del_end
    else:
        r.gnuc_beg = gnuc_del_end
        r.gnuc_end = gnuc_del_beg
    r.gnuc_range = '%d_%ddel' % (r.gnuc_beg, r.gnuc_end)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
    # print beg_codon, beg_codon.seq, end_codon, end_codon.seq

def nuc_mutation_del_coding_inframe(args, q, tpt, r):

    if q.beg.pos % 3 == 1:
        nuc_mutation_del_coding_inframe_inphase(args, q, tpt, r)
    else:
        nuc_mutation_del_coding_inframe_outphase(args, q, tpt, r)

def nuc_mutation_del_coding_frameshift(args, q, tpt, r):

    # assume frame-shift does not affect splicing
    r.natdelseq = tpt.seq[q.beg.pos-1:q.end.pos]
    if q.delseq and r.natdelseq != q.delseq: raise IncompatibleTranscriptError()
    beg_codon_index = (q.beg.pos + 2) / 3
    beg_codon_beg = beg_codon_index*3 - 2
    old_seq = tpt.seq[beg_codon_beg-1:]
    new_seq = tpt.seq[beg_codon_beg-1:q.beg.pos-1]+tpt.seq[q.end.pos:]
    if not old_seq: raise IncompatibleTranscriptError()
    ret = tpt.extend_taa_seq(beg_codon_index, old_seq, new_seq)
    if ret:
        taa_pos, taa_ref, taa_alt, termlen = ret
        r.taa_range = '%s%d%sfs*%s' % (taa_ref, taa_pos, taa_alt, termlen)
    else:
        r.taa_range = '(=)'
    r.tnuc_range = '%d_%ddel' % (q.beg.pos, q.end.pos)

    r.gnuc_beg, r.gnuc_end = tpt.tnuc_range2gnuc_range(q.beg.pos, q.end.pos)
    if r.gnuc_beg == r.gnuc_end:
        r.gnuc_range = '%ddel' % r.gnuc_beg
    else:
        r.gnuc_range = '%d_%ddel' % (r.gnuc_beg, r.gnuc_end)

    r.refdelseq = r.natdelseq if tpt.strand == '+' else reverse_complement(r.natdelseq)
    r.info = 'RefDelSeq=%s;NatDelSeq=%s' % (r.refdelseq, r.natdelseq)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)

def nuc_mutation_del_coding(args, q, tpt, r):

    """ assume no intron between q.beg and q.end """
    if (q.end.pos - q.beg.pos) % 3 == 2: # in-frame
        nuc_mutation_del_coding_inframe(args, q, tpt, r)
    else:   # frame-shift
        nuc_mutation_del_coding_frameshift(args, q, tpt, r)


def nuc_mutation_del(args, q, tpt):

    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    tpt.ensure_seq()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name
    r.muttype = 'del'

    # if q.beg.pos == q.end.pos+43241: # with respect to one exome boundary
    #     _nuc_mutation_del_intronic(args, q, tpt, r)
    # else:
    np = tpt.position_array()
    check_exon_boundary(np, q.beg)
    check_exon_boundary(np, q.end)

    gnuc_beg = tnuc2gnuc2(np, q.beg, tpt)
    gnuc_end = tnuc2gnuc2(np, q.end, tpt)
    r.gnuc_beg = min(gnuc_beg, gnuc_end)
    r.gnuc_end = max(gnuc_beg, gnuc_end)
    r.gnuc_range = '%d_%ddel' % (r.gnuc_beg, r.gnuc_end)
    tnuc_coding_beg = q.beg.pos if q.beg.tpos <= 0 else q.beg.pos+1
    tnuc_coding_end = q.end.pos if q.end.tpos >= 0 else q.end.pos-1
    gnuc_coding_beg = tnuc2gnuc(np, tnuc_coding_beg)
    gnuc_coding_end = tnuc2gnuc(np, tnuc_coding_end)

    refdelseq = faidx.refgenome.fetch_sequence(tpt.chrm, r.gnuc_beg, r.gnuc_end)
    natdelseq = refdelseq if tpt.strand == '+' else reverse_complement(refdelseq)
    if q.delseq and natdelseq != q.delseq:
        raise IncompatibleTranscriptError()

    reg = ''

    # if deletion affects coding region
    if tnuc_coding_beg <= tnuc_coding_end:
        q_coding = copy(q)
        q_coding.beg = Pos(tnuc_coding_beg)
        q_coding.end = Pos(tnuc_coding_end)

        # set deletion sequence in the coding region
        # [gnuc_beg] --- [gnuc_coding_beg] ===---=== [gnuc_coding_end] --- [gnuc_end]
        delseq_beg = abs(gnuc_coding_beg - gnuc_beg)
        delseq_end = len(q.delseq) + abs(gnuc_coding_end - gnuc_end)
        q_coding.delseq = q.delseq[delseq_beg:delseq_end]

        r_coding = Record()
        nuc_mutation_del_coding(args, q_coding, tpt, r_coding)
        r.taa_range = r_coding.taa_range
        reg += ',coding' if reg else 'coding'
    else:
        r.taa_range = ''

    if r.gnuc_beg == r.gnuc_end:
        r.tnuc_range = '%sdel' % (q.beg, )
    else:
        r.tnuc_range = '%s_%sdel' % (q.beg, q.end)

    if q.beg.tpos != 0 or q.end.tpos != 0:
        reg += ',intronic' if reg else 'intronic'
    r.reg = '%s (%s %s)' % (tpt.gene.name, tpt.strand, reg)
    r.info = 'RefDelSeq=%s;NatDelSeq=%s' % (refdelseq, natdelseq)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)

    return r

def _core_annotate_nuc_del(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = nuc_mutation_del(args, q, tpt)
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

def codon_mutation_del(args, q, tpt):
    
    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    tpt.ensure_seq()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name

    if q.beg*3 > len(tpt) or q.end*3 > len(tpt):
        raise IncompatibleTranscriptError('codon nonexistent')

    if q.delseq and tpt.taa_range2aa_seq(q.beg, q.end) != q.delseq:
        raise IncompatibleTranscriptError('unmatched reference')
    
    tnuc_beg = q.beg*3-2
    tnuc_end = q.end*3
    gnuc_beg, gnuc_end = tpt.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
    r.tnuc_range = '%d_%d' % (tnuc_beg, tnuc_end)
    r.gnuc_range = '%d_%d' % (gnuc_beg, gnuc_end)
    r.pos = '%d-%d (deletion)' % (gnuc_beg, gnuc_end)

    return r

def _core_annotate_codon_del(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = codon_mutation_del(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue
        except UnknownChromosomeError:
            continue
        r.muttype = 'del'
        r.taa_range = '%s%s_%s%sdel%s' % (q.beg_aa if q.beg_aa else '', str(q.beg),
                                          q.end_aa if q.end_aa else '', str(q.end),
                                          q.delseq)
        r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
        r.info = 'Uncertain' if r.info == '.' else (r.info+';Uncertain')
        r.format(q.op)
        found = True

    if not found:
        r = Record()
        r.taa_range = '%s%s_%s%sdel%s' % (q.beg_aa, str(q.beg), q.end_aa, str(q.end), q.delseq)
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)



