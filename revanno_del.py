from transcripts import *
from utils import *
from record import *
from copy import copy
from describe import *

def nuc_mutation_del_coding_inframe_inphase(args, q, t, r):

    """  deletion starts at the 1st base of the codon """

    beg_codon_index = (q.beg.pos + 2) / 3
    end_codon_index = (q.end.pos + 2) / 3
    if beg_codon_index == end_codon_index:
        beg_codon = t.cpos2codon(beg_codon_index)
        if not beg_codon: raise IncompatibleTranscriptError()
        end_codon = beg_codon
        r.taa_range = '%s%ddel' % (codon2aa(beg_codon.seq), beg_codon.index)
    else:
        beg_codon = t.cpos2codon(beg_codon_index)
        end_codon = t.cpos2codon(end_codon_index)
        if not beg_codon or not end_codon: raise IncompatibleTranscriptError()
        r.taa_range = '%s%d_%s%ddel' % (codon2aa(beg_codon.seq), beg_codon.index,
                                        codon2aa(end_codon.seq), end_codon.index)
    r.natdelseq = t.seq[q.beg.pos-1:q.end.pos]
    r.refdelseq = r.natdelseq if t.strand == '+' else reverse_complement(r.natdelseq)
    r.info = 'RefDelSeq=%s;NatDelSeq=%s' % (r.refdelseq, r.natdelseq)
    r.tnuc_range = '%d_%ddel' % (q.beg.pos, q.end.pos)
    if t.strand == '+':
        r.gnuc_beg = beg_codon.locs[0]
        r.gnuc_end = end_codon.locs[2]
    else:
        r.gnuc_beg = end_codon.locs[0]
        r.gnuc_end = beg_codon.locs[2]
    r.gnuc_range = '%d_%ddel' % (r.gnuc_beg, r.gnuc_end)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s, coding)' % (t.gene.name, t.strand)

def nuc_mutation_del_coding_inframe_outphase(args, q, t, r):

    """ deletion starts at the 2nd/3rd base of the codon """

    beg_codon_index = (q.beg.pos + 2) / 3
    end_codon_index = (q.end.pos + 2) / 3
    beg_codon_beg = beg_codon_index*3 - 2
    end_codon_end = end_codon_index*3 # 1 past the last codon
    # print q.beg.pos, q.end.pos
    # print beg_codon_index, end_codon_index
    # print beg_codon_beg, end_codon_end
    newcodonseq = t.seq[beg_codon_beg-1:q.beg.pos-1]+t.seq[q.end.pos:end_codon_end]
    if len(newcodonseq) != 3: raise IncompatibleTranscriptError()
    r.taa_alt = codon2aa(newcodonseq)
    # beg_codon_seq = t.seq[beg_codon_beg-1:beg_codon_beg+2]
    # end_codon_seq = t.seq[end_codon_end-3:end_codon_end]
    tnuc_delseq = t.seq[beg_codon_beg-1:end_codon_end]
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

    r.natdelseq = t.seq[q.beg.pos-1:q.end.pos]
    r.refdelseq = r.natdelseq if t.strand == '+' else reverse_complement(r.natdelseq)
    beg_codon = t.cpos2codon(beg_codon_index)
    end_codon = t.cpos2codon(end_codon_index)
    r.tnuc_range = '%d_%ddel' % (q.beg.pos, q.end.pos)
    r.info = 'RefDelSeq=%s;NatDelSeq=%s' % (r.refdelseq, r.natdelseq)
    gnuc_del_beg = reverse_tnuc_pos(beg_codon, q.beg.pos)
    gnuc_del_end = reverse_tnuc_pos(end_codon, q.end.pos)
    if t.strand == '+':
        r.gnuc_beg = gnuc_del_beg
        r.gnuc_end = gnuc_del_end
    else:
        r.gnuc_beg = gnuc_del_end
        r.gnuc_end = gnuc_del_beg
    r.gnuc_range = '%d_%ddel' % (r.gnuc_beg, r.gnuc_end)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s, coding)' % (t.gene.name, t.strand)
    # print beg_codon, beg_codon.seq, end_codon, end_codon.seq

def nuc_mutation_del_coding_inframe(args, q, t, r):

    if q.beg.pos % 3 == 1:
        nuc_mutation_del_coding_inframe_inphase(args, q, t, r)
    else:
        nuc_mutation_del_coding_inframe_outphase(args, q, t, r)

def nuc_mutation_del_coding_frameshift(args, q, t, r):

    # assume frame-shift does not affect splicing
    r.natdelseq = t.seq[q.beg.pos-1:q.end.pos]
    if q.delseq and r.natdelseq != q.delseq: raise IncompatibleTranscriptError()
    beg_codon_index = (q.beg.pos + 2) / 3
    beg_codon_beg = beg_codon_index*3 - 2
    old_seq = t.seq[beg_codon_beg-1:]
    new_seq = t.seq[beg_codon_beg-1:q.beg.pos-1]+t.seq[q.end.pos:]
    if not old_seq: raise IncompatibleTranscriptError()
    ret = t.extend_taa_seq(beg_codon_index, old_seq, new_seq)
    if ret:
        taa_pos, taa_ref, taa_alt, termlen = ret
        r.taa_range = '%s%d%sfs*%s' % (taa_ref, taa_pos, taa_alt, termlen)
    else:
        r.taa_range = '(=)'
    r.tnuc_range = '%d_%ddel' % (q.beg.pos, q.end.pos)

    r.gnuc_beg, r.gnuc_end = t.tnuc_range2gnuc_range(q.beg.pos, q.end.pos)
    if r.gnuc_beg == r.gnuc_end:
        r.gnuc_range = '%ddel' % r.gnuc_beg
    else:
        r.gnuc_range = '%d_%ddel' % (r.gnuc_beg, r.gnuc_end)

    r.refdelseq = r.natdelseq if t.strand == '+' else reverse_complement(r.natdelseq)
    r.info = 'RefDelSeq=%s;NatDelSeq=%s' % (r.refdelseq, r.natdelseq)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s, coding)' % (t.gene.name, t.strand)

def nuc_mutation_del_coding(args, q, t, r):

    """ assume no intron between q.beg and q.end """
    if (q.end.pos - q.beg.pos) % 3 == 2: # in-frame
        nuc_mutation_del_coding_inframe(args, q, t, r)
    else:   # frame-shift
        nuc_mutation_del_coding_frameshift(args, q, t, r)

def nuc_mutation_del(args, q, r, t, db):

    # check exon boundary
    t.ensure_position_array()
    check_exon_boundary(t.np, q.beg)
    check_exon_boundary(t.np, q.end)

    _gnuc_beg = t.tnuc2gnuc(q.beg)
    _gnuc_end = t.tnuc2gnuc(q.end)
    gnuc_beg = min(_gnuc_beg, _gnuc_end)
    gnuc_end = max(_gnuc_beg, _gnuc_end)

    gnuc_delseq = faidx.getseq(t.chrm, gnuc_beg, gnuc_end)
    tnuc_delseq = gnuc_delseq if t.strand == '+' else reverse_complement(gnuc_delseq)
    if q.delseq and tnuc_delseq != q.delseq:
        raise IncompatibleTranscriptError()

    # right-align
    gnuc_beg_r, gnuc_end_r = gnuc_roll_right_del(t.chrm, gnuc_beg, gnuc_end)
    gnuc_delseq_r = faidx.getseq(t.chrm, gnuc_beg_r, gnuc_end_r)
    r.gnuc_range = gnuc_del_id(t.chrm, gnuc_beg_r, gnuc_end_r)
    r.pos = '%d_%d' % (gnuc_beg_r, gnuc_end_r)

    # left-align
    gnuc_beg_l, gnuc_end_l = gnuc_roll_left_del(t.chrm, gnuc_beg, gnuc_end)
    gnuc_delseq_l = faidx.getseq(t.chrm, gnuc_beg_l, gnuc_end_l)
    r.append_info('left_align_gDNA=g.%s' % gnuc_del_id(t.chrm, gnuc_beg_l, gnuc_end_l))
    r.append_info('unaligned_gDNA=g.%s' % gnuc_del_id(t.chrm, gnuc_beg, gnuc_end))

    if t.strand == '+':
        c1l, p1l = t.gpos2codon(gnuc_beg_l)
        c2l, p2l = t.gpos2codon(gnuc_end_l)
        tnuc_delseq_l = gnuc_delseq_l
        c1r, p1r = t.gpos2codon(gnuc_beg_r)
        c2r, p2r = t.gpos2codon(gnuc_end_r)
        tnuc_delseq_r = gnuc_delseq_r
    else:
        c1l, p1l = t.gpos2codon(gnuc_end_r)
        c2l, p2l = t.gpos2codon(gnuc_beg_r)
        tnuc_delseq_l = reverse_complement(gnuc_delseq_r)
        c1r, p1r = t.gpos2codon(gnuc_end_l)
        c2r, p2r = t.gpos2codon(gnuc_beg_l)
        tnuc_delseq_r = reverse_complement(gnuc_delseq_l)

    # cDNA representation
    # right align
    r.tnuc_range = tnuc_del_id(p1r, p2r, tnuc_delseq_r)
    # left align
    r.append_info('left_align_cDNA=c.%s' % tnuc_del_id(p1l, p2l, tnuc_delseq_l))
    r.append_info('unalign_cDNA=c.%s' % tnuc_del_id(q.beg, q.end, tnuc_delseq))

    # tnuc_coding_beg = q.beg.pos if q.beg.tpos <= 0 else q.beg.pos+1
    # tnuc_coding_end = q.end.pos if q.end.tpos >= 0 else q.end.pos-1
    # gnuc_coding_beg = tnuc2gnuc(np, tnuc_coding_beg)
    # gnuc_coding_end = tnuc2gnuc(np, tnuc_coding_end)

    r.reg = describe_genic(args, t.chrm, gnuc_beg, gnuc_end, t, db)

    # if deletion affects coding region
    if t.transcript_type == 'protein_coding' and not same_intron(q.beg, q.end):
        expt = r.set_splice('lost')
        if not expt:
            c1, p1 = t.intronic_lean(q.beg, 'c_greater')
            c2, p2 = t.intronic_lean(q.end, 'c_smaller')

            if (gnuc_end - gnuc_beg + 1) % 3 == 0:
                del_coding_inframe(args, c1, c2, p1, p2, t, r)
            else:
                del_coding_frameshift(args, c1, c2, p1, p2, t, r)

    r.append_info('deletion_gDNA=%s;deletion_cDNA=%s' % (gnuc_delseq_r, tnuc_delseq_r))
    # if tnuc_coding_beg <= tnuc_coding_end:
    #     q_coding = copy(q)
    #     q_coding.beg = Pos(tnuc_coding_beg)
    #     q_coding.end = Pos(tnuc_coding_end)

    #     # set deletion sequence in the coding region
    #     # [gnuc_beg] --- [gnuc_coding_beg] ===---=== [gnuc_coding_end] --- [gnuc_end]
    #     delseq_beg = abs(gnuc_coding_beg - gnuc_beg)
    #     delseq_end = len(q.delseq) + abs(gnuc_coding_end - gnuc_end)
    #     q_coding.delseq = q.delseq[delseq_beg:delseq_end]

    #     r_coding = Record()
    #     nuc_mutation_del_coding(args, q_coding, t, r_coding)

    #     r.taa_range = r_coding.taa_range

    # if q.beg.tpos != 0 or q.end.tpos != 0:
    #     reg += ',intronic' if reg else 'intronic'
    # r.reg = '%s (%s %s)' % (t.gene.name, t.strand, reg)

def _core_annotate_nuc_del(args, q, tpts, db):

    found = False
    for t in tpts:
        if q.tpt and t.name != q.tpt:
            raise IncompatibleTranscriptError("Transcript name unmatched")

        r = Record()
        r.chrm = t.chrm
        r.tname = t.name
        r.gene = t.gene.name
        r.strand = t.strand

        try:
            t.ensure_seq()
            nuc_mutation_del(args, q, r, t, db)
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
        tnuc_del_id(q.beg, q.end, q.delseq)
        r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))
        r.format(q.op)

    return

def codon_mutation_del(args, q, t, db):
    
    if q.tpt and t.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    t.ensure_seq()

    r = Record()
    r.chrm = t.chrm
    r.tname = t.name

    if q.beg*3 > len(t) or q.end*3 > len(t):
        raise IncompatibleTranscriptError('codon nonexistent')

    if q.delseq and t.taa_range2aa_seq(q.beg, q.end) != q.delseq:
        raise IncompatibleTranscriptError('unmatched reference')
    
    tnuc_beg = q.beg*3-2
    tnuc_end = q.end*3
    gnuc_beg, gnuc_end = t.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
    r.tnuc_range = '%d_%d' % (tnuc_beg, tnuc_end)
    r.gnuc_range = '%d_%d' % (gnuc_beg, gnuc_end)
    r.pos = '%d-%d (deletion)' % (gnuc_beg, gnuc_end)

    return r

def _core_annotate_codon_del(args, q, tpts, db):

    found = False
    for t in tpts:
        try:
            r = codon_mutation_del(args, q, t, db)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue
        except UnknownChromosomeError:
            continue

        r.taa_range = '%s%s_%s%sdel%s' % (q.beg_aa if q.beg_aa else '', str(q.beg),
                                          q.end_aa if q.end_aa else '', str(q.end),
                                          q.delseq)
        r.reg = '%s (%s, coding)' % (t.gene.name, t.strand)
        r.info = 'Uncertain' if r.info == '.' else (r.info+';Uncertain')
        r.format(q.op)
        found = True

    if not found:
        r = Record()
        r.taa_range = '%s%s_%s%sdel%s' % (q.beg_aa, str(q.beg), q.end_aa, str(q.end), q.delseq)
        r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))

        r.format(q.op)



