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
        r.taa_range = '%s%ddel' % (standard_codon_table[beg_codon.seq], 
                                   beg_codon.index)
    else:
        beg_codon = tpt.cpos2codon(beg_codon_index)
        end_codon = tpt.cpos2codon(end_codon_index)
        if not beg_codon or not end_codon: raise IncompatibleTranscriptError()
        r.taa_range = '%s%d_%s%ddel' % (standard_codon_table[beg_codon.seq],
                                        beg_codon.index,
                                        standard_codon_table[end_codon.seq],
                                        end_codon.index)
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
    r.pos = '%s:%d-%d' % (r.chrm, r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)

def nuc_mutation_del_coding_inframe_outphase(args, q, tpt, r):

    """ deletion starts at the 2nd/3rd base of the codon """

    r.muttype = 'del'
    beg_codon_index = (q.beg.pos + 2) / 3
    end_codon_index = (q.end.pos + 2) / 3
    beg_codon_beg = beg_codon_index*3 - 2
    end_codon_end = end_codon_index*3 # 1 past the last codon
    # print q.beg.pos, q.end.pos
    # print beg_codon_index, end_codon_index
    # print beg_codon_beg, end_codon_end
    newcodonseq = tpt.seq[beg_codon_beg-1:q.beg.pos-1]+tpt.seq[q.end.pos:end_codon_end]
    if len(newcodonseq) != 3: raise IncompatibleTranscriptError()
    r.taa_alt = standard_codon_table[newcodonseq]
    beg_codon_seq = tpt.seq[beg_codon_beg-1:beg_codon_beg+2]
    end_codon_seq = tpt.seq[end_codon_end-3:end_codon_end]
    tnuc_delseq = tpt.seq[beg_codon_beg-1:end_codon_end]
    taa_delseq = ''
    for i in xrange(len(tnuc_delseq)/3):
        aa = standard_codon_table[tnuc_delseq[i*3:i*3+3]]
        taa_delseq += aa
    print taa_delseq, end_codon_seq, beg_codon_seq, r.taa_alt
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
    r.pos = '%s:%d-%d' % (r.chrm, r.gnuc_beg, r.gnuc_end)
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
    ret = extend_taa_seq(beg_codon_index, old_seq, new_seq, tpt)
    if ret:
        taa_pos, taa_ref, taa_alt, termlen = ret
        r.taa_range = '%s%d%sfs*%s' % (taa_ref, taa_pos, taa_alt, termlen)
    else:
        r.taa_range = 'synonymous'
    r.tnuc_range = '%d_%ddel' % (q.beg.pos, q.end.pos)

    r.gnuc_beg, r.gnuc_end = tpt.tnuc_range2gnuc_range(q.beg.pos, q.end.pos)
    if r.gnuc_beg == r.gnuc_end:
        r.gnuc_range = '%ddel' % r.gnuc_beg
    else:
        r.gnuc_range = '%d_%ddel' % (r.gnuc_beg, r.gnuc_end)

    r.refdelseq = r.natdelseq if tpt.strand == '+' else reverse_complement(r.natdelseq)
    r.info = 'RefDelSeq=%s;NatDelSeq=%s' % (r.refdelseq, r.natdelseq)
    r.pos = '%s:%d-%d' % (r.chrm, r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)

def nuc_mutation_del_coding(args, q, tpt, r):

    """ assume no intron between q.beg and q.end """
    if (q.end.pos - q.beg.pos) % 3 == 2: # in-frame
        nuc_mutation_del_coding_inframe(args, q, tpt, r)
    else:   # frame-shift
        nuc_mutation_del_coding_frameshift(args, q, tpt, r)

def _nuc_mutation_del_intronic(args, q, tpt, r):

    codon = tpt.cpos2codon((q.beg.pos+2)/3)
    if not codon: return None
    i = q.beg.pos - (codon.index-1)*3 - 1
    if tpt.strand == '-': i = 2-i
    if tpt.strand == '+':
        r.gnuc_beg = codon.locs[i] + q.beg.tpos
        r.gnuc_end = codon.locs[i] + q.end.tpos
        if q.beg.tpos > 0: j = i+1
        elif q.beg.tpos < 0: j = i
    else:
        r.gnuc_beg = codon.locs[i] - q.end.tpos
        r.gnuc_end = codon.locs[i] - q.beg.tpos
        if q.beg.tpos > 0: j = i
        elif q.beg.tpos < 0: j = i+1
    if j>0 and j<3 and codon.locs[j-1]+1 == codon.locs[j]:
        raise IncompatibleTranscriptError('Codon does not contain exon boundary')

    pl = []
    s = '-'.join(map(str, codon.locs[:j]))
    if s: pl.append(s)
    pl.append('(%d)-(%d)' % (r.gnuc_beg, r.gnuc_end))
    s = '-'.join(map(str, codon.locs[j:]))
    if s: pl.append(s)
    r.pos = '%s:%s' % (tpt.chrm, '-'.join(pl))

    r.reg = '%s (%s intronic)' % (tpt.gene.name, tpt.strand)
    r.gnuc_range = '%d_%ddel' % (r.gnuc_beg, r.gnuc_end)

    r.refdelseq = faidx.refgenome.fetch_sequence(tpt.chrm, r.gnuc_beg, r.gnuc_end)
    r.natdelseq = r.refdelseq if tpt.strand == '+' else reverse_complement(r.refdelseq)
    if q.delseq and r.natdelseq != q.delseq:
        raise IncompatibleTranscriptError()

    if q.gnuc_beg == r.gnuc_end:
        r.tnuc_range = '%sdel' % (q.beg, )
    r.tnuc_range = '%s_%sdel' % (q.beg, q.end)

    r.info = 'RefDelSeq=%s;NatDelSeq=%s' % (r.refdelseq, r.natdelseq)

def nuc_mutation_del_intronic(args, q, tpt, r):

    """ deletion occurs entirely in non-coding region
    need only find genomic location of the deletion """
    if q.beg.pos == q.end.pos+43241: # with respect to one exome boundary
        _nuc_mutation_del_intronic(args, q, tpt, r)
    else:
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

            r_coding = Record('del')
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

def nuc_mutation_del_mix(args, q, tpt, r):

    q_coding = copy(q)
    q_intronic = copy(q)
    if q.beg.tpos < 0 and q.end.tpos == 0:
        q_intronic.beg = q.beg
        q_intronic.end = Pos()
        q_intronic.end.pos = q.beg.pos
        q_intronic.end.tpos = -1
        q_coding.delseq = q.delseq[-q.beg.tpos:]
        q_intronic.delseq = q.delseq[:-q.beg.tpos]
    elif q.beg.tpos == 0 and q.end.tpos > 0:
        q_intronic.beg = Pos()
        q_intronic.beg.pos = q.end.pos
        q_intronic.beg.tpos = 1
        q_intronic.end = q.end
        q_coding.delseq = q.delseq[:-q.end.tpos]
        q_intronic.delseq = q.delseq[-q.end.tpos:]

    #print q_coding.delseq, q_intronic.delseq, q.delseq
    r_coding = Record('del')
    nuc_mutation_del_coding(args, q_coding, tpt, r_coding)
    r_intronic = Record('del')
    nuc_mutation_del_intronic(args, q_intronic, tpt, r_intronic)

    # combine record intronic record with coding record
    #r_intronic.format(q.op)
    #r_coding.format(q.op)
    #print r_intronic.gnuc_beg, r_intronic.gnuc_end
    #print r_coding.gnuc_beg, r_coding.gnuc_end
    if r_intronic.gnuc_beg == r_coding.gnuc_end + 1:
        r.gnuc_beg = r_coding.gnuc_beg
        r.gnuc_end = r_intronic.gnuc_end
        r.refdelseq = r_coding.refdelseq + r_intronic.refdelseq
    elif r_intronic.gnuc_end + 1 == r_coding.gnuc_beg:
        """ intronic -- coding """
        r.gnuc_beg = r_intronic.gnuc_beg
        r.gnuc_end = r_coding.gnuc_end
        r.refdelseq = r_intronic.refdelseq + r_coding.refdelseq
    r.natdelseq = r.refdelseq if tpt.strand == '+' else reverse_complement(r.refdelseq)
    r.gnuc_range = '%d_%ddel' % (r.gnuc_beg, r.gnuc_end)
    r.tnuc_range = '%s_%sdel' % (q.beg, q.end)
    r.taa_range = r_coding.taa_range
    r.pos = '%s:%d-%d' % (tpt.chrm, r.gnuc_beg, r.gnuc_end)
    r.reg = '%s (%s coding & intronic)' % (tpt.gene.name, tpt.strand)
    r.info = 'RefDelSeq=%s;NatDelSeq=%s' % (r.refdelseq, r.natdelseq)

def nuc_mutation_del(args, q, tpt):

    if q.tpt and tpt.name != q.tpt: raise IncompatibleTranscriptError("Transcript name unmatched")
    tpt.ensure_seq()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name
    r.muttype = 'del'

    nuc_mutation_del_intronic(args, q, tpt, r)
    # if q.beg.tpos == 0 and q.end.tpos == 0:
    #     nuc_mutation_del_coding(args, q, tpt, r)
    # elif q.beg.tpos != 0 and q.end.tpos != 0:
    #     nuc_mutation_del_intronic(args, q, tpt, r)
    # else:
    #     # one of the deletion start and end is in coding, the other in non-coding
    #     nuc_mutation_del_mix(args, q, tpt, r)

    return r

def _core_annotate_nuc_del(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = nuc_mutation_del(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
        found = True
        r.format(q.op)

    if not found:
        r = Record()
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

    return
