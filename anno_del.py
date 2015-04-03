from transcripts import *
from record import *
from describe import *
# from anno_reg import __annotate_reg_intergenic, _annotate_reg_gene_long_range

def taa_set_del(r, t, taa_beg, taa_end):

    i1r, i2r = t.taa_roll_right_del(taa_beg, taa_end)
    r.taa_range = t.taa_del_id(i1r, i2r)
    i1l, i2l = t.taa_roll_left_del(taa_beg, taa_end)
    r.append_info('left_align_protein=p.%s' % t.taa_del_id(i1l, i2l))
    r.append_info('unalign_protein=p.%s' % t.taa_del_id(taa_beg, taa_end))

def del_coding_inframe(args, cbeg, cend, pbeg, pend, q, t, r):

    if pbeg.pos % 3 == 1:       # in phase
        taa_set_del(r, t, cbeg.index, cend.index)
    else:                       # out-of-phase
        beg_codon_beg = cbeg.index*3-2
        end_codon_end = cend.index*3
        new_codon_seq = t.seq[beg_codon_beg-1:pbeg.pos-1] + \
                        t.seq[pend.pos:end_codon_end]
        if len(new_codon_seq) != 3:
            raise IncompatibleTranscriptError()
        r.taa_alt = codon2aa(new_codon_seq)
        tnuc_delseq = t.seq[beg_codon_beg-1:end_codon_end]
        taa_delseq = translate_seq(tnuc_delseq)
        # if taa_delseq[-1] == '*':
        if r.taa_alt == taa_delseq[-1]:
            # G100_S200delinsS becomes a pure deletion G100_D199del
            taa_set_del(r, t, cbeg.index, cend.index-1)
        elif r.taa_alt == taa_delseq[0]:
            # S100_G200delinsS becomes a pure deletion D101_G200del
            taa_set_del(r, t, cbeg.index+1, cend.index)
        else:
            r.taa_range = '%s%d_%s%ddelins%s' % (
                taa_delseq[0], cbeg.index,
                taa_delseq[-1], cend.index, r.taa_alt)

def del_coding_frameshift(args, cbeg, cend, pbeg, pend, q, t, r):

    """ assume frame-shift does not affect splicing """
    r.natdelseq = t.seq[pbeg.pos-1:pend.pos]
    if q.delseq and r.natdelseq != q.delseq:
        raise IncompatibleTranscriptError()

    cbeg_beg = cbeg.index*3 - 2
    old_seq = t.seq[cbeg_beg-1:]
    new_seq = t.seq[cbeg_beg-1:pbeg.pos-1]+t.seq[pend.pos:]
    if not old_seq:
        raise IncompatibleTranscriptError()

    ret = t.extend_taa_seq(cbeg.index, old_seq, new_seq)
    if ret:
        taa_pos, taa_ref, taa_alt, termlen = ret
        r.taa_range = '%s%d%sfs*%s' % (taa_ref, taa_pos, taa_alt, termlen)
    else:
        r.taa_range = '(=)'


def _annotate_del(args, q, db):

    normalize_reg(q)
    # for r in __annotate_del(args, q, db):
    #     r.format(q.op)

    for reg in describe(args, q, db):

        r = Record()
        r.reg = reg
        r.chrm = q.tok

        # right-align
        gnuc_beg_r, gnuc_end_r = gnuc_roll_right_del(q.tok, q.beg, q.end)
        r.gnuc_range = gnuc_del_id(q.tok, gnuc_beg_r, gnuc_end_r)

        # left-align
        gnuc_beg_l, gnuc_end_l = gnuc_roll_left_del(q.tok, q.beg, q.end)
        r.append_info('left_align_gDNA=g.%s' % gnuc_del_id(q.tok, gnuc_beg_l, gnuc_end_l))
        r.append_info('unaligned_gDNA=g.%s' % gnuc_del_id(q.tok, q.beg, q.end))

        if hasattr(reg, 't'):

            r.tname = reg.t.format()
            r.gene = reg.t.gene.name
            r.strand = reg.t.strand
            
            # whole gene deletion
            if q.end > t.cds_end-2 and q.beg < t.cds_beg + 2:
                r.append_info('whole_gene_deletion')

            # loss of start codon
            # TransVar took a simplistic approach, as long as
            # the deletion hit start codon, annotation is labeled as a start loss
            if q.beg <= t.cds_beg + 2 and q.end >= t.cds_beg:
                r.append_info('start_loss')

            # loss of stop codon
            if q.beg <= t.cds_end and q.end >= t.cds_end - 2:
                r.append_info('stop_loss')

            
            if t.strand == '+':
                c1, p1 = reg.t.gpos2codon(q.beg, args)
                c2, p2 = reg.t.gpos2codon(q.end, args)
            else:
                c1, p1 = reg.t.gpos2codon(q.end, args)
                c2, p2 = reg.t.gpos2codon(q.beg, args)

            # cDNA representation
            p1r, p2r = t.tnuc_roll_right_del(pbeg.pos, pend.pos)
            r.tnuc_range = t.tnuc_del_id(p1r, p2r)
            # left-aligned cDNA identifier
            p1l, p2l = t.tnuc_roll_left_del(pbeg.pos, pend.pos)
            r.append_info('left_align_cDNA=c.%s' % t.tnuc_del_id(p1l, p2l))
            r.append_info('unalign_cDNA=c.%s' % t.tnuc_del_id(pbeg.pos, pend.pos))

            if (q.end - q.beg + 1) % 3 == 0:
                del_coding_inframe(args, c1, c2, p1, p2, q, t, r)
            else:
                del_coding_frameshift(args, c1, c2, p1, p2, q, t, r)
                # t.ensure_seq()
                # alt_seq = t.seq[pbeg.included_plus()-1:pend.included_minus()]

        r.format(q.op)
