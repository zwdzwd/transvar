import faidx
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

def del_coding_inframe(args, c1, c2, p1, p2, q, t, r):

    if p1.pos % 3 == 1:       # in phase
        taa_set_del(r, t, c1.index, c2.index)
    else:                       # out-of-phase

        if len(c1.seq) != 3 or len(c2.seq) != 3:
            if len(t.seq) % 3 != 0:
                r.append_info("truncated_refseq_at_boundary_(start_codon_seq_%s_and_end_codon_seq_%s)" % (c1.seq, c2.seq))
                return
            raise IncompatibleTranscriptError()
        
        beg_codon_beg = c1.index*3-2
        end_codon_end = c2.index*3
        new_codon_seq = t.seq[beg_codon_beg-1:p1.pos-1] + \
                        t.seq[p2.pos:end_codon_end]

        if len(new_codon_seq) != 3:
            sys.stderr.write(q.op+'\t'+t.gene.name+'\t'+t.transcript_type+'\n')
            err_print(p1)
            err_print(p2)
            err_print(len(t.seq))
            err_print(c1.seq)
            err_print(c2.seq)
            err_print(len(t.seq) % 3)
            err_print(t.seq[-10:])
            if (len(t.seq)%3 != 0):
                r.append_info('truncated_refseq_at_boundary_(codon_seq_%s)' % c1.seq)
            raise IncompatibleTranscriptError('new_codon_seq: %s' % new_codon_seq)

        r.taa_alt = codon2aa(new_codon_seq)
        tnuc_delseq = t.seq[beg_codon_beg-1:end_codon_end]
        taa_delseq = translate_seq(tnuc_delseq)
        # if taa_delseq[-1] == '*':
        if r.taa_alt == taa_delseq[-1]:
            # G100_S200delinsS becomes a pure deletion G100_D199del
            taa_set_del(r, t, c1.index, c2.index-1)
        elif r.taa_alt == taa_delseq[0]:
            # S100_G200delinsS becomes a pure deletion D101_G200del
            taa_set_del(r, t, c1.index+1, c2.index)
        else:
            r.taa_range = '%s%d_%s%ddelins%s' % (
                taa_delseq[0], c1.index,
                taa_delseq[-1], c2.index, r.taa_alt)

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

        gnuc_delseq = faidx.getseq(q.tok, q.beg, q.end)

        # right-align
        gnuc_beg_r, gnuc_end_r = gnuc_roll_right_del(q.tok, q.beg, q.end)
        gnuc_delseq_r = faidx.getseq(q.tok, gnuc_beg_r, gnuc_end_r)
        r.gnuc_range = gnuc_del_id(q.tok, gnuc_beg_r, gnuc_end_r)

        # left-align
        gnuc_beg_l, gnuc_end_l = gnuc_roll_left_del(q.tok, q.beg, q.end)
        gnuc_delseq_l = faidx.getseq(q.tok, gnuc_beg_l, gnuc_end_l)
        r.append_info('left_align_gDNA=g.%s' % gnuc_del_id(q.tok, gnuc_beg_l, gnuc_end_l))
        r.append_info('unaligned_gDNA=g.%s' % gnuc_del_id(q.tok, q.beg, q.end))

        if hasattr(reg, 't'):

            t = reg.t
            r.tname = t.format()
            r.gene = t.gene.name
            r.strand = t.strand
            
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
                c1, p1 = t.gpos2codon(q.beg, args)
                c2, p2 = t.gpos2codon(q.end, args)
                tnuc_delseq = gnuc_delseq
                c1l, p1l = t.gpos2codon(gnuc_beg_l, args)
                c2l, p2l = t.gpos2codon(gnuc_end_l, args)
                tnuc_delseq_l = gnuc_delseq_l
                c1r, p1r = t.gpos2codon(gnuc_beg_r, args)
                c2r, p2r = t.gpos2codon(gnuc_end_r, args)
                tnuc_delseq_r = gnuc_delseq_r
            else:
                c1, p1 = t.gpos2codon(q.end, args)
                c2, p2 = t.gpos2codon(q.beg, args)
                tnuc_delseq = reverse_complement(gnuc_delseq)
                c1l, p1l = t.gpos2codon(gnuc_end_r, args)
                c2l, p2l = t.gpos2codon(gnuc_beg_r, args)
                tnuc_delseq_l = reverse_complement(gnuc_delseq_r)
                c1r, p1r = t.gpos2codon(gnuc_end_l, args)
                c2r, p2r = t.gpos2codon(gnuc_beg_l, args)
                tnuc_delseq_r = reverse_complement(gnuc_delseq_l)

            # cDNA representation
            # right align
            r.tnuc_range = t.tnuc_del_id(p1r, p2r, tnuc_delseq_r)
            # left align
            r.append_info('left_align_cDNA=c.%s' % t.tnuc_del_id(p1l, p2l, tnuc_delseq_l))
            r.append_info('unalign_cDNA=c.%s' % t.tnuc_del_id(p1.pos, p2.pos, tnuc_delseq))

            if t.transcript_type == 'protein_coding' and not same_intron(p1, p2):
                expt = False
                if hasattr(reg, 'splice_donors') and reg.splice_donors:
                    r.append_info('donor_splice_site_on_exon_%d_lost' % reg.splice_donors[0])
                    expt = True
                    
                if hasattr(reg, 'splice_acceptors') and reg.splice_acceptors:
                    r.append_info('acceptor_splice_site_on_exon_%d_lost' % reg.splice_acceptors[0])
                    expt = True

                if hasattr(reg, 'splice_both') and reg.splice_both:
                    r.append_info('whole_exon_[%s]_lost' % ','.join(map(str,reg.splice_both)))
                    expt = True

                if hasattr(reg, 'cross_start') and reg.cross_start:
                    r.append_info('cds_start_(%s:%d)_lost' % (t.chrm, t.cds_beg))
                    expt = True

                if hasattr(reg, 'cross_end') and reg.cross_end:
                    r.append_info('cds_end_(%s:%d)_lost' % (t.chrm, t.cds_end))
                    expt = True

                # print q.beg, t.cds_end, q.end
                # print hasattr(reg, 'cross_end')
                # if hasattr(reg, 'cross_end'):
                #     print reg.cross_end
                if not expt:
                    if (q.end - q.beg + 1) % 3 == 0:
                        del_coding_inframe(args, c1, c2, p1, p2, q, t, r)
                    else:
                        del_coding_frameshift(args, c1, c2, p1, p2, q, t, r)
                        # t.ensure_seq()
                        # alt_seq = t.seq[pbeg.included_plus()-1:pend.included_minus()]

        r.format(q.op)
