import faidx
from transcripts import *
from record import *
from describe import *
# from anno_reg import __annotate_reg_intergenic, _annotate_reg_gene_long_range

def _annotate_del(args, q, db):

    normalize_reg(q)
    # for r in __annotate_del(args, q, db):
    #     r.format(q.op)

    gnuc_delseq = faidx.getseq(q.tok, q.beg, q.end)
    warning = None
    if q.delseq and q.delseq != gnuc_delseq:
        warning = "invalid_deletion_seq_%s_(expect_%s)" % (gnuc_delseq, q.delseq)
        err_print("Warning: %s invalid deletion sequence %s (expect %s), maybe wrong reference?" % (q.op, gnuc_delseq, q.delseq))

    # right-align
    gnuc_beg_r, gnuc_end_r = gnuc_roll_right_del(q.tok, q.beg, q.end)
    gnuc_delseq_r = faidx.getseq(q.tok, gnuc_beg_r, gnuc_end_r)

    # left-align
    gnuc_beg_l, gnuc_end_l = gnuc_roll_left_del(q.tok, q.beg, q.end)
    gnuc_delseq_l = faidx.getseq(q.tok, gnuc_beg_l, gnuc_end_l)

    for reg in describe(args, q, db):

        r = Record()
        r.reg = reg
        r.chrm = q.tok
        if warning is not None:
            r.append_info(warning)

        r.gnuc_range = gnuc_del_id(q.tok, gnuc_beg_r, gnuc_end_r)
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
                c1, p1 = t.gpos2codon(q.beg)
                c2, p2 = t.gpos2codon(q.end)
                tnuc_delseq = gnuc_delseq
                c1l, p1l = t.gpos2codon(gnuc_beg_l)
                c2l, p2l = t.gpos2codon(gnuc_end_l)
                tnuc_delseq_l = gnuc_delseq_l
                c1r, p1r = t.gpos2codon(gnuc_beg_r)
                c2r, p2r = t.gpos2codon(gnuc_end_r)
                tnuc_delseq_r = gnuc_delseq_r
            else:
                c1, p1 = t.gpos2codon(q.end)
                c2, p2 = t.gpos2codon(q.beg)
                tnuc_delseq = reverse_complement(gnuc_delseq)
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
            r.append_info('unalign_cDNA=c.%s' % tnuc_del_id(p1.pos, p2.pos, tnuc_delseq))

            if t.transcript_type == 'protein_coding' and not same_intron(p1, p2):

                expt = r.set_splice('lost')
                if not expt:
                    if (q.end - q.beg + 1) % 3 == 0:
                        del_coding_inframe(args, c1, c2, p1, p2, t, r)
                    else:
                        del_coding_frameshift(args, c1, c2, p1, p2, t, r)
                        # t.ensure_seq()
                        # alt_seq = t.seq[pbeg.included_plus()-1:pend.included_minus()]

        r.format(q.op)
