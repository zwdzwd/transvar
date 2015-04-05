from transcripts import *
from record import *
from describe import *
# from anno_reg import __annotate_reg_intergenic

def ins_gene_coding_inframe_inphase(t, r, c, p, insseq):

    taa_insseq = translate_seq(insseq)
    if taa_insseq[-1] == '*':
        return ins_gene_coding_frameshift(t, r, c, p, insseq)

    # a pure insertion
    taa_set_ins(r, t, c.index, taa_insseq)
    

def ins_gene_coding_inframe_outphase(t, r, c, p, insseq):

    """ one codon *** split into *XX and X** """
    codon_beg = c.index*3-2
    codon_end = c.index*3
    codon_subseq1 = t.seq[codon_beg-1:p.pos]
    codon_subseq2 = t.seq[p.pos:codon_end]

    try:
        alt_seq = codon_subseq1+insseq+codon_subseq2
        taa_altseq = translate_seq(alt_seq)
        if taa_altseq[-1] == '*':
            return ins_gene_coding_frameshift(t, r, c, p, alt_seq)

        taa_ref = codon2aa(c.seq)
        if taa_ref == taa_altseq[0]:
            # SdelinsSH becomes a pure insertion
            # [current_codon]_[codon_after]insH
            taa_set_ins(r, t, c.index, taa_altseq[1:])

        elif taa_ref == taa_altseq[-1]:
            # SdelinsHS becomes a pure insertion
            # [codon_before]_[current_codon]insH
            taa_set_ins(r, t, c.index-1, taa_altseq[:-1])
        else:
            r.taa_range = '%s%ddelins%s' % (taa_ref, c.index, taa_altseq)
    except IncompatibleTranscriptError:
        r.append_info("truncated_refseq_at_boundary")


def ins_gene_coding_inframe(t, r, c, p, insseq):

    if p.pos % 3 == 0:
        ins_gene_coding_inframe_inphase(t, r, c, p, insseq)
    else:
        ins_gene_coding_inframe_outphase(t, r, c, p, insseq)

def ins_gene_coding_frameshift(t, r, c, p, insseq):

    cbeg_beg = c.index * 3 - 2
    old_seq = t.seq[cbeg_beg-1:]
    new_seq = t.seq[cbeg_beg-1:p.pos]+insseq+t.seq[p.pos:]
    ret = t.extend_taa_seq(c.index, old_seq, new_seq)
    if ret:
        taa_pos, taa_ref, taa_alt, termlen = ret
        if taa_alt == '*':      # a substitution into stop codon
            r.taa_range = '%s%d*' % (taa_ref, taa_pos)
        else:
            r.taa_range = '%s%d%sfs*%s' % (taa_ref, taa_pos, taa_alt, termlen)
    else:
        r.taa_range = '(=)'

def _annotate_ins(args, q, db):

    for reg in describe(args, q, db):

        r = Record()
        r.reg = reg
        r.chrm = q.tok
        r.pos = q.pos

        if q.insseq:
            # right align
            pos_r, gnuc_insseq_r = gnuc_roll_right_ins(q.tok, q.pos, q.insseq)
            r.gnuc_range = '%dins%s' % (pos_r, gnuc_insseq_r)
            
            # left align
            pos_l, gnuc_insseq_l = gnuc_roll_left_ins(q.tok, q.pos, q.insseq)
            r.append_info('left_align_gDNA=g.%dins%s' % (pos_l, gnuc_insseq_l))
            r.append_info('unalign_gDNA=g.%dins%s' % (q.pos, q.insseq))
        else:
            r.gnuc_range = '%dins' % (q.pos,)

        if hasattr(reg, 't'):

            t = reg.t

            r.tname = t.format()
            r.gene = t.gene.name
            r.strand = t.strand

            if t.strand == '+':
                c1, p1 = t.gpos2codon(q.pos)
                c2, p2 = t.gpos2codon(q.pos+1)
                tnuc_insseq = q.insseq
                c1l, p1l = t.gpos2codon(pos_l)
                c2l, p2l = t.gpos2codon(pos_l+1)
                tnuc_insseq_l = gnuc_insseq_l
                c1r, p1r = t.gpos2codon(pos_r)
                c2r, p2r = t.gpos2codon(pos_r+1)
                tnuc_insseq_r = gnuc_insseq_r
            else:
                c1, p1 = t.gpos2codon(q.pos+1)
                c2, p2 = t.gpos2codon(q.pos)
                tnuc_insseq = reverse_complement(q.insseq)
                c1l, p1l = t.gpos2codon(pos_r+1)
                c2l, p2l = t.gpos2codon(pos_r)
                tnuc_insseq_l = reverse_complement(gnuc_insseq_r)
                c1r, p1r = t.gpos2codon(pos_l+1)
                c2r, p2r = t.gpos2codon(pos_l)
                tnuc_insseq_r = reverse_complement(gnuc_insseq_l)

                
            r.tnuc_range = '%s_%sins%s' % (p1r, p2r, tnuc_insseq_r)
            r.append_info('left_align_cDNA=c.%s_%sins%s' % (p1l, p2l, tnuc_insseq_l))
            r.append_info('unalign_cDNA=c.%s_%sins%s' % (p1, p2, tnuc_insseq))

            # TODO: check if insertion hits start codon

            # infer protein level mutation if in cds
            if reg.cds and not reg.splice: # this skips insertion that occur to sites next to donor or acceptor splicing site.
                if len(tnuc_insseq) % 3 == 0:
                    ins_gene_coding_inframe(t, r, c1, p1, tnuc_insseq)
                else:
                    ins_gene_coding_frameshift(t, r, c1, p1, tnuc_insseq)

        r.format(q.op)
