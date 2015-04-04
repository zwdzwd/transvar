from transcripts import *
from record import *
from describe import *
# from anno_reg import __annotate_reg_intergenic

def taa_set_ins(r, t, index, taa_insseq):
    i1r, taa_insseq_r = t.taa_roll_right_ins(index, taa_insseq)
    try:
        r.taa_range = t.taa_ins_id(i1r, taa_insseq_r)
        i1l, taa_insseq_l = t.taa_roll_left_ins(index, taa_insseq)
        r.append_info('left_align_protein=p.%s' % t.taa_ins_id(i1l, taa_insseq_l))
        r.append_info('unalign_protein=p.%s' % t.taa_ins_id(index, taa_insseq))
    except IncompatibleTranscriptError:
        r.append_info("truncated_refseq_at_boundary")

def tnuc_set_ins(r, t, p1, tnuc_insseq):

    # note that intronic indel are NOT re-aligned,
    # because they are anchored with respect to exon boundaries.
    p1r, tnuc_insseq_r = t.tnuc_roll_right_ins(p1, tnuc_insseq)
    r.tnuc_range = '%d_%dins%s' % (p1r, p1r+1, tnuc_insseq_r)
    p1l, tnuc_insseq_l = t.tnuc_roll_left_ins(p1, tnuc_insseq)
    r.append_info('left_align_cDNA=c.%d_%dins%s' % (p1l, p1l+1, tnuc_insseq_l))
    r.append_info('unalign_cDNA=c.%s_%sins%s' % (p1, p1+1, tnuc_insseq))

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
                c, p = t.gpos2codon(q.pos, args)
                tnuc_insseq = q.insseq
            else:
                c, p = t.gpos2codon(q.pos+1, args)
                tnuc_insseq = reverse_complement(q.insseq)

            tnuc_set_ins(r, t, p.pos, tnuc_insseq)

            # TODO: check if insertion hits start codon

            # infer protein level mutation if in cds
            if reg.cds and not reg.splice: # this skips insertion that occur to sites next to donor or acceptor splicing site.
                if len(tnuc_insseq) % 3 == 0:
                    ins_gene_coding_inframe(t, r, c, p, tnuc_insseq)
                else:
                    ins_gene_coding_frameshift(t, r, c, p, tnuc_insseq)

        r.format(q.op)
