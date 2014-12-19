from transcripts import *
from record import *
from anno_reg import __annotate_reg_intergenic

def taa_set_ins(r, t, index, taa_insseq):

    i1r, taa_insseq_r = t.taa_roll_right_ins(index, taa_insseq)
    r.taa_range = t.taa_ins_id(i1r, taa_insseq_r)

    i1l, taa_insseq_l = t.taa_roll_left_ins(index, taa_insseq)
    r.append_info('LEFTALNP=p.%s' % t.taa_ins_id(i1l, taa_insseq_l))
    r.append_info('UALNP=p.%s' % t.taa_ins_id(index, taa_insseq))

def tnuc_set_ins(r, t, p1, tnuc_insseq):

    # note that intronic indel are NOT re-aligned,
    # because they are anchored with respect to exon boundaries.
    p1r, tnuc_insseq_r = t.tnuc_roll_right_ins(p1, tnuc_insseq)
    r.tnuc_range = '%d_%dins%s' % (p1r, p1r+1, tnuc_insseq_r)
    p1l, tnuc_insseq_l = t.tnuc_roll_left_ins(p1, tnuc_insseq)
    r.append_info('LEFTALNC=c.%d_%dins%s' % (p1l, p1l+1, tnuc_insseq_l))
    r.append_info('UALNC=c.%s_%sins%s' % (p1, p1+1, tnuc_insseq))

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
    

def ins_gene_coding_inframe(t, r, c, p, insseq):

    if p.pos % 3 == 0:
        ins_gene_coding_inframe_inphase(t, r, c, p, insseq)
    else:
        ins_gene_coding_inframe_outphase(t, r, c, p, insseq)

def ins_gene_coding_frameshift(t, r, c, p, insseq):

    cbeg_beg = c.index * 3 - 2
    old_seq = t.seq[cbeg_beg-1:]
    new_seq = t.seq[cbeg_beg-1:p.pos]+insseq+t.seq[p.pos:]
    ret = extend_taa_seq(c.index, old_seq, new_seq, t)
    if ret:
        taa_pos, taa_ref, taa_alt, termlen = ret
        if taa_alt == '*':      # a substitution into stop codon
            r.taa_range = '%s%d*' % (taa_ref, taa_pos)
        else:
            r.taa_range = '%s%d%sfs*%s' % (taa_ref, taa_pos, taa_alt, termlen)
    else:
        r.taa_range = '(=)'

def _ins_gene(args, q, t):

    r = Record()
    r.chrm = t.chrm
    r.tname = t.name
    r.pos = '%d' % q.pos

    if q.insseq:
        pos_r, gnuc_insseq_r = gnuc_roll_right_ins(t.chrm, q.pos, q.insseq)
        r.gnuc_range = '%dins%s' % (pos_r, gnuc_insseq_r)
        pos_l, gnuc_insseq_l = gnuc_roll_left_ins(t.chrm, q.pos, q.insseq)
        r.append_info('LEFTALNG=g.%dins%s' % (pos_l, gnuc_insseq_l))
        r.append_info('UALNG=g.%dins%s' % (q.pos, q.insseq))
    else:
        r.gnuc_range = '%dins%s' % (q.pos, q.insseq)


    # c1 and c2 might be the same, rg1 and rg2 might be the same
    if t.strand == '+':
        c1, p1, rg1 = t.gpos2codon(q.pos)
        tnuc_insseq = q.insseq
    else:
        c1, p1, rg1 = t.gpos2codon(q.pos+1)
        tnuc_insseq = reverse_complement(q.insseq)
    r.reg = '%s (%s, %s)' % (t.gene.name, t.strand, rg1.format())

    # assert p1.tpos == p2.tpos == 0
    tnuc_set_ins(r, t, p1.pos, tnuc_insseq)

    # TODO: if insertion hits start codon
    # infer protein level mutation if in cds
    if rg1.cds and not rg1.splice: # this skips insertion that occur to sites next to donor or acceptor splicing site.
        if len(tnuc_insseq) % 3 == 0:
            ins_gene_coding_inframe(t, r, c1, p1, tnuc_insseq)
        else:
            ins_gene_coding_frameshift(t, r, c1, p1, tnuc_insseq)

    return r

def __annotate_ins(args, q, db):

    tpts = [t for t in db.get_transcripts(q.tok, q.pos)]
    if tpts:
        if args.longest:
            tpts.sort(key=lambda t:len(t), reverse=True)
            tpts = tpts[:1]
        for t in tpts:
            yield _ins_gene(args, q, t)
    else:
        yield __annotate_reg_intergenic(args, db, q.tok, q.pos, q.pos)

def _annotate_ins(args, q, db):
    for r in __annotate_ins(args, q, db):
        r.format(q.op)



