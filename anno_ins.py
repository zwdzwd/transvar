from transcripts import *
from record import *
from anno_reg import __annotate_reg_intergenic

def ins_gene_coding_inframe_inphase(t, r, c, p, insseq):

    tnuc_insseq = t.gnuc_seq2tnuc(q.insseq)
    taa_insseq = translate_seq(tnuc_insseq)
    if taa_insseq[-1] == '*':
        return ins_gene_coding_frameshift(t, r, c, p, insseq)

    # a pure insertion
    c1 = t.cpos2codon((q.pos.pos+2)/3)
    c2 = t.cpos2codon((q.pos.pos+3)/3)
    r.taa_range = '%s%d_%s%dins%s' % (codon2aa(c1.seq), c1.index, 
                                      codon2aa(c2.seq), c2.index, 
                                      taa_insseq)

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
        taa_insseq = taa_altseq[1:]
        c = t.taa_left_align_insertion_g(c.index, taa_insseq)
        c2 = t.cpos2codon(c.index+1)
        r.taa_range = '%s%d_%s%dins%s' % (
            c.aa(), c.index, c2.aa(), c2.index, taa_insseq)
    elif taa_ref == taa_altseq[-1]:
        # SdelinsHS becomes a pure insertion
        # [codon_before]_[current_codon]insH
        taa_insseq = taa_altseq[:-1]
        c = t.taa_left_align_insertion_g(c.index-1, taa_insseq)
        c2 = t.cpos2codon(c.index+1)
        r.taa_range = '%s%d_%s%dins%s' % (
            c.aa(), c.index, c2.aa(), c2.index, taa_insseq)
    else:
        r.taa_range = '%s%ddelins%s' % (taa_ref, c.index, taa_altseq)

def ins_gene_coding_inframe(t, r, c, p, insseq):

    if p.pos % 3 == 0:
        ins_gene_coding_inframe_inphase(t, r, c, p, insseq)
    else:
        ins_gene_coding_inframe_outphase(t, r, c, p, insseq)

def ins_gene_coding_frameshift(t, r, c, p, insseq):

    cbeg_beg = c.index * 3 - 2
    old_seq = tpt.seq[cbeg_beg-1:]
    new_seq = tpt.seq[cbeg_beg-1:p.pos]+insseq+tpt.seq[p.pos:]
    ret = extend_taa_seq(c.index, old_seq, new_seq, t)
    if ret:
        taa_pos, taa_ref, taa_alt, termlen = ret
        r.taa_range = '%s%d%sfs*%s' % (taa_ref, taa_pos, taa_alt, termlen)
    else:
        r.taa_range = '(=)'

def _ins_gene(args, q, t):

    r = Record()
    r.chrm = t.chrm
    r.tname = t.name
    r.pos = '%d' % q.pos
    r.gnuc_range = '%dins%s' % (q.pos, q.insseq)

    # c1 and c2 might be the same, rg1 and rg2 might be the same
    if t.strand == '+':
        c1, p1, rg1 = t.gpos2codon(q.tok, q.pos)
        c2, p2, rg2 = t.gpos2codon(q.tok, q.pos+1)
        insseq = q.insseq
    else:
        c1, p1, rg1 = t.gpos2codon(q.tok, q.pos+1)
        c2, p2, rg2 = t.gpos2codon(q.tok, q.pos)
        insseq = reverse_complement(q.insseq)

    r.tnuc_range = '%s_%sins%s' % (p1, p2, insseq)
    r.reg = '%s (%s, %s)' % (t.gene.name, t.strand, rg1.format())

    # TODO: if insertion hits start codon
    # infer protein level mutation if in cds
    if rg1.cds and not rg1.splice: # this skips insertion that occur to sites next to donor or acceptor splicing site.
        if len(insseq) % 3 == 0:
            ins_gene_coding_inframe(t, r, c1, p1, insseq)
        else:
            ins_gene_coding_frameshift(t, r, c1, p1, insseq)

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



