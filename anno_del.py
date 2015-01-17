from transcripts import *
from record import *
from anno_reg import __annotate_reg_intergenic, _annotate_reg_gene_long_range

def taa_set_del(r, t, taa_beg, taa_end):

    i1r, i2r = t.taa_roll_right_del(taa_beg, taa_end)
    r.taa_range = t.taa_del_id(i1r, i2r)
    i1l, i2l = t.taa_roll_left_del(taa_beg, taa_end)
    r.append_info('LEFTALNP=p.%s' % t.taa_del_id(i1l, i2l))
    r.append_info('UALNP=p.%s' % t.taa_del_id(taa_beg, taa_end))

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


def _annotate_del_single_gene(args, q, t):

    r = Record()
    r.chrm = t.chrm
    r.tname = t.name
    r.pos = '%d-%d' % (q.beg, q.end)

    # genomic annotation
    gnuc_beg_r, gnuc_end_r = gnuc_roll_right_del(q.tok, q.beg, q.end)
    r.gnuc_range = gnuc_del_id(q.tok, gnuc_beg_r, gnuc_end_r)
    gnuc_beg_l, gnuc_end_l = gnuc_roll_left_del(q.tok, q.beg, q.end)
    r.append_info('LEFTALNG=g.%s' % 
                  gnuc_del_id(q.tok, gnuc_beg_l, gnuc_end_l))
    r.append_info('UALNG=g.%s' % 
                  gnuc_del_id(q.tok, q.beg, q.end))

    if q.end > t.cds_end-2 and q.beg < t.cds_beg + 2:
        # whole gene gets deleted.
        r.reg = '%s (%s, WholeGeneDeletion)' % (t.gene.name, t.strand)
    elif q.beg <= t.cds_beg + 2 and q.end >= t.cds_beg:
        # loss of start codon
        # TransVar took a simplistic approach, as long as
        # the deletion hit start codon, annotation is labeled as a start loss
        r.reg = '%s (%s, StartLoss)' % (t.gene.name, t.strand)
    elif q.beg <= t.cds_end and q.end >= t.cds_end - 2:
        r.reg = '%s (%s, StopLoss)' % (t.gene.name, t.strand)
    else:                       # in CDS
        if t.strand == '+':
            cbeg, pbeg, rg_beg = t.gpos2codon(q.beg)
            cend, pend, rg_end = t.gpos2codon(q.end)
        else:
            cbeg, pbeg, rg_beg = t.gpos2codon(q.end)
            cend, pend, rg_end = t.gpos2codon(q.beg)

        # cDNA representation
        p1r, p2r = t.tnuc_roll_right_del(pbeg.pos, pend.pos)
        r.tnuc_range = t.tnuc_del_id(p1r, p2r)
        # left-aligned cDNA identifier
        p1l, p2l = t.tnuc_roll_left_del(pbeg.pos, pend.pos)
        r.append_info('LEFTALNC=c.%s' % t.tnuc_del_id(p1l, p2l))
        r.append_info('UALNC=c.%s' % t.tnuc_del_id(pbeg.pos, pend.pos))

        # r.append_info('BEGCodon=%s' % '-'.join(map(str, cbeg.locs)))
        # r.append_info('ENDCodon=%s' % '-'.join(map(str, cend.locs)))
        if rg_beg.format() == rg_end.format():
            r.append_info('REG=%s' % rg_beg.format())
        else:
            r.append_info('BEGREG=%s' % rg_beg.format())
            r.append_info('ENDREG=%s' % rg_end.format())

        # if beg and end are both in introns, still regarded as a
        # valid splicing

        # if one of the beg and end is in intron, the other in exon
        # not valid splicng
        if ((rg_beg.intronic and rg_end.exonic) or
            (rg_beg.exonic and rg_end.intronic)):
            r.reg = '%s (%s, Intronic;Exonic)' % (t.gene.name, t.strand)
            r.append_info('DisruptedSplicing')
        if (rg_beg.cds and rg_end.cds and
            rg_beg.exon == rg_end.exon): # coding sequence deletion
            r.reg = '%s (%s, Coding)' % (t.gene.name, t.strand)
            if (q.end - q.beg + 1) % 3 == 0: # in-frame deletion
                del_coding_inframe(args, cbeg, cend, pbeg, pend, q, t, r)
            else:               # frame-shift deletion
                del_coding_frameshift(args, cbeg, cend, pbeg, pend, q, t, r)
                t.ensure_seq()
                alt_seq = t.seq[pbeg.included_plus()-1:pend.included_minus()]

    return r

def _annotate_del_gene(args, q, db):

    tpts = [t for t in db.get_transcripts(q.tok, q.beg, q.end)]
    if tpts:
        if args.longest:
            tpts.sort(key=lambda t: len(t), reverse=True)
            tpts = tpts[:1]
        
        genes = list(set([t.gene for t in tpts]))
        if len(genes) == 1:
            for t in tpts:
                yield _annotate_del_single_gene(args, q, t)
        else:
            for r in _annotate_reg_gene_long_range(args, q, tpts, genes, db):
                yield r

def __annotate_del(args, q, db):
    # if annotation is in the coding region
    gene_found = False
    for r in _annotate_del_gene(args, q, db):
        yield r
        gene_found = True

    if not gene_found:
        yield __annotate_reg_intergenic(args, db, q.tok, q.beg, q.end)

def _annotate_del(args, q, db):

    normalize_reg(q)
    for r in __annotate_del(args, q, db):
        r.format(q.op)
