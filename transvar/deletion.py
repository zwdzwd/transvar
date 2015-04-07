from transcripts import *
from utils import *
from record import *
from copy import copy
from describe import *

# TODO: refactor left-right align

def _annotate_deletion_cdna(args, q, r, t, db):

    # check exon boundary
    t.check_exon_boundary(q.beg)
    t.check_exon_boundary(q.end)

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
    r.pos = '%d-%d' % (gnuc_beg_r, gnuc_end_r)

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

def annotate_deletion_cdna(args, q, tpts, db):

    found = False
    for t in tpts:
        if q.tpt and t.name != q.tpt:
            raise IncompatibleTranscriptError("Transcript name unmatched")

        r = Record()
        r.chrm = t.chrm
        r.tname = t.format()
        r.gene = t.gene.name
        r.strand = t.strand

        try:
            t.ensure_seq()
            _annotate_deletion_cdna(args, q, r, t, db)
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

def annotate_deletion_protein(args, q, tpts, db):

    found = False
    for t in tpts:
        try:
            if q.tpt and t.name != q.tpt:
                raise IncompatibleTranscriptError("Transcript name unmatched")
            t.ensure_seq()

            r = Record()
            r.chrm = t.chrm
            r.tname = t.format()
            r.gene = t.gene.name
            r.strand = t.strand

            if q.beg*3 > len(t) or q.end*3 > len(t):
                raise IncompatibleTranscriptError('codon nonexistent')

            if q.delseq and t.taa_range2aa_seq(q.beg, q.end) != q.delseq:
                raise IncompatibleTranscriptError('unmatched reference')

            tnuc_beg = q.beg*3-2
            tnuc_end = q.end*3
            gnuc_beg, gnuc_end = t.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
            r.tnuc_range = '%d_%d' % (tnuc_beg, tnuc_end)
            r.gnuc_range = '%d_%d' % (gnuc_beg, gnuc_end)
            r.pos = '%d-%d' % (gnuc_beg, gnuc_end)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue
        except UnknownChromosomeError:
            continue

        taa_set_del(r, t, q.beg, q.end)
        r.reg = describe_genic(args, t.chrm, gnuc_beg, gnuc_end, t, db)
        r.append_info('imprecise')
        r.format(q.op)
        found = True

    if not found:
        r = Record()
        r.taa_range = t.taa_del_id(q.beg, q.end)
        r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))

        r.format(q.op)

def annotate_deletion_gdna(args, q, db):

    normalize_reg(q)

    gnuc_delseq = faidx.getseq(q.tok, q.beg, q.end)
    warning = None
    if q.delseq and q.delseq != gnuc_delseq:
        warning = "invalid_deletion_seq_%s_(expect_%s)" % (gnuc_delseq, q.delseq)
        err_warn("%s invalid deletion sequence %s (expect %s), maybe wrong reference?" % (q.op, gnuc_delseq, q.delseq))

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
        r.pos = '%d-%d' % (gnuc_beg_r, gnuc_end_r)
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

        r.format(q.op)


### add taa feature in deletion ###

def taa_set_del(r, t, taa_beg, taa_end):

    i1r, i2r = t.taa_roll_right_del(taa_beg, taa_end)
    r.taa_range = t.taa_del_id(i1r, i2r)
    i1l, i2l = t.taa_roll_left_del(taa_beg, taa_end)
    r.append_info('left_align_protein=p.%s' % t.taa_del_id(i1l, i2l))
    r.append_info('unalign_protein=p.%s' % t.taa_del_id(taa_beg, taa_end))

def del_coding_inframe(args, c1, c2, p1, p2, t, r):

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
            err_print(t.gene.name+'\t'+t.transcript_type)
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

def del_coding_frameshift(args, cbeg, cend, pbeg, pend, t, r):

    """ assume frame-shift does not affect splicing """
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

