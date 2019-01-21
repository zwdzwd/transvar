"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Wanding Zhou, Tenghui Chen, Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from .transcripts import *
from .utils import *
from .record import *
from copy import copy
from .describe import *
from .proteinseqs import *

# TODO: refactor left-right align
class GNucDeletion():

    def __init__(self, chrm, gnuc_beg, gnuc_end):
        self.chrm = chrm

        # unaligned
        self.gnuc_beg = gnuc_beg
        self.gnuc_end = gnuc_end
        self.gnuc_delseq = faidx.getseq(chrm, gnuc_beg, gnuc_end)

        # right-aligned
        self.gnuc_beg_r, self.gnuc_end_r = gnuc_roll_right_del(chrm, gnuc_beg, gnuc_end)
        self.gnuc_delseq_r = faidx.getseq(chrm, self.gnuc_beg_r, self.gnuc_end_r) # TODO: this retrieval can be saved

        # left-aligned
        self.gnuc_beg_l, self.gnuc_end_l = gnuc_roll_left_del(chrm, gnuc_beg, gnuc_end)
        # TODO protect against self.gnuc_beg_l == 0
        seq = faidx.getseq(chrm, self.gnuc_beg_l-1, self.gnuc_end_l)
        self.gnuc_delseq_l = seq[1:]
        self.gnuc_leftbase_l = seq[0] # this is important for VCF output

        self.tnuc = False

    def compute_tnuc(self, t, p1=None, p2=None):

        self.tnuc = True
        if t.strand == '+':
            self.tnuc_delseq = self.gnuc_delseq
            if p1 is None:
                self.c1, self.p1 = t.gpos2codon(self.gnuc_beg)
            else:
                self.p1 = p1

            if p2 is None:
                self.c2, self.p2 = t.gpos2codon(self.gnuc_end)
            else:
                self.p2 = p2

            c1l, self.p1l = t.gpos2codon(self.gnuc_beg_l)
            c2l, self.p2l = t.gpos2codon(self.gnuc_end_l)
            self.tnuc_delseq_l = self.gnuc_delseq_l
            c1r, self.p1r = t.gpos2codon(self.gnuc_beg_r)
            c2r, self.p2r = t.gpos2codon(self.gnuc_end_r)
            self.tnuc_delseq_r = self.gnuc_delseq_r
        else:
            self.tnuc_delseq = reverse_complement(self.gnuc_delseq)
            if p1 is None:
                self.c1, self.p1 = t.gpos2codon(self.gnuc_end)
            else:
                self.p1 = p1

            if p2 is None:
                self.c2, self.p2 = t.gpos2codon(self.gnuc_beg)
            else:
                self.p2 = p2

            c1l, self.p1l = t.gpos2codon(self.gnuc_end_r)
            c2l, self.p2l = t.gpos2codon(self.gnuc_beg_r)
            self.tnuc_delseq_l = reverse_complement(self.gnuc_delseq_r)
            c1r, self.p1r = t.gpos2codon(self.gnuc_end_l)
            c2r, self.p2r = t.gpos2codon(self.gnuc_beg_l)
            self.tnuc_delseq_r = reverse_complement(self.gnuc_delseq_l)

    def set_record(self, r, args):
        r.pos = '%d-%d' % (self.gnuc_beg_r, self.gnuc_end_r)

        r.gnuc_range = gnuc_del_id(self.chrm, self.gnuc_beg_r, self.gnuc_end_r, args)

        # optional output
        if args.gseq:
            r.vcf_pos = self.gnuc_beg_l - 1
            r.vcf_ref = self.gnuc_leftbase_l + self.gnuc_delseq_l
            r.vcf_alt = self.gnuc_leftbase_l

        r.append_info('left_align_gDNA=g.%s' % gnuc_del_id(
            self.chrm, self.gnuc_beg_l, self.gnuc_end_l, args))
        r.append_info('unaligned_gDNA=g.%s' % gnuc_del_id(
            self.chrm, self.gnuc_beg, self.gnuc_end, args))

        if self.tnuc:
            r.tnuc_range = tnuc_del_id(self.p1r, self.p2r, args, self.tnuc_delseq_r)
            r.append_info('left_align_cDNA=c.%s' % tnuc_del_id(
                self.p1l, self.p2l, args, self.tnuc_delseq_l))
            r.append_info('unalign_cDNA=c.%s' % tnuc_del_id(
                self.p1, self.p2, args, self.tnuc_delseq))

def _annotate_deletion_cdna(args, q, r, t, db):

    # check exon boundary
    t.check_exon_boundary(q.beg)
    t.check_exon_boundary(q.end)

    _gnuc_beg = t.tnuc2gnuc(q.beg)
    _gnuc_end = t.tnuc2gnuc(q.end)

    gnuc_beg = min(_gnuc_beg, _gnuc_end)
    gnuc_end = max(_gnuc_beg, _gnuc_end)

    # set_deletion_id(args, t.chrm, gnuc_beg, gnuc_end, t, tnuc_delseq0=q.delseq)
    gnd = GNucDeletion(t.chrm, gnuc_beg, gnuc_end)
    gnd.compute_tnuc(t, p1=q.beg, p2=q.end)

    if q.delseq and q.delseq != gnd.tnuc_delseq and not args.ignore:
        raise IncompatibleTranscriptError()

    gnd.set_record(r, args)

    r.reg = describe_genic(args, t.chrm, gnuc_beg, gnuc_end, t, db)

    # if deletion affects coding region
    if t.transcript_type == 'protein_coding' and not same_intron(q.beg, q.end):
        if not r.set_splice('lost', csqn_action="Deletion"):
            c1, p1 = t.intronic_lean(q.beg, 'c_greater')
            c2, p2 = t.intronic_lean(q.end, 'c_smaller')

            if (gnuc_end - gnuc_beg + 1) % 3 == 0:
                del_coding_inframe(args, c1, c2, p1, p2, t, r)
            else:
                del_coding_frameshift(args, c1, c2, p1, p2, t, r)
    else:
        r.set_csqn_byreg("Deletion")

    return r

def annotate_deletion_cdna(args, q, tpts, db):

    # found = False
    records = []
    for t in tpts:
        if q.tpt and t.name != q.tpt:
            raise IncompatibleTranscriptError("Transcript name unmatched")
        # if args.strictversion and t.version != q.tpt.version:
        #     raise IncompatibleTranscriptError("Transcript name unmatched")

        r = Record(is_var=True)
        r.chrm = t.chrm
        r.tname = t.format()
        r.gene = t.gene_name
        r.strand = t.strand

        try:
            t.ensure_seq()
            _annotate_deletion_cdna(args, q, r, t, db)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue

        # found = True
        records.append(r)

    format_records(records, q.op, args)
    # format_one(r, rs, q.op, args)
    # format_all(rs, q.op, args)

    # if not found:
    #     r = Record(is_var=True)
    #     tnuc_del_id(q.beg, q.end, args, q.delseq)
    #     r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))
    #     r.format(q.op)

    return records

def annotate_deletion_protein(args, q, tpts, db):

    # found = False
    # rs = []
    records = []
    for t in tpts:
        try:
            if q.tpt and t.name != q.tpt:
                raise IncompatibleTranscriptError("Transcript name unmatched")
            t.ensure_seq()

            r = Record(is_var=True)
            r.chrm = t.chrm
            r.tname = t.format()
            r.gene = t.gene_name
            r.strand = t.strand

            if q.end*3 > t.cdslen():
                raise IncompatibleTranscriptError('codon nonexistent')

            if q.delseq and t.taa_range2aa_seq(q.beg, q.end) != q.delseq:
                raise IncompatibleTranscriptError('unmatched reference')

            tnuc_beg = q.beg*3-2
            tnuc_end = q.end*3
            gnuc_beg, gnuc_end = t.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
            gnd = GNucDeletion(t.chrm, gnuc_beg, gnuc_end)
            gnd.compute_tnuc(t)
            gnd.set_record(r, args)

            # r.tnuc_range = '%d_%ddel' % (tnuc_beg, tnuc_end)
            # r.gnuc_range = '%d_%ddel' % (gnuc_beg, gnuc_end)
            # r.pos = '%d-%d' % (gnuc_beg, gnuc_end)
            r.csqn.append("InFrameDeletion")
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue

        taa_set_del(r, t, q.beg, q.end, args)
        r.reg = describe_genic(args, t.chrm, gnuc_beg, gnuc_end, t, db)
        r.append_info('imprecise')
        # found = True
        records.append(r)
        # format_one(r, rs, q.op, args)

    # format_all(rs, q.op, args)
    format_records(records, q.op, args)
    # if not found:
    #     r = Record(is_var=True)
    #     r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))
    #     r.format(q.op)
    return records

def annotate_deletion_gdna(args, q, db):

    normalize_reg(q)

    gnd = GNucDeletion(q.tok, q.beg, q.end)

    warning = None
    if q.delseq and q.delseq != gnd.gnuc_delseq:
        warning = "invalid_deletion_seq_%s_(expect_%s)" % (gnd.gnuc_delseq, q.delseq)
        err_warn("%s invalid deletion sequence %s (expect %s), maybe wrong reference?" % (q.op, gnd.gnuc_delseq, q.delseq))

    # rs = []
    records = []
    for reg in describe(args, q, db):

        r = Record(is_var=True)
        r.reg = reg
        r.chrm = q.tok
        if warning is not None:
            r.append_info(warning)

        # r.gnuc_range = gnuc_del_id(q.tok, gnuc_beg_r, gnuc_end_r, args)
        # r.pos = '%d-%d' % (gnuc_beg_r, gnuc_end_r)
        # r.append_info('left_align_gDNA=g.%s' % gnuc_del_id(q.tok, gnuc_beg_l, gnuc_end_l, args))
        # r.append_info('unaligned_gDNA=g.%s' % gnuc_del_id(q.tok, q.beg, q.end, args))

        db.query_dbsnp_range(r, q.beg, q.end, '')
        if hasattr(reg, 't'):

            t = reg.t
            gnd.compute_tnuc(t)
            r.tname = t.format()
            r.gene = t.gene_name
            r.strand = t.strand

            # whole gene deletion
            if q.end > t.cds_end-2 and q.beg < t.cds_beg + 2:
                r.append_info('whole_gene_deletion')

            # # loss of start codon
            # # TransVar took a simplistic approach, as long as
            # # the deletion hit start codon, annotation is labeled as a start loss
            # if q.beg <= t.cds_beg + 2 and q.end >= t.cds_beg:
            #     r.append_info('start_loss')

            # # loss of stop codon
            # if q.beg <= t.cds_end and q.end >= t.cds_end - 2:
            #     r.append_info('stop_loss')

            gnd.set_record(r, args)
            if t.transcript_type == 'protein_coding' and not same_intron(gnd.p1, gnd.p2):
                if not r.set_splice('lost', csqn_action="Deletion"):
                    if (q.end - q.beg + 1) % 3 == 0:
                        del_coding_inframe(args, gnd.c1, gnd.c2, gnd.p1, gnd.p2, t, r)
                    else:
                        del_coding_frameshift(args, gnd.c1, gnd.c2, gnd.p1, gnd.p2, t, r)
            else:
                r.set_csqn_byreg("Deletion")
        else:
            gnd.set_record(r, args)
            r.set_csqn_byreg("Deletion")
        records.append(r)

    format_records(records, q.op, args)
    # format_one(r, rs, q.op, args)
    # format_all(rs, q.op, args)
    return records

### add taa feature in deletion ###

def taa_del_id(t, taa_beg, taa_end, args):

    if taa_beg == taa_end:
        s = '%s%ddel%s' % (aaf(t.cpos2aa(taa_beg), args), taa_beg, aaf(t.taa2aa(taa_beg), args))
    else:
        taa_del_len = taa_end - taa_beg + 1
        if taa_del_len > args.seqmax and args.seqmax >= 0:
            taa_delrep = str(taa_del_len)
        else:
            taa_delrep = aaf(t.taa_range2aa_seq(taa_beg, taa_end), args)
        s = '%s%d_%s%ddel%s' % (
            aaf(t.cpos2aa(taa_beg), args), taa_beg,
            aaf(t.cpos2aa(taa_end), args), taa_end,
            taa_delrep)

    return s

def taa_set_del(r, t, taa_beg, taa_end, args):

    i1r, i2r = t.taa_roll_right_del(taa_beg, taa_end)
    r.taa_range = taa_del_id(t, i1r, i2r, args)
    i1l, i2l = t.taa_roll_left_del(taa_beg, taa_end)
    r.append_info('left_align_protein=p.%s' %
                  taa_del_id(t, i1l, i2l, args))
    r.append_info('unalign_protein=p.%s' %
                  taa_del_id(t, taa_beg, taa_end, args))
    variant_protein_seq_del(r, t, args, i1r, i2r)

def del_coding_inframe(args, c1, c2, p1, p2, t, r):

    if p1.pos % 3 == 1:       # in phase
        r.csqn.append('InFrameDeletion')
        taa_set_del(r, t, c1.index, c2.index, args)
    else:                       # out-of-phase

        if len(c1.seq) != 3 or len(c2.seq) != 3:
            if len(t.seq) % 3 != 0: # truncated transcript sequence
                raise IncompatibleTranscriptError("truncated_transcript_sequence_at_boundary_(start_codon_seq_%s_and_end_codon_seq_%s)" % (c1.seq, c2.seq))
            else:               # unknown error
                raise IncompatibleTranscriptError("unknown_error_causing_codon_seq_not_multiplicative_of_3")

        beg_codon_beg = c1.index*3-2
        end_codon_end = c2.index*3
        new_codon_seq = t.seq[beg_codon_beg-1:p1.pos-1] + \
                        t.seq[p2.pos:end_codon_end]

        if len(new_codon_seq) != 3:
            # err_print(t.gene_name+'\t'+t.transcript_type)
            # err_print(p1)
            # err_print(p2)
            # err_print(len(t.seq))
            # err_print(c1.seq)
            # err_print(c2.seq)
            # err_print(len(t.seq) % 3)
            # err_print(t.seq[-10:])
            # if (len(t.seq)%3 != 0):
            #     r.append_info('truncated_refseq_at_boundary_(codon_seq_%s)' % c1.seq)
            raise IncompatibleTranscriptError('unknown_error_causing_codon_seq_not_multiplicative_of_3_%s' % new_codon_seq)

        taa_alt = codon2aa(new_codon_seq)
        tnuc_delseq = t.seq[beg_codon_beg-1:end_codon_end]
        taa_delseq = translate_seq(tnuc_delseq)
        # if taa_delseq[-1] == '*':
        if taa_alt == taa_delseq[-1]:
            # G100_S200delinsS becomes a pure deletion G100_D199del
            r.csqn.append("InFrameDeletion")
            taa_set_del(r, t, c1.index, c2.index-1, args)
        elif taa_alt == taa_delseq[0]:
            # S100_G200delinsS becomes a pure deletion D101_G200del
            r.csqn.append('InFrameDeletion')
            taa_set_del(r, t, c1.index+1, c2.index, args)
        elif new_codon_seq == '*':
            r.csqn.append("Nonsense")
            r.taa_range = "%s%d%s" % (aaf(taa_delseq[0], args), c1.index, taa_alt)
            variant_protein_seq_sub(r, t, args, c1.index, c1.index, taa_alt)
        else:
            r.csqn.append('MultiAAMissense')
            r.taa_range = '%s%d_%s%ddelins%s' % (
                aaf(taa_delseq[0], args), c1.index,
                aaf(taa_delseq[-1], args), c2.index, taa_alt)
            variant_protein_seq_sub(r, t, args, c1.index, c2.index, taa_alt)

def del_coding_frameshift(args, cbeg, cend, pbeg, pend, t, r):

    """ assume frame-shift does not affect splicing """
    cbeg_beg = cbeg.index*3 - 2
    old_seq = t.seq[cbeg_beg-1:]
    new_seq = t.seq[cbeg_beg-1:pbeg.pos-1]+t.seq[pend.pos:]
    if not old_seq:
        raise IncompatibleTranscriptError("invalid_cDNA_position_%d;expect_[0_%d]" % cbeg_beg, len(t.seq))

    aae = t.extend_taa_seq(cbeg.index, old_seq, new_seq)
    if aae:
        r.taa_range = aae.format(args)
        r.csqn.append("Frameshift")
        variant_protein_seq_fs(r, t, aae, args)
    else: # rare chance when stop codon seen before difference
        r.taa_range = '(=)'
        r.csqn.append("Synonymous")
