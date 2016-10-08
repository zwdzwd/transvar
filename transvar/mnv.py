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
from __future__ import division
from .err import *
from .record import *
from .transcripts import *
from .describe import *
from .insertion import taa_set_ins, annotate_insertion_gdna
from .deletion import taa_set_del, annotate_deletion_gdna
from .snv import annotate_snv_gdna
from .proteinseqs import *

def annotate_mnv_cdna(args, q, tpts, db):

    records = []
    for t in tpts:
        try:

            if q.tpt and t.name != q.tpt:
                raise IncompatibleTranscriptError("unmatched_transcript_name_%s;expect_%s" % (q.tpt, t.name))
            t.ensure_seq()

            r = Record(is_var=True)
            r.chrm = t.chrm
            r.tname = t.format()
            r.gene = t.gene_name
            r.strand = t.strand

            t.check_exon_boundary(q.beg)
            t.check_exon_boundary(q.end)

            _gnuc_beg = t.tnuc2gnuc(q.beg)
            _gnuc_end = t.tnuc2gnuc(q.end)
            gnuc_beg = min(_gnuc_beg, _gnuc_end)
            gnuc_end = max(_gnuc_beg, _gnuc_end)
            r.pos = '%d-%d' % (gnuc_beg, gnuc_end)

            gnuc_refseq = faidx.getseq(t.chrm, gnuc_beg, gnuc_end)
            tnuc_refseq = reverse_complement(gnuc_refseq) if t.strand == '-' else gnuc_refseq
            gnuc_altseq = reverse_complement(q.altseq) if t.strand == '-' else q.altseq
            if q.refseq and tnuc_refseq != q.refseq:
                raise IncompatibleTranscriptError('reference_unmatched_%s_expect_%s' % (q.refseq, tnuc_refseq))

            r.gnuc_range = nuc_set_mnv(gnuc_beg, gnuc_end, gnuc_refseq, gnuc_altseq)
            r.tnuc_range = nuc_set_mnv(q.beg, q.end, tnuc_refseq, q.altseq)

            r.reg = describe_genic(args, t.chrm, gnuc_beg, gnuc_end, t, db)
            if (not r.set_splice('lost', 'BlockSubstitution')):
                if t.transcript_type == 'protein_coding' and r.reg.entirely_in_cds():
                    try:
                        tnuc_mnv_coding(t, q.beg.pos, q.end.pos, q.altseq, r, args)
                    except IncompatibleTranscriptError as e:
                        r.append_info(e.message)
                        # _beg, _end, _seqlen = inst
                        # r.append_info('mnv_(%s-%s)_at_truncated_refseq_of_length_%d' % (_beg, _end, _seqlen))
                else:
                    r.csqn.append(r.reg.csqn()+"BlockSubstitution")

        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue

        records.append(r)
    format_records(records, q.op, args)
    return records

def annotate_mnv_protein(args, q, tpts, db):

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

            tnuc_beg = q.beg*3-2
            tnuc_end = q.end*3
            gnuc_beg, gnuc_end = t.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
            tnuc_refseq = t.seq[tnuc_beg-1:tnuc_end]
            gnuc_refseq = reverse_complement(tnuc_refseq) if t.strand == '-' else tnuc_refseq
            taa_refseq = translate_seq(tnuc_refseq)
            if q.beg_aa and q.beg_aa != taa_refseq[0]:
                raise IncompatibleTranscriptError('beginning reference amino acid unmatched')
            if q.end_aa and q.end_aa != taa_refseq[-1]:
                raise IncompatibleTranscriptError('ending reference amino acid unmatched')
            if q.refseq and taa_refseq != q.refseq:
                raise IncompatibleTranscriptError('reference sequence unmatched')
            # reverse translate
            tnuc_altseq = []
            cdd_altseq = []
            for aa in q.altseq:
                tnuc_altseq.append(aa2codon(aa)[0])
                cdd_altseq.append('/'.join(aa2codon(aa)))
            tnuc_altseq = ''.join(tnuc_altseq)
            gnuc_altseq = reverse_complement(tnuc_altseq) if t.strand == '-' else tnuc_altseq
            r.tnuc_range = nuc_set_mnv(tnuc_beg, tnuc_end, tnuc_refseq, tnuc_altseq)
            r.gnuc_range = nuc_set_mnv(gnuc_beg, gnuc_end, gnuc_refseq, gnuc_altseq)
            r.pos = '%d-%d' % (gnuc_beg, gnuc_end)
            r.csqn.append("MultiAAMissense")
            if len(cdd_altseq) <= 2:
                r.append_info('candidate_alternative_sequence=%s' % ('+'.join(cdd_altseq), ))
            else:
                r.append_info('%d_CandidatesOmitted' % (aaseq_redundancy(q.altseq)))

        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue
        r.taa_range = '%s%s_%s%sdelins%s' % (
            aaf(q.beg_aa, args), str(q.beg), aaf(q.end_aa, args), str(q.end), aaf(q.altseq, args)) # q.refseq,
        r.reg = RegCDSAnno(t)
        r.reg.from_taa_range(q.beg, q.end)

        records.append(r)

    format_records(records, q.op, args)
    return records

    # if not found:
    #     r = Record(is_var=True)
    #     r.taa_range = '%s%s_%s%sdelins%s' % (
    #         aaf(q.beg_aa, args), str(q.beg), aaf(q.end_aa, args), str(q.end), aaf(q.altseq, args)) # q.refseq,
    #     r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))

    #     r.format(q.op)

def decompose_mut(q):
    from . import ssw
    # print q.altseq
    # print q.refseq
    aln = ssw.ssw_aln(q.altseq, q.refseq)
    # b = 'GGGGGGGGGCGTACCCTGGAG'
    # a = 'GCTACCCAGGAG'
    # aln = ssw.ssw_aln(a, b)

    # if aln.rbeg == 0 and aln.qbeg == 0:
    rpos = 0
    qpos = 0
    for ct, cl in aln.cigar:
        if ct == 0:
            _altseq = ''
            _refseq = ''
            _coords = []
            for i in range(cl):
                if q.altseq[qpos+i] == q.refseq[rpos+i]:
                    _altseq += ' '
                    _refseq += ' '
                else:
                    if i == 0 or q.altseq[qpos+i-1] == q.refseq[rpos+i-1]:
                        _coords.append(i)
                    _altseq += q.altseq[qpos+i]
                    _refseq += q.refseq[rpos+i]
            ss = list(zip(_coords, _altseq.strip().split(), _refseq.strip().split()))
            for j, alt, ref in ss:
                if len(ref) == 1:
                    qq = QuerySNV()
                    qq.pos = q.beg+rpos+j
                    qq.ref = ref
                    qq.alt = alt
                    yield qq
                else:
                    qq = QueryMNV()
                    qq.beg = q.beg+rpos+j
                    qq.end = q.beg+rpos+j+len(ref)-1
                    qq.refseq = ref
                    qq.altseq = alt
                    yield qq
            rpos += cl
            qpos += cl
        elif ct == 1:
            qq = QueryINS()
            qq.pos = q.beg+rpos-1
            qq.insseq = q.altseq[qpos:qpos+cl]
            yield qq
            qpos += cl
        elif ct == 2:
            qq = QueryDEL()
            qq.beg = q.beg+rpos
            qq.end = q.beg+rpos+cl-1
            qq.delseq = q.refseq[rpos:rpos+cl]
            yield qq
            rpos += cl
        elif ct == 4:
            qq = QueryMNV()
            qq.beg = q.beg+rpos
            qq.end = q.beg+rpos+cl-1
            qq.refseq = q.refseq[rpos:rpos+cl]
            qq.altseq = q.altseq[qpos:qpos+cl]
            yield qq
            rpos += cl
            qpos += cl

def annotate_mnv_gdna(args, q, db):

    # check reference sequence
    gnuc_refseq = faidx.refgenome.fetch_sequence(q.tok, q.beg, q.end)
    if q.refseq and gnuc_refseq != q.refseq:

        r = Record(is_var=True)
        r.chrm = q.tok
        r.pos = '%d-%d' % (q.beg, q.end)
        r.info = "invalid_reference_seq_%s_(expect_%s)" % (q.refseq, gnuc_refseq)
        r.format(q.op)
        err_print("Warning: %s invalid reference %s (expect %s), maybe wrong reference?" % (q.op, q.refseq, gnuc_refseq))
        return

    else:                       # make sure q.refseq exists
        q.refseq = gnuc_refseq

    if args.haplotype:
        from . import anno
        for qq in decompose_mut(q):
            qq.op = q.op
            qq.tok = q.tok
            anno._main_core_(args, qq, db, 'g')
        return

    gnuc_altseq = q.altseq
    gnuc_refseq, gnuc_altseq, head_trim, tail_trim = double_trim(gnuc_refseq, gnuc_altseq)
    q.beg += head_trim
    q.end -= tail_trim

    if q.beg == q.end and len(gnuc_altseq) == 1:
        q.pos = q.beg
        q.ref = gnuc_refseq
        q.alt = gnuc_altseq
        annotate_snv_gdna(args, q, db)
        return

    if len(gnuc_refseq) == 0:
        if len(gnuc_altseq) > 0:
            q.pos = q.beg-1
            q.insseq = gnuc_altseq
            annotate_insertion_gdna(args, q, db)
            return

    if len(gnuc_altseq) == 0:
        if len(gnuc_refseq) > 0:
            q.delseq = gnuc_refseq
            annotate_deletion_gdna(args, q, db)
            return

    records = []
    for reg in describe(args, q, db):

        r = Record(is_var=True)
        r.reg = reg
        r.chrm = q.tok
        r.pos = '%d-%d' % (q.beg, q.end)
        r.gnuc_range = nuc_set_mnv(q.beg, q.end, gnuc_refseq, gnuc_altseq)

        db.query_dbsnp_range(r, q.beg, q.end, gnuc_altseq)
        if hasattr(reg, 't'):

            t = reg.t
            r.tname = t.format()
            r.gene = t.gene_name
            r.strand = t.strand

            c1, p1 = t.gpos2codon(q.beg)
            c2, p2 = t.gpos2codon(q.end)
            if t.strand == '+':
                tnuc_beg = p1
                tnuc_end = p2
                tnuc_refseq = gnuc_refseq
                tnuc_altseq = gnuc_altseq
            else:
                tnuc_beg = p2
                tnuc_end = p1
                tnuc_refseq = reverse_complement(gnuc_refseq)
                tnuc_altseq = reverse_complement(gnuc_altseq)
            r.tnuc_range = nuc_set_mnv(tnuc_beg, tnuc_end, tnuc_refseq, tnuc_altseq)

            if not r.set_splice('lost', 'BlockSubstitution'):
                if r.reg.t.transcript_type == 'protein_coding' and r.reg.entirely_in_cds():
                    try:
                        _, tnuc_beg_adj = t.intronic_lean(tnuc_beg, 'c_greater')
                        _, tnuc_end_adj = t.intronic_lean(tnuc_end, 'c_smaller')
                        tnuc_mnv_coding(t, tnuc_beg_adj.pos, tnuc_end_adj.pos, tnuc_altseq, r, args)
                    except IncompatibleTranscriptError as e:
                        r.append_info(e.message)
                        # wrap_exception(e, q, args)
                        # if len(inst) == 3:
                        #     _beg, _end, _seqlen = inst
                        #     r.append_info('mnv_(%s-%s)_at_truncated_refseq_of_length_%d' % (_beg, _end, _seqlen))
                        # else:
                        #     raise inst
                else:
                    r.csqn.append(r.reg.csqn()+"BlockSubstitution")

        elif isinstance(reg, RegSpanAnno):

            r.csqn.append(reg.csqn()+"BlockSubstitution")
            tnames = []
            strands = []
            genes = []
            if hasattr(reg.b1, 't'):
                if reg.b1.t.name not in tnames:
                    tnames.append(reg.b1.t.name)
                    strands.append(reg.b1.t.strand)
                    genes.append(reg.b1.t.gene_name)

            if hasattr(reg.b2, 't'):
                if reg.b2.t.name not in tnames:
                    tnames.append(reg.b2.t.name)
                    strands.append(reg.b2.t.strand)
                    genes.append(reg.b2.t.gene_name)

            if tnames:
                r.tname = ','.join(tnames)
            if strands:
                r.strand = ','.join(strands)
            if genes:
                r.gene = ','.join(genes)

        records.append(r)

    format_records(records, q.op, args)
    return records

def tnuc_mnv_coding(t, beg, end, altseq, r, args):

    if (len(altseq) - (end-beg+1)) % 3 == 0: # in frame

        # beg and end are integer tnuc positions
        # altseq follows the tnuc (cDNA) order
        # set taa range

        beg_codon_index = (beg + 2) // 3
        end_codon_index = (end + 2) // 3

        beg_codon_beg = beg_codon_index*3 - 2
        end_codon_end = end_codon_index*3 # 1 past the last codon

        old_seq = t.seq[beg_codon_beg-1:end_codon_end]
        new_seq = t.seq[beg_codon_beg-1:beg-1]+altseq+t.seq[end:end_codon_end]

        if beg_codon_index == end_codon_index:
            r.append_info('codon_cDNA=%s' %
                          '-'.join(map(str, list(range(beg_codon_beg, beg_codon_beg+3)))))
        else:
            r.append_info('begin_codon_cDNA=%s' %
                          '-'.join(map(str, list(range(beg_codon_beg, beg_codon_beg+3)))))
            r.append_info('end_codon_cDNA=%s' %
                          '-'.join(map(str, list(range(end_codon_end-2, end_codon_end+1)))))

        if len(old_seq) % 3 != 0:
            # raise IncompatibleTranscriptError(beg, end, len(t.seq))
            raise IncompatibleTranscriptError(
                'mnv_[%s_%s]_at_truncated_refseq_of_length_%d' % (beg, end, len(t.seq)))

        old_taa_seq = translate_seq(old_seq)
        new_taa_seq = translate_seq(new_seq)
        if old_taa_seq == new_taa_seq:
            r.csqn.append("Synonymous")
            r.taa_range = '(=)'
            return

        # block substitution in nucleotide level may end up
        # an insertion or deletion on the protein level
        old_taa_seq1, new_taa_seq1, head_trim, tail_trim = double_trim(old_taa_seq, new_taa_seq)
        if not old_taa_seq1:
            _beg_index = beg_codon_index + head_trim - 1
            _end_index = beg_codon_index + head_trim
            r.csqn.append("InFrameInsertion")
            taa_set_ins(r, t, _beg_index, new_taa_seq1, args)
            return

        if not new_taa_seq1:
            r.csqn.append("InFrameDeletion")
            taa_set_del(r, t, beg_codon_index+head_trim,
                        end_codon_index-tail_trim, args)
            return

        if len(old_taa_seq1) == 1:
            if len(new_taa_seq1) == 1:
                r.csqn.append("Missense")
                taa_pos = beg_codon_index + head_trim
                taa_alt = aaf(new_taa_seq1, args)
                r.taa_range = '%s%d%s' % (
                    aaf(old_taa_seq1[0], args), taa_pos, taa_alt)
            else:
                r.csqn.append("MultiAAMissense")
                taa_pos = beg_codon_index + head_trim
                taa_alt = aaf(new_taa_seq1, args)
                r.taa_range = '%s%ddelins%s' % (
                    aaf(old_taa_seq1[0], args), taa_pos, taa_alt)
            variant_protein_seq_sub(r, t, args, taa_pos, taa_pos, taa_alt)
            return

        r.csqn.append("MultiAAMissense")
        taa_beg = beg_codon_index + head_trim
        taa_end = end_codon_index - tail_trim
        taa_alt = aaf(new_taa_seq1, args)
        r.taa_range = '%s%d_%s%ddelins%s' % (
            aaf(old_taa_seq1[0], args), taa_beg,
            aaf(old_taa_seq1[-1], args), taa_end, taa_alt)
        variant_protein_seq_sub(r, t, args, taa_beg, taa_end, taa_alt)

    else:                   # frame-shift

        beg_codon_index = (beg + 2) // 3
        beg_codon_beg = beg_codon_index * 3 - 2
        old_seq = t.seq[beg_codon_beg-1:]
        new_seq = t.seq[beg_codon_beg-1:beg-1]+altseq+t.seq[end:]

        aae = t.extend_taa_seq(beg_codon_index, old_seq, new_seq)
        if aae:
            r.csqn.append("Frameshift")
            r.taa_range = aae.format(args)
            variant_protein_seq_fs(r, t, aae, args)
        else:
            r.csqn.append("Synonymous")
            r.taa_range = '(=)'

def nuc_set_mnv(beg, end, refseq, altseq):

    if beg == end:
        if len(altseq) == 1:
            return '%s%s>%s' % (beg, refseq, altseq)
        else:
            return '%sdelins%s' % (beg, altseq)
    else:
        return '%s_%sdelins%s' % (beg, end, altseq)
