""" annotate block substitution,
assume this does not affect splice site """
from err import *
from record import *
from transcripts import *

def nuc_mutation_mnv_coding_inframe(args, q, tpt, r):

    beg_codon_index = (q.beg.pos + 2) / 3
    end_codon_index = (q.end.pos + 2) / 3
    
    beg_codon_beg = beg_codon_index*3 - 2
    end_codon_end = end_codon_index*3 # 1 past the last codon

    old_seq = tpt.seq[beg_codon_beg-1:end_codon_end]
    new_seq = tpt.seq[beg_codon_beg-1:q.beg.pos-1]+q.altseq+tpt.seq[q.end.pos:end_codon_end]
    old_taa_seq = translate_seq(old_seq)
    new_taa_seq = translate_seq(new_seq)
    if old_taa_seq == new_taa_seq:
        r.taa_range = '(=)'
    else:
        # block substitution in nucleotide level may end up
        # an insertion or deletion on the protein level
        old_taa_seq1, new_taa_seq1, head_trim, tail_trim = double_trim(old_taa_seq, new_taa_seq)
        if not old_taa_seq1:
            r.taa_range = '%s%d_%s%dins%s' % (
                beg_codon_index + head_trim, old_taa_seq1[0],
                beg_codon_index + head_trim + 1, old_taa_seq1[1], 
                new_taa_seq1)
        elif not new_taa_seq1:
            if len(old_taa_seq1) == 1:
                r.taa_range = '%s%ddel' % (old_taa_seq1[0], beg_codon_index + head_trim)
            else:
                r.taa_range = '%s%d_%s%ddel' % (
                    old_taa_seq1[0], beg_codon_index + head_trim, 
                    old_taa_seq1[1], end_codon_index - tail_trim)
        else:
            if len(old_taa_seq1) == 1:
                if len(new_taa_seq1) == 1:
                    r.taa_range = '%s%d%s' % (
                        old_taa_seq1[0], beg_codon_index + head_trim, new_taa_seq1)
                else:
                    r.taa_range = '%s%ddelins%s' % (
                        old_taa_seq1[0], beg_codon_index + head_trim, new_taa_seq1)
            else:
                r.taa_range = '%s%d_%s%ddelins%s' % (
                    old_taa_seq1[0], beg_codon_index + head_trim,
                    old_taa_seq1[-1], end_codon_index - tail_trim, new_taa_seq1)

    return

def nuc_mutation_mnv_coding_frameshift(orgs, q, tpt, r):

    beg_codon_index = (q.beg.pos + 2) / 3
    beg_codon_beg = beg_codon_index * 3 - 2
    old_seq = tpt.seq[beg_codon_beg-1:]
    new_seq = tpt.seq[beg_codon_beg-1:q.beg.pos-1]+q.altseq+tpt.seq[q.end.pos:]
    print old_seq[:10], new_seq[:10]
    ret = extend_taa_seq(beg_codon_beg, old_seq, new_seq, tpt)
    if ret:
        taa_pos, taa_ref, taa_alt, termlen = ret
        r.taa_range = '%s%d%sfs*%s' % (taa_ref, taa_pos, taa_alt, termlen)
    else:
        r.taa_range = '(=)'

    return

def nuc_mutation_mnv_coding(args, q, tpt, r):

    if (len(q.altseq) - (q.end.pos-q.beg.pos+1))  % 3 == 0:
        nuc_mutation_mnv_coding_inframe(args, q, tpt, r)
    else:
        nuc_mutation_mnv_coding_frameshift(args, q, tpt, r)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)

def nuc_mutation_mnv(args, q, tpt):

    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    tpt.ensure_seq()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name
    r.muttype = 'mnv'

    np = tpt.position_array()
    check_exon_boundary(np, q.beg)
    check_exon_boundary(np, q.end)

    gnuc_beg = tnuc2gnuc2(np, q.beg, tpt)
    gnuc_end = tnuc2gnuc2(np, q.end, tpt)
    r.gnuc_beg = min(gnuc_beg, gnuc_end)
    r.gnuc_end = max(gnuc_beg, gnuc_end)
    tnuc_coding_beg = q.beg.pos if q.beg.tpos <= 0 else q.beg.pos+1 # base after intron
    tnuc_coding_end = q.end.pos if q.end.tpos >= 0 else q.end.pos-1 # base before intron
    gnuc_coding_beg = tnuc2gnuc(np, tnuc_coding_beg)
    gnuc_coding_end = tnuc2gnuc(np, tnuc_coding_end)

    r.refrefseq = faidx.refgenome.fetch_sequence(tpt.chrm, r.gnuc_beg, r.gnuc_end)
    r.natrefseq = reverse_complement(r.refrefseq) if tpt.strand == '-' else r.refrefseq
    r.refaltseq = reverse_complement(q.altseq) if tpt.strand == '-' else q.altseq
    if q.refseq and r.natrefseq != q.refseq:
        raise IncompatibleTranscriptError()

    r.gnuc_range = '%d_%d%s>%s' % (r.gnuc_beg, r.gnuc_end, r.refrefseq, r.refaltseq)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.tnuc_range = '%s_%s%s>%s' % (q.beg, q.end, r.natrefseq, q.altseq)

    reg = ''

    if region_in_exon(np, q.beg, q.end):
        nuc_mutation_mnv_coding(args, q, tpt, r)
    elif not region_in_intron(np, q.beg, q.end): # if block mutation occurs across splice site
        raise UnImplementedError('Block mutation occurs across splice site')

    return r

def _core_annotate_nuc_mnv(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = nuc_mutation_mnv(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
        found = True
        r.format(q.op)

    if not found:
        r = Record()
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

    return
