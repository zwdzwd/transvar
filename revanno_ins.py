""" in-frame and frameshift must be separated for coding sequence """
from transcripts import *
from utils import *
from record import *

def nuc_mutation_ins_coding_inframe_inphase(args, q, tpt, r):

    """ insertion is after the 3rd base of a codon """

    # check stop codon
    taa_insseq = ''
    for i in xrange(len(q.insseq)/3):
        if codon2aa(q.insseq[i*3:i*3+3]) == '*':
            return nuc_mutation_ins_coding_frameshift(args, q, tpt, r)
        taa_insseq += codon2aa(q.insseq[i*3:i*3+3])

    # otherwise, a pure insertion
    c1 = tpt.cpos2codon((q.pos.pos+2)/3)
    c2 = tpt.cpos2codon((q.pos.pos+3)/3)
    if not c1 or not c2: raise IncompatibleTranscriptError()
    r.taa_range = '%s%d_%s%dins%s' % (codon2aa(c1.seq), c1.index, 
                                      codon2aa(c2.seq), c2.index, 
                                      taa_insseq)
    refinsseq = q.insseq if tpt.strand == '+' else reverse_complement(q.insseq)
    gnuc_beg, gnuc_end = tpt.tnuc_range2gnuc_range(q.pos.pos, q.pos.pos+1)
    r.gnuc_range = '%d_%dins%s' % (gnuc_beg, gnuc_end, refinsseq)
    r.tnuc_range = '%d_%dins%s' % (q.pos.pos, q.pos.pos+1, q.insseq)
    r.pos = '%d-(ins)-%d' % (gnuc_beg, gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
    r.info = 'NatInsSeq=%s;RefInsSeq=%s;Phase=0' % (q.insseq, refinsseq)

def nuc_mutation_ins_coding_inframe_outphase(args, q, tpt, r):

    """ insertion is after 1st or 2nd base of a codon """

    codon_index = (q.pos.pos+2)/3
    codon = tpt.cpos2codon(codon_index)
    if not codon: raise IncompatibleTranscriptError()

    codon_beg = codon_index*3-2
    codon_end = codon_index*3
    # 0, 1,2 indicating insertion happen after 3rd, 1st or 2nd base of the codon
    phase = q.pos.pos - codon_beg + 1
    codon_subseq1 = tpt.seq[codon_beg-1:q.pos.pos]
    codon_subseq2 = tpt.seq[q.pos.pos:codon_end]
    new_seq = codon_subseq1+q.insseq+codon_subseq2
    taa_insseq = ''
    for i in xrange(len(new_seq)/3):
        if codon2aa(new_seq[i*3:i*3+3]) == '*':
            return nuc_mutation_ins_coding_frameshift(args, q, tpt, r)
        taa_insseq += codon2aa(new_seq[i*3:i*3+3])

    if not codon: raise IncompatibleTranscriptError()
    taa_ref = codon2aa(codon.seq)
    #print codon_beg, len(tpt.seq), new_seq, codon_subseq1, codon_subseq2
    if taa_ref == taa_insseq[0]:
        # SdelinsSH becomes a pure insertion [current_codon]_[codon_after]insH
        taa_ref_after = codon2aa(tpt.seq[codon.index*3:codon.index*3+3])
        r.taa_range = '%s%d_%s%dins%s' % (taa_ref, codon.index,
                                          taa_ref_after, codon.index+1, taa_insseq[1:])
    elif taa_ref == taa_insseq[-1]:
        # SdelinsHS becomes a pure insertion [codon_before]_[current_codon]insH
        taa_ref_before = codon2aa(tpt.seq[codon.index*3-6:codon.index*3-3])
        r.taa_range = '%s%d_%s%dins%s' % (taa_ref_before, codon.index-1,
                                          taa_ref, codon.index, taa_insseq[:-1])
    else:
        r.taa_range = '%s%ddelins%s' % (taa_ref, codon.index, taa_insseq)
    refinsseq = q.insseq if tpt.strand == '+' else reverse_complement(q.insseq)
    gnuc_beg, gnuc_end = tpt.tnuc_range2gnuc_range(q.pos.pos, q.pos.pos+1)
    r.gnuc_range = '%d_%dins%s' % (gnuc_beg, gnuc_end, refinsseq)
    r.tnuc_range = '%d_%dins%s' % (q.pos.pos, q.pos.pos+1, q.insseq)
    r.pos = '%d-(ins)-%d' % (gnuc_beg, gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
    r.info = 'NatInsSeq=%s(%s)%s;RefInsSeq=%s;Phase=%d' % (codon_subseq1, q.insseq, codon_subseq2, refinsseq, phase)

def nuc_mutation_ins_coding_inframe(args, q, tpt, r):

    # if stop codon, represented as a deletion
    if q.pos.pos % 3 == 0:
        nuc_mutation_ins_coding_inframe_inphase(args, q, tpt, r)
    else:
        nuc_mutation_ins_coding_inframe_outphase(args, q, tpt, r)

def nuc_mutation_ins_coding_frameshift(args, q, tpt, r):

    beg_codon_index = (q.pos.pos + 2) / 3
    beg_codon_beg = beg_codon_index*3 - 2
    if beg_codon_beg+3 > len(tpt.seq): raise IncompatibleTranscriptError()
    #print beg_codon_beg, len(tpt.seq)
    #print tpt.name, tpt.chrm, tpt.exons
    #print tpt.seq
    old_seq = tpt.seq[beg_codon_beg-1:]
    new_seq = tpt.seq[beg_codon_beg-1:q.pos.pos]+q.insseq+tpt.seq[q.pos.pos:]
    ret = extend_taa_seq(beg_codon_index, old_seq, new_seq, tpt)
    if ret:
        taa_pos, taa_ref, taa_alt, termlen = ret
        r.taa_range = '%s%d%sfs*%s' % (taa_ref, taa_pos, taa_alt, termlen)
    else:
        r.taa_range = '(=)'
    r.tnuc_range = '%d_%dins%s' % (q.pos.pos, q.pos.pos+1, q.insseq)
    gnuc_beg, gnuc_end = tpt.tnuc_range2gnuc_range(q.pos.pos, q.pos.pos+1)
    refinsseq = q.insseq if tpt.strand == '+' else reverse_complement(q.insseq)
    r.gnuc_range = '%d_%dins%s' % (gnuc_beg, gnuc_end, refinsseq)
    r.info = 'NatInsSeq=%s;RefInsSeq=%s' % (q.insseq, refinsseq)
    r.pos = '%d-%d' % (gnuc_beg, gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)

def nuc_mutation_ins_coding(args, q, tpt, r):

    """ assuming insertion does not affect splicing """
    if len(q.insseq) % 3 == 0:
        nuc_mutation_ins_coding_inframe(args, q, tpt, r)
    else:
        nuc_mutation_ins_coding_frameshift(args, q, tpt, r)

    return

def nuc_mutation_ins_intronic(args, q, tpt, r):

    """ assuming insertion does not affect splicing """

    codon = tpt.cpos2codon((q.pos.pos+2)/3)
    if not codon: raise IncompatibleTranscriptError("No codon")
    #print q.pos.pos, codon.index
    #print codon.locs
    i = q.pos.pos - (codon.index-1)*3 - 1
    if tpt.strand == '-': i = 2-i
    #print i
    if tpt.strand == '+':
        gnuc_beg = codon.locs[i] + q.pos.tpos
        gnuc_end = codon.locs[i] + q.pos.tpos + 1
        if q.pos.tpos > 0: j = i+1
        elif q.pos.tpos < 0: j = i
        refinsseq = q.insseq
    else:
        gnuc_beg = codon.locs[i] - q.pos.tpos - 1
        gnuc_end = codon.locs[i] - q.pos.tpos
        if q.pos.tpos > 0: j = i
        elif q.pos.tpos < 0: j = i+1
        refinsseq = reverse_complement(q.insseq)
    #print gnuc_beg, gnuc_end
    #print j
    if j>0 and j<3 and codon.locs[j-1]+1 == codon.locs[j]:
        raise IncompatibleTranscriptError('Codon does not contain exon boundary')

    pl = []
    s = '-'.join(map(str, codon.locs[:j]))
    if s: pl.append(s)
    if not (j>0 and codon.locs[j-1] == gnuc_beg): pl.append('(%d)' % gnuc_beg)
    pl.append('(ins)')
    if not (j<3 and codon.locs[j] == gnuc_end): pl.append('(%d)' % gnuc_end)
    s = '-'.join(map(str, codon.locs[j:]))
    if s: pl.append(s)
    r.pos = '-'.join(pl)
    r.reg = '%s (%s, intronic)' % (tpt.gene.name, tpt.strand)
    r.gnuc_range = '%d_%dins%s' % (gnuc_beg, gnuc_end, refinsseq)
    r.tnuc_range = '%d%d_%d%dins%s' % (q.pos.pos, q.pos.tpos, q.pos.pos, q.pos.tpos+1, q.insseq)
    r.info = 'RefInsSeq=%s;NatInsSeq=%s' % (refinsseq, q.insseq)

def nuc_mutation_ins(args, q, tpt):
    
    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    tpt.ensure_seq()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name
    r.muttype = 'ins'

    if q.pos.tpos == 0:
        nuc_mutation_ins_coding(args, q, tpt, r)
    else:
        nuc_mutation_ins_intronic(args, q, tpt, r)

    return r

def _core_annotate_nuc_ins(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = nuc_mutation_ins(args, q, tpt)
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
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

    return


def codon_mutation_ins(args, q, tpt):
    
    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError("Transcript name unmatched")
    tpt.ensure_seq()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name

    if q.beg*3 > len(tpt) or q.end*3 > len(tpt):
        raise IncompatibleTranscriptError('codon nonexistent')

    tnuc_beg = q.beg*3-2
    tnuc_end = q.end*3
    if hasattr(q, 'beg_aa') and q.beg_aa and q.beg_aa != tpt.taa2aa(q.beg):
        raise IncompatibleTranscriptError('Unmatched reference amino acid')
    if hasattr(q, 'end_aa') and q.end_aa and q.end_aa != tpt.taa2aa(q.end):
        raise IncompatibleTranscriptError('Unmatched reference amino acid')
    gnuc_beg, gnuc_end = tpt.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
    r.tnuc_range = '(%d-%d)ins%d' % (tnuc_beg-1, tnuc_end, len(q.insseq)*3)
    r.gnuc_range = '(%d-%d)ins%d' % (gnuc_beg-1, gnuc_end, len(q.insseq)*3)
    r.pos = '%d-%d (insertion)' % (gnuc_beg, gnuc_end)

    return r

def _core_annotate_codon_ins(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = codon_mutation_ins(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue
        except UnknownChromosomeError:
            continue
        r.muttype = 'ins'
        r.taa_range = '%s%s_%s%sins%s' % (q.beg_aa, str(q.beg), q.end_aa, str(q.end), q.insseq)
        r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
        r.info = 'Uncertain' if r.info == '.' else (r.info+';Uncertain')
        r.format(q.op)
        found = True

    if not found:
        r = Record()
        r.taa_range = '%s%s_%s%sins%s' % (q.beg_aa, str(q.beg), q.end_aa, str(q.end), q.insseq)
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)


