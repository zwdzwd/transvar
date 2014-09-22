
from transcripts import *
from utils import *
from record import *

def nuc_mutation_ins_coding_inframe_inphase(args, q, tpt, r):

    """ insertion is after the 3rd base of a codon """

    # check stop codon
    taa_insseq = ''
    for i in xrange(len(q.insseq)/3):
        if standard_codon_table[q.insseq[i*3:i*3+3]] == '*':
            nuc_mutation_ins_coding_outframe(args, q, tpt, r)
            return
        taa_insseq += standard_codon_table[q.insseq[i*3:i*3+3]]

    # otherwise, a pure insertion
    c1 = tpt.cpos2codon((q.pos.pos+2)/3)
    c2 = tpt.cpos2codon((q.pos.pos+3)/3)
    if not c1 or not c2: raise IncompatibleTranscriptError()
    r.taa_range = '%s%d_%s%dins%s' % (standard_codon_table[c1.seq], c1.index, 
                                      standard_codon_table[c2.seq], c2.index, 
                                      taa_insseq)

def nuc_mutation_ins_coding_inframe_outphase(args, q, tpt, r):

    """ insertion is after 1st or 2nd base of a codon """

    codon = tpt.cpos2codon((q.pos.pos+2)/3)
    if not codon: raise IncompatibleTranscriptError()

    i = q.pos.pos % 3
    c1.seq 

def nuc_mutation_ins_coding_inframe(args, q, tpt, r):

    # if stop codon, represented as a deletion
    if q.pos.pos % 3 == 0:
        nuc_mutation_ins_coding_inframe_inphase(args, q, tpt, r)
    else:
        nuc_mutation_ins_coding_inframe_outphase(args, q, tpt, r)

def nuc_mutation_ins_coding_outframe(args, q, tpt, r):

    beg_codon_index = (q.pos.pos + 2) / 3
    beg_codon_beg = beg_codon_index*3 - 2
    old_seq = tpt.seq[beg_codon_beg-1:]
    new_seq = tpt.seq[beg_codon_beg-1:q.pos.pos]+q.insseq+tpt.seq[q.pos.pos:]
    taa_pos, taa_ref, taa_alt, termlen = extend_taa_seq(beg_codon_index, old_seq, new_seq)
    r.taa_range = '%s%d%sfs*%d' % (taa_ref, taa_pos, taa_alt, termlen)
    r.tnuc_range = '%d_%dins%s' % (q.pos.pos, q.pos.pos+1, q.insseq)
    gnuc_beg, gnuc_end = tpt.tnuc_range2gnuc_range(q.pos.pos, q.pos.pos+1)
    refinsseq = q.insseq if tpt.strand == '+' else reverse_complement(q.insseq)
    r.gnuc_range = '%d_%dins%s' % (gnuc_beg, gnuc_end, refinsseq)
    r.info = 'NatInsSeq=%s;RefInsSeq=%s' % (q.insseq, refinsseq)
    r.pos = '%s:%d-%d' % (r.chrm, gnuc_beg, gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)

    return


def nuc_mutation_ins_coding(args, q, tpt, r):

    """ assuming insertion does not affect splicing """
    if len(q.insseq) % 3 == 0:
        nuc_mutation_ins_coding_inframe(args, q, tpt, r)
    else:
        nuc_mutation_ins_coding_outframe(args, q, tpt, r)

    return

def nuc_mutation_ins(args, q, tpt):
    
    if q.tpt and tpt.name != q.tpt: raise IncompatibleTranscriptError("Transcript name unmatched")
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
        found = True
        r.format(q.op)

    if not found:
        r = Record()
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

    return
