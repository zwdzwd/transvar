from transcripts import *
from utils import *
from record import *

def nuc_mutation_del_coding_inframe_inphase(args, q, tpt, r):
    """  deletion starts from the 1st base of the codon """

    beg_codon_index = (q.beg.pos + 2) / 3
    end_codon_index = (q.end.pos + 2) / 3
    if beg_codon_index == end_codon_index:
        beg_codon = tpt.cpos2codon(beg_codon_index)
        end_codon = beg_codon
        r.taa_range = '%s%ddel' % (standard_codon_table[beg_codon.seq], 
                                   beg_codon.index)
    else:
        beg_codon = tpt.cpos2codon(beg_codon_index)
        end_codon = tpt.cpos2codon(end_codon_index)
        r.taa_range = '%s%d_%s%ddel' % (standard_codon_table[beg_codon.seq],
                                        beg_codon.index,
                                        standard_codon_table[end_codon.seq],
                                        end_codon.index)
    natdelseq = tpt.seq[q.beg.pos-1:q.end.pos]
    r.info = 'NatDelSeq=%s' % natdelseq
    r.info += ';RefDelSeq=%s' % (natdelseq if tpt.strand == '+' else reverse_complement(natdelseq), )
    r.tnuc_range = '%d_%ddel' % (q.beg.pos, q.end.pos)
    if tpt.strand == '+':
        gnuc_beg = beg_codon.locs[0]
        gnuc_end = end_codon.locs[2]
    else:
        gnuc_beg = end_codon.locs[0]
        gnuc_end = beg_codon.locs[2]
    r.gnuc_range = '%d_%ddel' % (gnuc_beg, gnuc_end)
    r.pos = '%s:%d-%d' % (r.chrm, gnuc_beg, gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)

def nuc_mutation_del_coding_inframe_outphase(args, q, tpt, r):

    """ deletion starts from the 2nd/3rd base of the codon """

    r.muttype = 'delins'
    beg_codon_index = (q.beg.pos + 2) / 3
    end_codon_index = (q.end.pos + 2) / 3
    beg_codon_beg = beg_codon_index*3 - 2
    end_codon_end = end_codon_index*3
    # print q.beg.pos, q.end.pos
    # print beg_codon_index, end_codon_index
    # print beg_codon_beg, end_codon_end
    newcodonseq = tpt.seq[beg_codon_beg-1:q.beg.pos-1]+tpt.seq[q.end.pos:end_codon_end]
    r.taa_alt = standard_codon_table[newcodonseq]
    beg_codon_seq = tpt.seq[beg_codon_beg:beg_codon_beg+3]
    end_codon_seq = tpt.seq[end_codon_end-3:end_codon_end]
    r.taa_range = '%s%d_%s%d' % (standard_codon_table[beg_codon_seq], beg_codon_index, 
                                 standard_codon_table[end_codon_seq], end_codon_index)
    r.taa_range += 'delins%s' % r.taa_alt
    beg_codon = tpt.cpos2codon(beg_codon_index)
    end_codon = tpt.cpos2codon(end_codon_index)
    r.tnuc_range = '%d_%ddel' % (q.beg.pos, q.end.pos)
    natdelseq = tpt.seq[q.beg.pos-1:q.end.pos]
    r.info = 'NatDelSeq=%s' % natdelseq
    r.info += ';RefDelSeq=%s' % (natdelseq if tpt.strand == '+' else reverse_complement(natdelseq), )
    gnuc_del_beg = reverse_tnuc_pos(beg_codon, q.beg.pos)
    gnuc_del_end = reverse_tnuc_pos(end_codon, q.end.pos)
    if tpt.strand == '+':
        gnuc_beg = gnuc_del_beg
        gnuc_end = gnuc_del_end
    else:
        gnuc_beg = gnuc_del_end
        gnuc_end = gnuc_del_beg
    r.gnuc_range = '%d_%ddel' % (gnuc_beg, gnuc_end)
    r.pos = '%s:%d-%d' % (r.chrm, gnuc_beg, gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
    # print beg_codon, beg_codon.seq, end_codon, end_codon.seq
    
def nuc_mutation_del_coding_inframe(args, q, tpt, r):

    if q.beg.pos % 3 == 1:
        nuc_mutation_del_coding_inframe_inphase(args, q, tpt, r)
    else:
        nuc_mutation_del_coding_inframe_outphase(args, q, tpt, r)

def nuc_mutation_del_coding_frameshift(args, q, tpt, r):

    # assume frame-shift does not affect splicing
    beg_codon_index = (q.beg.pos + 2) / 3
    beg_codon_beg = beg_codon_index*3 - 2
    old_seq = tpt.seq[beg_codon_beg-1:]
    new_seq = tpt.seq[beg_codon_beg-1:q.beg.pos-1]+tpt.seq[q.end.pos:]
    taa_pos = None
    termlen = None
    for i in xrange(len(new_seq)/3):
        taa_ref_run = standard_codon_table[old_seq[3*i:3*i+3]]
        taa_alt_run = standard_codon_table[new_seq[3*i:3*i+3]]
        # print i, old_seq[3*i:3*i+3], new_seq[3*i:3*i+3], taa_ref_run, taa_alt_run, taa_pos
        if taa_pos == None and taa_ref_run != taa_alt_run:
            taa_pos = i
            taa_ref = taa_ref_run
            taa_alt = taa_alt_run
        if taa_alt_run == '*':
            if taa_pos == None:
                err_die('Terminating codon encountered before difference.', __name__)
                return None
            termlen = i + 1 - taa_pos
            break
    if termlen == None:
        err_die('No terminating codon before the end of the new transcript.', __name__)
        return None
    taa_pos += beg_codon_index
    r.taa_range = '%s%d%sfs*%d' % (taa_ref, taa_pos, taa_alt, termlen)
    r.tnuc_range = '%d_%ddel' % (q.beg.pos, q.end.pos)
    gnuc_beg, gnuc_end = tpt.tnuc_range2gnuc_range(q.beg.pos, q.end.pos)
    r.gnuc_range = '%ddel' % gnuc_beg if gnuc_beg == gnuc_end else '%d_%ddel' % (gnuc_beg, gnuc_end)
    natdelseq = tpt.seq[q.beg.pos-1:q.end.pos]
    r.info = 'NatDelSeq=%s' % natdelseq
    r.info += ';RefDelSeq=%s' % (natdelseq if tpt.strand == '+' else reverse_complement(natdelseq), )
    r.pos = '%s:%d-%d' % (r.chrm, gnuc_beg, gnuc_end)
    r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
    

def nuc_mutation_del_coding(args, q, tpt, r):

    if (q.end.pos - q.beg.pos) % 3 == 2: # in-frame
        nuc_mutation_del_coding_inframe(args, q, tpt, r)
    else:   # frame-shift
        nuc_mutation_del_coding_frameshift(args, q, tpt, r)

def nuc_mutation_del_intronic(args, q, tpt, r):

    # deletion occurs entirely in non-coding region
    # need only find genomic location of the deletion
    if q.beg.pos == q.end.pos: # with respect to one exome boundary
        codon = tpt.cpos2codon((q.beg.pos+2)/3)
        if not codon: return None
        i = q.beg.pos - (codon.index-1)*3 - 1
        print q.beg.pos, q.end.pos, codon
        if q.beg.tpos > 0:
            if tpt.strand == '+':
                if codon.locs[i] + 1 == codon.locs[i+1]: return False
                gnuc_beg = codon.locs[i] + q.beg.tpos
                gnuc_end = codon.locs[i] + q.end.tpos
                pl = []
                s = '-'.join(map(str, codon.locs[:i+1]))
                if s: pl.append(s)
                pl.append('(%d-%d)' % (gnuc_beg, gnuc_end))
                s = '-'.join(map(str, codon.locs[i+1:]))
                if s: pl.append(s)
                r.pos = '-'.join(pl)

            elif tpt.strand == '-':
                ir = 2-i
                if codon.locs[ir-1]+1 == codon.locs[ir]: return False
                gnuc_beg = codon.locs[ir] - q.end.tpos
                gnuc_end = codon.locs[ir] - q.beg.tpos
                pl = []
                s = '-'.join(map(str, codon.locs[:ir]))
                if s: pl.append(s)
                pl.append('(%d-%d)' % (gnuc_beg, gnuc_end))
                s = '-'.join(map(str, codon.locs[ir:]))
                if s: pl.append(s)
                r.pos = '-'.join(pl)

        elif q.beg.tpos < 0:
            if tpt.strand == '+':
                if codon.locs[i-1]+1 == codon.locs[i]: return False
                gnuc_beg = codon.locs[i] + q.beg.tpos
                gnuc_end = codon.locs[i] + q.end.tpos
                pl = []
                s = '-'.join(map(str, codon.locs[:i]))
                if s: pl.append(s)
                pl.append('(%d-%d)' % (gnuc_beg, gnuc_end))
                s = '-'.join(map(str, codon.locs[i:]))
                if s: pl.append(s)
                r.pos = '-'.join(pl)

            elif tpt.strand == '-':
                ir = 2-i
                if codon.locs[ir]+1 == codon.locs[ir+1]: return False
                gnuc_beg = codon.locs[ir] - q.end.tpos
                gnuc_end = codon.locs[ir] - q.beg.tpos
                pl = []
                s = '-'.join(map(str, codon.locs[:ir+1]))
                if s: pl.append(s)
                pl.append('(%d-%d)' % (gnuc_beg, gnuc_end))
                s = '-'.join(map(str, codon.locs[ir+1:]))
                if s: pl.append(s)
                r.pos = '-'.join(pl)

        r.reg = '%s (%s intronic)' % (tpt.gene.name, tpt.strand)
        r.gnuc_range = '%d_%ddel' % (gnuc_beg, gnuc_end)
        r.tnuc_range = '%d%d_%d%ddel' % (q.beg.pos, q.beg.tpos, q.end.pos, q.end.tpos)
        refdelseq = faidx.refgenome.fetch_sequence(tpt.chrm, gnuc_beg, gnuc_end)
        if q.delseq:
            if tpt.strand == '+' and refdelseq != q.delseq:
                raise IncompatibleTranscriptError()
            elif tpt.strand == '-' and refdelseq != reverse_complement(q.delseq):
                raise IncompatibleTranscriptError()
        natdelseq = refdelseq if tpt.strand == '+' else reverse_complement(refdelseq)
        r.info = 'RefDelSeq=%s;NatDelSeq=%s' % (refdelseq, natdelseq)
    else:
        err_die('Non-coding deletion range. not implemented yet')


def nuc_mutation_del(args, q, tpt):

    if q.tpt and tpt.name != q.tpt: return None
    tpt.ensure_seq()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name
    r.muttype = 'del'

    if q.beg.tpos == 0 and q.end.tpos == 0:
        nuc_mutation_del_coding(args, q, tpt, r)
    elif q.beg.tpos != 0 and q.end.tpos != 0:
        nuc_mutation_del_intronic(args, q, tpt, r)
    else:
        # one of the deletion start and end is in coding, the other in non-coding
        err_die('Mixing coding and non-coding, not implemented yet')

    return r

def _core_annotate_nuc_del(args, q, tpts):

    found = False
    for tpt in tpts:
        r = nuc_mutation_del(args, q, tpt)
        if r:
            found = True
            r.format(q.op)

    if not found:
        r = Record()
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

    return
