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

from transcripts import *
from utils import *
from record import *
from err import *
from snv import annotate_snv_protein

import itertools
alphabets = 'ACGT'

def fuzzy_match_deletion(t, codon, q, args):

    matches = {}                # map right-aligned gDNA identifier to detailed information
    for ds in [1,2]:            # probe deletion of length 1,2
        i = codon.index*3
        for j in xrange(i, max(0, i-10), -1):
            jb = j/3*3
            old_seq = t.seq[jb:]
            new_seq = t.seq[jb:j]+t.seq[j+ds:]
            tnuc_delseq = t.seq[j:j+ds]
            aae = t.extend_taa_seq(j/3+1, old_seq, new_seq)
            if aae:
                if (q.ref == aae.taa_ref and ((not q.alt) or q.alt == aae.taa_alt)
                    and q.stop_index == aae.termlen and q.pos == aae.taa_pos):
                    t.ensure_position_array()
                    if t.strand == '+':
                        gnuc_beg, gnuc_end = t.np[j], t.np[j+ds-1]
                        gnuc_delseq = t.seq[j:j+ds]
                    else:
                        gnuc_beg, gnuc_end = t.np[j+ds-1], t.np[j]
                        gnuc_delseq = reverse_complement(t.seq[j:j+ds])
                    gnuc_beg_r, gnuc_end_r = gnuc_roll_right_del(t.chrm, gnuc_beg, gnuc_end)
                    gnuc_delseq_r = faidx.getseq(t.chrm, gnuc_beg_r, gnuc_end_r)
                    gnuc_id_r = gnuc_del_id(t.chrm, gnuc_beg_r, gnuc_end_r, args, gnuc_delseq=gnuc_delseq_r)
                    if gnuc_id_r not in matches:
                        # compute gDNA id left aligned
                        gnuc_beg_l, gnuc_end_l = gnuc_roll_left_del(t.chrm, gnuc_beg, gnuc_end)
                        gnuc_delseq_l = faidx.getseq(t.chrm, gnuc_beg_l, gnuc_end_l)
                        # compute cDNA level id, left and right aligned
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
                        matches[gnuc_id_r] = (
                            gnuc_del_id(t.chrm, gnuc_beg_l, gnuc_end_l, args, gnuc_delseq_l),
                            tnuc_del_id(p1r, p2r, args, tnuc_delseq_r),
                            tnuc_del_id(p1l, p2l, args, tnuc_delseq_l))
                        
    return matches

class FuzzyInsMatch():

    def __init__(self):
        self.insseq = ''
        self.termlen = -1
        self.tnuc_pos = -1

    def __repr__(self):
        return '<FMatch: %dfs*%dins%s>' % (self.tnuc_pos, self.termlen, self.insseq)

def fuzzy_match_insertion_aa_change(t, j, ins_len, q):
    """
    given insertion length and insertion location
    find all the insertion sequence that match taa_alt
    try to match the insertion length just for early stop

    j is the tnuc location of insertion
    ins_len is the nucleotide length of insertion
    """
    jb = j/3*3
    old_seq = t.seq[jb:]
    match_seq = []
    termlen_match = False
    for _insseq in itertools.product(alphabets, repeat=ins_len):
        insseq = ''.join(_insseq)
        new_seq = t.seq[jb:j]+insseq+t.seq[j:]
        aae = t.extend_taa_seq(j/3+1, old_seq, new_seq)
        if (aae and 
            (((not q.alt) or aae.taa_alt == q.alt) and q.ref == aae.taa_ref
             and q.pos == aae.taa_pos and aae.termlen <= q.stop_index)):
            m = FuzzyInsMatch()
            m.insseq = insseq
            m.termlen = aae.termlen
            match_seq.append(m)
            if m.termlen == q.stop_index:
                termlen_match = True

    # early stop when termlen also matches
    if termlen_match:
        match_seq = [m for m in match_seq if m.termlen == q.stop_index]
    return match_seq

def fuzzy_match_insertion_scan_loc(t, codon_index, ins_len, q):

    """ given insertion length, scan different insertion location """
    i = codon_index*3
    matched_seqs = []
    # prime 10 insertion location before the channged aa
    for j in xrange(i, max(0,i-10), -1):
        matched_seq = fuzzy_match_insertion_aa_change(t, j, ins_len, q)
        if matched_seq:
            for m in matched_seq:
                m.tnuc_pos = j
            matched_seqs.extend(matched_seq)

    return matched_seqs

def fuzzy_match_insertion_update_match(matches, t, insseq, j):

    t.ensure_position_array()

    # if this is the first match and there is ambiguous base 'N'
    # instantiate this by a normal base 'C' (since there 'C' is
    # not involved in the stop codon and still keeping the 'N' in
    # the candidate
    # make sure that there is at least 1 output with no 'N'
    if len(matches) == 0 and 'N' in insseq:
        fuzzy_match_insertion_update_match(matches, t, insseq.replace('N','C'), j)

    gnuc_insseq = insseq if t.strand == '+' else reverse_complement(insseq)
    # note that t.np is 0-based
    gnuc_ins = gnuc_set_ins_core(t.chrm, t.np[j-1] if t.strand == '+' else t.np[j], gnuc_insseq)
    tnuc_ins = tnuc_set_ins_core(gnuc_ins, t)
    matches[gnuc_ins.right_align()] = (
        gnuc_ins.left_align(),
        tnuc_ins.right_align(),
        tnuc_ins.left_align())

def fuzzy_match_insertion(t, codon, q):

    matches = {}

    # prime the first changed aa, this involves
    # at most 5 base non-free insertion (those do not
    #  interfere with reference bases to form new codons)
    mseqs = []
    for ds in xrange(1,6):
        if len(mseqs) > 3000:
            break
        exact_match_found = False
        mseq = fuzzy_match_insertion_scan_loc(t, codon.index, ds, q)
        if mseq and [m for m in mseq if m.termlen == q.stop_index]:
            for m in mseq:
                if m.termlen == q.stop_index:
                    fuzzy_match_insertion_update_match(matches, t, m.insseq, m.tnuc_pos)
                    exact_match_found = True
        mseqs.extend(mseq)
        if exact_match_found:
            break

    # early stop
    if matches:
        return matches

    # cannot match first aa change
    # TODO: enlarge the search span for locations instead of giving up
    if not mseqs:
        return {}

    max_termlen = max(m.termlen for m in mseqs)
    freefill_len = q.stop_index - max_termlen
    # determine additional insertion leading to the stop codon
    if freefill_len < 10:
        for i, m in enumerate(mseqs):
            if m.termlen == max_termlen:
                k = 3-m.tnuc_pos%3
                filled_insseq = m.insseq[:k]+'N'*(freefill_len*3)+m.insseq[k:]
                # print i
                fuzzy_match_insertion_update_match(matches, t, filled_insseq, m.tnuc_pos)

    return matches

def _annotate_frameshift(args, q, t):

    if q.alt and q.alt not in reverse_codon_table:
        err_warn('unknown alternative: %s, ignore alternative' % q.alt)
        q.alt = ''

    # when there's a transcript specification
    if q.tpt and t.name != q.tpt:
        raise IncompatibleTranscriptError('transcript id unmatched')

    t.ensure_seq()

    if (q.pos <= 0 or q.pos > t.cdslen()):
        raise IncompatibleTranscriptError('codon nonexistent')
    codon = t.cpos2codon(q.pos)
    if not codon:
        raise IncompatibleTranscriptError('codon nonexistent')

    # skip if reference amino acid is given
    # and codon sequence does not generate reference aa
    # codon.seq is natural sequence
    if q.ref and codon.seq not in aa2codon(q.ref):
        raise IncompatibleTranscriptError('reference amino acid unmatched')

    r = Record(is_var=True)
    r.chrm = t.chrm
    r.tname = t.format()
    r.gene = t.gene_name
    r.strand = t.strand
    r.pos = '-'.join(map(str, codon.locs))

    # backward_probe_deletion
    # print codon.index
    # print t.seq

    matches = fuzzy_match_deletion(t, codon, q, args)
    if not matches:
        matches = fuzzy_match_insertion(t, codon, q)
    if matches:
        gmatches = sorted(matches.keys())
        chosen = None

        # get a match without 'N'
        for i, gnuc_range in enumerate(gmatches):
            if 'N' not in gnuc_range:
                r.gnuc_range = gnuc_range
                chosen = i
                break

        gnuc_id_l, tnuc_id_r, tnuc_id_l = matches[r.gnuc_range]
        r.tnuc_range = tnuc_id_r
        r.append_info('left_align_cDNA=c.%s' % tnuc_id_l)
        r.append_info('left_align_gDNA=g.%s' % gnuc_id_l)
        if len(matches) > 1:
            cands = []
            for k in xrange(0,min(args.nc,len(gmatches))):
                if k == chosen:
                    continue
                gnuc_id_r = gmatches[k]
                gnuc_id_l, tnuc_id_r, tnuc_id_l = matches[gnuc_id_r]
                cands.append('g.%s/c.%s/g.%s/c.%s' % (gnuc_id_r, tnuc_id_r, gnuc_id_l, tnuc_id_l))
            r.append_info('candidates=%s' % ','.join(cands))
            if len(matches) > args.nc:
                r.append_info('%d_CandidatesOmitted' % (len(matches)-args.nc))
    else:
        tnuc_beg = (codon.index-1)*3+1
        tnuc_end = (q.pos+q.stop_index)*3 # may be out of bound
        r.tnuc_range = '(%d_%d)' % (tnuc_beg, tnuc_end)
        if tnuc_end >= t.cdslen():
            t.ensure_position_array()
            gnuc_beg = t._tnuc2gnuc(tnuc_beg)
            gnuc_end = t.cds_end + tnuc_end - t.cdslen()
        else:
            gnuc_beg, gnuc_end = t.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
        r.gnuc_range = '(%d_%d)' % (gnuc_beg, gnuc_end)
        r.append_info('imprecise')

    return r


def annotate_frameshift(args, q, tpts, db):

    if q.alt == '*':
        annotate_snv_protein(args, q, tpts, db)
        return

    found = False
    rs = []
    for t in tpts:
        try:
            r = _annotate_frameshift(args, q, t)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue

        r.taa_range = format_fs(q, args)
        r.csqn.append("Frameshift")
        r.reg = RegCDSAnno(t)
        r.reg.from_taa_range(q.pos, q.pos+(q.stop_index if q.stop_index>=0 else 0))
        found = True
        format_one(r, rs, q, args)

    format_all(rs, q, args)

    if not found:
        r = Record(is_var=True)
        r.taa_range = format_fs(q, args)
        r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))
        r.format(q.op)


def format_fs(q, args):
    
    if q.stop_index >= 0:
        return '%s%d%sfs*%d' % (aaf(q.ref, args), q.pos, aaf(q.alt, args), q.stop_index)
    else:
        return '%s%d%sfs' % (aaf(q.ref, args), q.pos, aaf(q.alt, args))
