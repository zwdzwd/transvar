from transcripts import *
from utils import *
from record import *

def codon_mutation_fs(args, q, tpt):

    if q.alt and q.alt not in reverse_codon_table:
        sys.stderr.write("Unknown alternative: %s, ignore alternative.\n" % q.alt)
        q.alt = ''

    # when there's a transcript specification
    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError('transcript id unmatched')

    tpt.ensure_seq()

    if (q.pos <= 0 or q.pos > len(tpt)):
        raise IncompatibleTranscriptError('codon nonexistent')
    codon = tpt.cpos2codon(q.pos)
    if not codon:
        raise IncompatibleTranscriptError('codon nonexistent')

    # if hasattr(q, 'beg_aa') and q.beg_aa and q.beg_aa != codon2aa(codon.seq):
    #     raise IncompatibleTranscriptError('Unmatched reference amino acid')

    # skip if reference amino acid is given
    # and codon sequence does not generate reference aa
    # codon.seq is natural sequence
    if q.ref and codon.seq not in aa2codon(q.ref):
        raise IncompatibleTranscriptError('reference amino acid unmatched')

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name
    r.pos = '-'.join(map(str, codon.locs))
    tnuc_beg = (codon.index-1)*3+1
    tnuc_end = (q.pos+q.stop_index)*3 # may be out of bound
    r.tnuc_range = '%d-%d' % (tnuc_beg, tnuc_end)
    if tnuc_end >= len(tpt):
        np = tpt.position_array()
        gnuc_beg = tnuc2gnuc(np, tnuc_beg)
        gnuc_end = tpt.cds_end + tnuc_end - len(tpt)
    else:
        gnuc_beg, gnuc_end = tpt.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
    r.gnuc_range = '%d-%d' % (gnuc_beg, gnuc_end)
    r.info = "RoughEstimateFromFrameShift"

    return r


def _core_annotate_codon_fs(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = codon_mutation_fs(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
        except UnknownChromosomeError:
            continue
        r.taa_range = '%s%d%sfs*%d' % (q.ref, q.pos, q.alt, q.stop_index)
        r.reg = '%s (%s, coding)' % (tpt.gene.name, tpt.strand)
        r.format(q.op)
        found = True

    if not found:
        r = Record('FS')
        r.taa_range = '%s%d%sfs*%d' % (q.ref, q.pos, q.alt, q.stop_index)
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

