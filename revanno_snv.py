from transcripts import *
from utils import *
from record import *

def nuc_mutation_snv_coding(r, tpt, codon, q):

    r.reg = '%s (%s coding)' % (tpt.gene.name, tpt.strand)
    r.pos = '-'.join(map(str, codon.locs))
    r.tnuc_pos = q.cpos()
    r.tnuc_ref = q.ref
    r.tnuc_alt = q.alt

    if (q.ref and q.ref != tpt.seq[q.cpos()-1]):
        raise IncompatibleTranscriptError('SNV ref not matched')
    if codon.strand == '+':
        r.gnuc_pos = codon.locs[0]+(q.cpos()-1)%3
        r.gnuc_ref = q.ref if q.ref else codon.seq[(q.cpos()-1)%3]
        if q.alt: r.gnuc_alt = q.alt
    else:
        r.gnuc_pos = codon.locs[-1]-(q.cpos()-1)%3
        r.gnuc_ref = complement(q.ref if q.ref else codon.seq[(q.cpos()-1)%3])
        if q.alt: r.gnuc_alt = complement(q.alt)

    r.taa_ref = codon2aa(codon.seq)
    r.taa_pos = codon.index
    if not q.alt:
        r.taa_alt = ''
    else:
        mut_seq = list(codon.seq[:])
        mut_seq[(q.cpos()-1) % 3] = q.alt
        r.taa_alt = codon2aa(''.join(mut_seq))
        r.append_info('NCodonSeq=%s;NAltCodonSeq=%s' % (codon.seq, ''.join(mut_seq)))

def nuc_mutation_snv_intronic(r, tpt, codon, q):

    r.reg = '%s (%s intronic)' % (tpt.gene.name, tpt.strand)
    i = q.cpos() - (codon.index-1)*3 - 1
    np = tpt.position_array()
    check_exon_boundary(np, q.pos)
    r.gnuc_pos = tnuc2gnuc2(np, q.pos, tpt)
    r.gnuc_ref = faidx.refgenome.fetch_sequence(tpt.chrm, r.gnuc_pos, r.gnuc_pos)
    r.pos = r.gnuc_pos
    r.gnuc_ref = faidx.refgenome.fetch_sequence(tpt.chrm, r.gnuc_pos, r.gnuc_pos)
    if tpt.strand == '+':
        if q.ref and r.gnuc_ref != q.ref: raise IncompatibleTranscriptError()
        r.gnuc_alt = q.alt if q.alt else ''
    else:
        if q.ref and r.gnuc_ref != complement(q.ref): raise IncompatibleTranscriptError()
        r.gnuc_alt = complement(q.alt) if q.alt else ''
    r.tnuc_pos = q.pos
    r.tnuc_ref = r.gnuc_ref if tpt.strand == '+' else complement(r.gnuc_ref)
    r.tnuc_alt = q.alt

def nuc_mutation_snv(args, q, tpt):

    if q.tpt and tpt.name != q.tpt:
        raise IncompatibleTranscriptError('transcript id unmatched')
    tpt.ensure_seq()

    if (q.cpos() <= 0 or q.cpos() > len(tpt)):
        raise IncompatibleTranscriptError()
    codon = tpt.cpos2codon((q.cpos()+2)/3)
    if not codon:
        raise IncompatibleTranscriptError()

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name
    if tpt.gene.dbxref:
        r.info = 'DBXref=%s' % tpt.gene.dbxref

    if q.pos.tpos == 0:                # coding region
        nuc_mutation_snv_coding(r, tpt, codon, q)
    else:          # coordinates are with respect to the exon boundary
        nuc_mutation_snv_intronic(r, tpt, codon, q)

    return r

def _core_annotate_nuc_snv(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r = nuc_mutation_snv(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
        except UnknownChromosomeError:
            continue
        found = True
        r.format(q.op)

    if not found:
        r = Record()
        r.tnuc_pos = q.pos
        r.tnuc_ref = q.ref
        r.tnuc_alt = q.alt
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

    return

def codon_mutation_snv(args, q, tpt):

    """ find all the mutations given a codon position, yield records """

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

    # skip if reference amino acid is given
    # and codon sequence does not generate reference aa
    # codon.seq is natural sequence
    if q.ref and codon.seq not in aa2codon(q.ref):
        raise IncompatibleTranscriptError('reference amino acid unmatched')

    r = Record()
    r.chrm = tpt.chrm
    r.tname = tpt.name
    if tpt.gene.dbxref:
        r.info = 'DBXref=%s' % tpt.gene.dbxref
    r.pos = '-'.join(map(str, codon.locs))

    # if alternative amino acid is given
    # filter the target mutation set to those give
    # the alternative aa
    
    if q.alt:
        tgt_codon_seqs = [x for x in aa2codon(q.alt) if x != codon.seq]
        diffs = [codondiff(x, codon.seq) for x in tgt_codon_seqs]
        diffinds = sorted(range(len(diffs)), key=lambda i: len(diffs[i]))

        # guessed mutation
        gi = diffinds[0]        # guessed diff index
        gdiff = diffs[gi]       # guessed diff
        gtgtcodonseq = tgt_codon_seqs[gi]
        if len(gdiff) == 1:
            nrefbase = codon.seq[gdiff[0]]
            naltbase = gtgtcodonseq[gdiff[0]]

            r.tnuc_pos = (codon.index-1)*3 + 1 + gdiff[0]
            r.tnuc_ref = nrefbase
            r.tnuc_alt = naltbase
            if codon.strand == '+':
                r.gnuc_ref = nrefbase
                r.gnuc_alt = naltbase
                r.gnuc_pos = codon.locs[gdiff[0]]
            else:
                r.gnuc_ref = complement(nrefbase)
                r.gnuc_alt = complement(naltbase)
                r.gnuc_pos = codon.locs[2-gdiff[0]]
        else:
            tnuc_beg = (codon.index-1)*3 + 1 + gdiff[0]
            tnuc_end = (codon.index-1)*3 + 1 + gdiff[-1]
            tnuc_ref = codon.seq[gdiff[0]:gdiff[-1]+1]
            tnuc_alt = gtgtcodonseq[gdiff[0]:gdiff[-1]+1]
            r.tnuc_range = '%d_%d%s>%s' % (tnuc_beg, tnuc_end, tnuc_ref, tnuc_alt)
            if codon.strand == '+':
                r.gnuc_range = '%d_%d%s>%s' % (codon.locs[gdiff[0]], codon.locs[gdiff[-1]],
                                               tnuc_ref, tnuc_alt)
            else:
                r.gnuc_range = '%d_%d%s>%s' % (codon.locs[2-gdiff[-1]],
                                               codon.locs[2-gdiff[0]],
                                               reverse_complement(tnuc_ref), 
                                               reverse_complement(tnuc_alt))
        # candidate mutations
        cdd_snv_muts = []
        cdd_mnv_muts = []
        for i in diffinds:
            if i == gi: continue
            diff = diffs[i]
            tgtcodonseq = tgt_codon_seqs[i]
            if len(diff) == 1:
                nrefbase = codon.seq[diff[0]]
                naltbase = tgtcodonseq[diff[0]]
                tnuc_pos = (codon.index-1)*3 + 1 + diff[0]
                tnuc_tok = 'c.%d%s>%s' % (tnuc_pos, nrefbase, naltbase)
                if codon.strand ==  '+':
                    gnuc_tok  = '%s:g.%d%s>%s' % (tpt.chrm, codon.locs[diff[0]],
                                                  nrefbase, naltbase)
                else:
                    gnuc_tok = '%s:g.%d%s>%s' % (tpt.chrm, codon.locs[2-diff[0]],
                                                 complement(nrefbase), complement(naltbase))
                cdd_snv_muts.append(gnuc_tok)
            else:
                tnuc_beg = (codon.index-1)*3 + 1 + diff[0]
                tnuc_end = (codon.index-1)*3 + 1 + diff[-1]
                tnuc_ref = codon.seq[diff[0]:diff[-1]+1]
                tnuc_alt = tgtcodonseq[diff[0]:diff[-1]+1]
                tnuc_tok = 'c.%d_%d%s>%s' % (tnuc_beg, tnuc_end, tnuc_ref, tnuc_alt)
                if codon.strand == '+':
                    gnuc_tok = '%s:g.%d_%d%s>%s' % (tpt.chrm, 
                                                    codon.locs[diff[0]],
                                                    codon.locs[diff[-1]],
                                                    tnuc_ref, tnuc_alt)
                else:
                    gnuc_tok = '%s:g.%d_%d%s>%s' % (tpt.chrm,
                                                    codon.locs[2-diff[-1]],
                                                    codon.locs[2-diff[0]],
                                                    reverse_complement(tnuc_ref),
                                                    reverse_complement(tnuc_alt))
                cdd_mnv_muts.append(gnuc_tok)
                                                    
        r.append_info('NCodonSeq=%s;NCddSeqs=%s' % (codon.seq, ','.join(tgt_codon_seqs)))
        if cdd_snv_muts:
            r.append_info('CddSNVMuts=%s' % ','.join(cdd_snv_muts))
        if cdd_mnv_muts:
            r.append_info('CddMNVMuts=%s' % ','.join(cdd_mnv_muts))
        if args.dbsnp_fh:
            dbsnps = []
            for entry in args.dbsnp_fh.fetch(tpt.chrm, codon.locs[0]-3, codon.locs[0]):
                f = entry.split('\t')
                dbsnps.append('%s(%s:%s%s>%s)' % (f[2], f[0], f[1], f[3], f[4]))
            if dbsnps:
                r.append_info('DBSNP=%s' % ','.join(dbsnps))
    else:
        r.gnuc_range = '%d_%d' % (codon.locs[0], codon.locs[2])
        r.tnuc_range = '%d_%d' % ((codon.index-1)*3+1, (codon.index-1)*3+3)

    return r, codon

def __core_annotate_codon_snv(args, q):
    for tpt in q.gene.tpts:
        try:
            r, c = codon_mutation_snv(args, q, tpt)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue
        except UnknownChromosomeError:
            continue
        yield tpt, c

def _core_annotate_codon_snv(args, q, tpts):

    found = False
    for tpt in tpts:
        try:
            r, c = codon_mutation_snv(args, q, tpt)
        except IncompatibleTranscriptError as e:
            continue
        except SequenceRetrievalError as e:
            continue
        except UnknownChromosomeError as e:
            sys.stderr.write(str(e))
            continue
        r.taa_pos = q.pos
        r.taa_ref = q.ref
        r.taa_alt = q.alt
        r.reg = '%s (%s, coding, exon %s)' % (
            tpt.gene.name, tpt.strand,
            tpt.tnuc_range2exon_inds(c.index*3-2, c.index*3))
        r.format(q.op)
        found = True

    if not found:
        r = Record()
        r.taa_pos = q.pos
        r.taa_ref = q.ref
        r.taa_alt = q.alt
        r.info = 'NoValidTranscriptFound'
        r.format(q.op)

