"""
The MIT License

Copyright (c) 2016
Wanding Zhou (zhouwanding@gmail.com)

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
from .transcripts import *
from .utils import *
from .record import *
from .describe import *
from .err import *
from .proteinseqs import *

#############################
###### cDNA annotation ######
#############################

def annotate_snv_cdna(args, q, tpts, db):

    """annotate cDNA SNV

    Args:
        args: (argparse.Namespace): command line arguments
        q (record.QuerySNV): query of a SNV
        tpts (a list of transcript.Transcript): transcripts to consider
        db (annoDB.AnnoDB): annotation database

    Returns:
        records (list of record.Record): a list of records
    """

    found = False
    records = []
    gene_name = ""
    for t in tpts:
        qpos = t.tnuc_resolve_pos(q.pos)
        try:
            if q.tpt and t.name != q.tpt:
                raise IncompatibleTranscriptError(
                    'transcript_id_unmatched_%s;expect_%s' % (q.tpt, t.name))
            t.ensure_seq()

            if (qpos.pos <= 0 or qpos.pos > t.cdslen()):
                raise IncompatibleTranscriptError(
                    'invalid_cDNA_position_%d;expect_[0_%d]' % (qpos.pos, t.cdslen()))
            codon = t.cpos2codon((qpos.pos+2)//3)
            if not codon:
                raise IncompatibleTranscriptError('invalid_cDNA_position_%d' % qpos.pos)

            r = Record(is_var=True)
            r.chrm = t.chrm
            r.tname = t.format()
            r.gene = t.gene_name
            gene_name = t.gene_name
            r.strand = t.strand

            r.gnuc_pos = t.tnuc2gnuc(qpos)
            r.gnuc_ref = faidx.refgenome.fetch_sequence(t.chrm, r.gnuc_pos, r.gnuc_pos)
            if t.strand == '+':
                if q.ref and r.gnuc_ref != q.ref:
                    raise IncompatibleTranscriptError(
                        'invalid_reference_%s;expect_%s' % (q.ref, r.gnuc_ref))
                r.gnuc_alt = q.alt if q.alt else ''
            else:
                if q.ref and r.gnuc_ref != complement(q.ref):
                    raise IncompatibleTranscriptError(
                        'invalid_reference_%s;expect_%s' % (q.ref, r.gnuc_ref))
                r.gnuc_alt = complement(q.alt) if q.alt else ''

            r.tnuc_pos = q.pos
            r.tnuc_ref = r.gnuc_ref if t.strand == '+' else complement(r.gnuc_ref)
            r.tnuc_alt = q.alt

            db.query_dbsnp(r, r.gnuc_pos, r.gnuc_ref, r.gnuc_alt)
            r.reg = describe_genic_site(args, t.chrm, r.gnuc_pos, t, db)
            # coding region
            if qpos.tpos == 0 and t.transcript_type == 'protein_coding':

                # Incompatible transcript
                if (q.ref and q.ref != t.seq[qpos.pos-1]):
                    raise IncompatibleTranscriptError(
                        'invalid_reference_%s;expect_%s' % (q.ref, t.seq[qpos.pos-1]))

                r.taa_ref = aaf(codon2aa(codon.seq), args)
                r.taa_pos = codon.index
                if not q.alt:
                    r.taa_alt = ''
                else:
                    mut_seq = list(codon.seq[:])
                    mut_seq[(qpos.pos-1) % 3] = q.alt
                    r.taa_alt = aaf(codon2aa(''.join(mut_seq)), args)
                    if r.taa_ref != r.taa_alt:
                        if r.taa_alt == '*':
                            r.csqn.append('Nonsense')
                        else:
                            r.csqn.append('Missense')
                    elif r.taa_alt:
                        r.csqn.append('Synonymous')
                    r.append_info(
                        'reference_codon=%s;alternative_codon=%s' % (
                            codon.seq, ''.join(mut_seq)))

                    variant_protein_seq_sub(
                        r, t, args, r.taa_pos, r.taa_pos, r.taa_alt)

            else:  # coordinates are with respect to the exon boundary
                r.csqn.append(r.reg.csqn()+"SNV")
                t.check_exon_boundary(qpos)

        except IncompatibleTranscriptError as e:
            continue
        except SequenceRetrievalError as e:
            continue
        records.append(r)

    format_records(records, q.op, args)
    return records


################################
###### protein annotation ######
################################

def annotate_snv_protein(args, q, tpts, db):

    """annotate protein SNV

    Args:
        args: (argparse.Namespace): command line arguments
        q (record.QuerySNV): query of a SNV
        tpts (a list of transcript.Transcript): transcripts to consider
        db (annoDB.AnnoDB): annotation database

    Returns:
        records (list of record.Record): a list of records
    """
    records = []
    for t in tpts:
        try:
            r, c = _annotate_snv_protein(args, q, t, db)
        except IncompatibleTranscriptError as e:
            continue
        except SequenceRetrievalError as e:
            err_warn(e)
            continue

        r.gene = t.gene_name
        r.strand = t.strand
        # alternative should be based on actual alternative aa determined
        set_taa_snv(r, q.pos, q.ref, r.taa_alt, args)
        r.reg = RegCDSAnno(t, c)

        records.append(r)

    format_records(records, q.op, args)
    return records


def _annotate_snv_protein(args, q, t, db):

    """ find all the mutations given a codon position, yield records """

    if q.alt and q.alt not in reverse_codon_table:
        err_warn('unknown alternative: %s, ignore alternative.' % q.alt)
        q.alt = ''

    # when there's a transcript specification
    if q.tpt and t.name != q.tpt:
        raise IncompatibleTranscriptError('transcript_id_unmatched_%s;expect_%s' % (q.tpt, t.name))

    t.ensure_seq()

    if (q.pos <= 0 or q.pos > t.cdslen()):
        raise IncompatibleTranscriptError('codon_nonexistent_%d;expect_[0:%d]' % (q.pos, t.cdslen()))
    codon = t.cpos2codon(q.pos)
    if not codon:
        raise IncompatibleTranscriptError('codon_nonexistent_%d;expect_[0:%d]' % (q.pos, t.cdslen()))

    # skip if reference amino acid is given
    # and codon sequence does not generate reference aa
    # codon.seq is natural sequence
    if q.ref and codon.seq not in aa2codon(q.ref):
        raise IncompatibleTranscriptError('invalid_reference_amino_acid_%s;expect_%s' % (q.ref, codon2aa(codon.seq)))

    r = Record(is_var=True)
    r.chrm = t.chrm
    r.tname = t.format()

    # if alternative amino acid is given
    # filter the target mutation set to those give
    # the alternative aa

    r.taa_alt = q.alt
    if q.alt:
        # here we assume a variant when an alternative aa is given,
        # even it might be the same as the reference
        tgt_codon_seqs = [x for x in aa2codon(q.alt) if x != codon.seq]
        # tgt_codon_seqs = aa2codon(q.alt)
        diffs = [codondiff(x, codon.seq) for x in tgt_codon_seqs]
        diffinds = sorted(list(range(len(diffs))), key=lambda i: len(diffs[i]))

        # guessed mutation
        gi = diffinds[0]        # guessed diff index
        gdiff = diffs[gi]       # guessed diff
        gtgtcodonseq = tgt_codon_seqs[gi]
        r.taa_alt = codon2aa(gtgtcodonseq)
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
            r.tnuc_range = '%d_%ddel%sins%s' % (tnuc_beg, tnuc_end, tnuc_ref, tnuc_alt)
            if codon.strand == '+':
                r.gnuc_range = '%d_%ddel%sins%s' % (codon.locs[gdiff[0]], codon.locs[gdiff[-1]],
                                                    tnuc_ref, tnuc_alt)
            else:
                r.gnuc_range = '%d_%ddel%sins%s' % (codon.locs[2-gdiff[-1]],
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
                    gnuc_tok  = '%s:g.%d%s>%s' % (t.chrm, codon.locs[diff[0]],
                                                  nrefbase, naltbase)
                else:
                    gnuc_tok = '%s:g.%d%s>%s' % (t.chrm, codon.locs[2-diff[0]],
                                                 complement(nrefbase), complement(naltbase))
                cdd_snv_muts.append(gnuc_tok)
            else:
                tnuc_beg = (codon.index-1)*3 + 1 + diff[0]
                tnuc_end = (codon.index-1)*3 + 1 + diff[-1]
                tnuc_ref = codon.seq[diff[0]:diff[-1]+1]
                tnuc_alt = tgtcodonseq[diff[0]:diff[-1]+1]
                tnuc_tok = 'c.%d_%ddel%sins%s' % (tnuc_beg, tnuc_end, tnuc_ref, tnuc_alt)
                if codon.strand == '+':
                    gnuc_tok = '%s:g.%d_%ddel%sins%s' % (t.chrm,
                                                         codon.locs[diff[0]],
                                                         codon.locs[diff[-1]],
                                                         tnuc_ref, tnuc_alt)
                else:
                    gnuc_tok = '%s:g.%d_%ddel%sins%s' % (t.chrm,
                                                         codon.locs[2-diff[-1]],
                                                         codon.locs[2-diff[0]],
                                                         reverse_complement(tnuc_ref),
                                                         reverse_complement(tnuc_alt))
                cdd_mnv_muts.append(gnuc_tok)

        r.append_info('reference_codon=%s;candidate_codons=%s' % (codon.seq, ','.join(tgt_codon_seqs)))
        if cdd_snv_muts:
            r.append_info('candidate_snv_variants=%s' % ','.join(cdd_snv_muts))
        if cdd_mnv_muts:
            r.append_info('candidate_mnv_variants=%s' % ','.join(cdd_mnv_muts))

        db.query_dbsnp_codon(r, codon, q.alt if q.alt else None)
        db.query_feature(r, r.chrm, codon.locs[0], codon.locs[2])
    else:
        r.gnuc_range = '%d_%d' % (codon.locs[0], codon.locs[2])
        r.tnuc_range = '%d_%d' % ((codon.index-1)*3+1, (codon.index-1)*3+3)

    return r, codon

# used by codon search
def __core_annotate_codon_snv(args, q, db):
    for t in q.gene.tpts:
        try:
            r, c = _annotate_snv_protein(args, q, t, db)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue
        yield t, c

#############################
###### gDNA annotation ######
#############################

def annotate_snv_gdna(args, q, db):

    """annotate gDNA SNV

    Args:
        args: (argparse.Namespace): command line arguments
        q (record.QuerySNV): query of a SNV
        db (annodb.AnnoDB): annotation database

    Returns:
        records (list of record.Record): a list of records
    """

    gnuc_ref = faidx.refgenome.fetch_sequence(q.tok, q.pos, q.pos)
    if not q.ref:
        q.ref = gnuc_ref

    if gnuc_ref != q.ref:
        wrap_exception(Exception("invalid_reference_base_%s;expect_%s" % (q.ref, gnuc_ref)))

    records = []
    for reg in describe(args, q, db):

        # skip if transcript ID does not match
        if q.tpt and hasattr(reg, 't') and reg.t.name != q.tpt:
            continue

        r = Record(is_var=True)
        r.reg = reg
        r.chrm = q.tok
        r.gnuc_pos = q.pos
        r.pos = r.gnuc_pos
        r.gnuc_ref = gnuc_ref
        r.gnuc_alt = q.alt if q.alt else ''
        db.query_dbsnp(r, q.pos, q.ref, q.alt if q.alt else None)

        if hasattr(reg, 't'):
            r = annotate_snv_gdna_trannscript(reg, r, q, args)
        else:
            r.csqn.append(r.reg.csqn()+"SNV")

        records.append(r)

    format_records(records, q.op, args)
    return records

def annotate_snv_gdna_trannscript(reg, r, q, args):

    """annotate gDNA SNV at transcript region

    Args:
        reg (RegSpanAnno): span region
        q (record.QuerySNV): query of SNV
        r (record.Record): record to be updated
        args (argparse.Namespace): command line arguments

    Return:
        r (record.Record): record
    """

    c,p = reg.t.gpos2codon(q.pos)

    r.tname = reg.t.format()
    r.gene = reg.t.gene_name
    r.strand = reg.t.strand
    r.tnuc_pos = p

    if c.strand == '+':
        r.tnuc_ref = r.gnuc_ref
        r.tnuc_alt = r.gnuc_alt
    else:
        r.tnuc_ref = complement(r.gnuc_ref)
        r.tnuc_alt = complement(r.gnuc_alt) if r.gnuc_alt else ''

    if r.gnuc_ref == r.gnuc_alt:
        csqn_action = "Synonymous"
    else:
        csqn_action = "SNV"

    if not r.set_splice("mutated", csqn_action):
        if p.tpos == 0 and reg.t.transcript_type=='protein_coding':
            if c.seq in standard_codon_table:
                r.taa_ref = aaf(standard_codon_table[c.seq], args)
                r.taa_pos = c.index

                if args.aacontext>0 and r.taa_ref:
                    aa1 = aaf(reg.t.taa_range2aa_seq(
                        c.index-args.aacontext if c.index>=args.aacontext else 0, c.index-1), args)
                    aa2 = aaf(reg.t.taa_range2aa_seq(c.index+1, c.index+args.aacontext), args)
                    r.append_info('aacontext=%s[%s]%s' % (aa1, r.taa_ref, aa2))

                if q.alt:
                    if c.strand == '+':
                        alt_seq = set_seq(c.seq, c.locs.index(q.pos), q.alt)
                    else:
                        alt_seq = set_seq(c.seq, 2-c.locs.index(q.pos), complement(q.alt))

                    r.taa_alt = aaf(codon2aa(alt_seq), args)
                    if r.taa_alt != r.taa_ref:
                        if r.taa_alt == '*':
                            r.csqn.append('Nonsense')
                        else:
                            r.csqn.append('Missense')
                    elif r.taa_alt:
                        r.csqn.append('Synonymous')
                    variant_protein_seq_sub(
                        r, reg.t, args, r.taa_pos, r.taa_pos, r.taa_alt)
            else:
                r.append_info('truncated_refseq_at_boundary_(codon_seq_%s_codon_index_%d_protein_length_%d)'
                              % (c.seq, c.index, reg.t.cdslen()//3))

            r.append_info('codon_pos=%s' % (c.locformat(),))
            r.append_info('ref_codon_seq=%s' % c.seq)
        else:
            r.csqn.append(r.reg.csqn()+"SNV")

    return r

def set_taa_snv(r, pos, ref, alt, args):

    r.taa_pos = pos
    r.taa_ref = aaf(ref, args)
    r.taa_alt = aaf(alt, args)
    if r.taa_ref != r.taa_alt:
        if r.taa_alt == '*':
            r.csqn.append('Nonsense')
        elif r.taa_alt:
            r.csqn.append('Missense')
        else:
            r.csqn.append('Unclassified')
    elif r.taa_ref:
        r.csqn.append('Synonymous')
