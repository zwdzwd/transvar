"""
The MIT License

Copyright (c) 2016
Wanding Zhou, Van Andel Research Institute

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

def annotate_region_cdna_transcript1(args, q, t, db):

    """Annotate cDNA region based on one transcript

    Args:
        args (argparse.Namespace): command line arguments
        q (record.QueryREG): query of region
        t (transcript.Transcript): transcript
        db (annodb.AnnoDB): annotation database

    Returns:
        r (record.Record): record

    Raises:
        IncompatibleTranscriptError: transcript is incompatible
    """

    ## checks
    # check transcript name if it is given
    if q.tpt and t.name != q.tpt:
        raise IncompatibleTranscriptError('Transcript name unmatched')
    # check q.beg and q.end is a valid Pos w.r.t exon boundaries
    t.check_exon_boundary(q.beg)
    t.check_exon_boundary(q.end)

    # transcript info
    r = Record()
    r.chrm = t.chrm
    r.tname = t.format()
    r.gene = t.gene_name
    r.strand = t.strand

    # region
    gnuc_beg = t.tnuc2gnuc(q.beg)
    gnuc_end = t.tnuc2gnuc(q.end)
    r.gnuc_beg = min(gnuc_beg, gnuc_end)
    r.gnuc_end = max(gnuc_beg, gnuc_end)
    r.reg = describe_genic(args, t.chrm, r.gnuc_beg, r.gnuc_end, t, db)
    expt = r.set_splice('included')

    # reference
    r.refrefseq = faidx.refgenome.fetch_sequence(t.chrm, r.gnuc_beg, r.gnuc_end)
    r.natrefseq = reverse_complement(r.refrefseq) if t.strand == '-' else r.refrefseq
    if q.refseq and r.natrefseq != q.refseq:
        raise IncompatibleTranscriptError()

    # g-syntax
    if r.gnuc_beg != r.gnuc_end:
        r.gnuc_range = '%d_%d%s' % (r.gnuc_beg, r.gnuc_end, r.refrefseq)
    else:
        r.gnuc_range = '%d%s' % (r.gnuc_beg, r.refrefseq)

    # c-syntax
    if q.beg != q.end:
        r.tnuc_range = '%s_%s%s' % (q.beg, q.end, r.natrefseq)
    else:
        r.tnuc_range = '%s%s' % (q.beg, r.natrefseq)

    # p-syntax
    if hasattr(r.reg, 'cover_cds') and r.reg.cover_cds:
        c1, p1 = t.intronic_lean(q.beg, 'c_greater')
        c2, p2 = t.intronic_lean(q.end, 'c_smaller')

        if c1.index == c2.index:
            r.taa_pos = c1.index
            r.taa_ref = aaf(codon2aa(c1.seq), args)
        else:
            taa_ref = aaf(translate_seq(t.seq[c1.index*3-3:c2.index*3]), args)
            r.taa_range = '%d_%d%s' % (c1.index, c2.index, aaf(taa_ref, args))

    return r

def annotate_region_cdna(args, q, tpts, db):

    """Annotate based on all transcripts

    Print or return records

    Args:
        args (argparse.Namespace): command line arguments
        q (record.QueryREG): query of region
        tpts (list of transcript.Transcript): transcripts
        db (annodb.AnnoDB): annotation database

    Returns:
        records (list of record.Record): a list of records

    Raises:
        IncompatibleTranscriptError: transcript is incompatible
    """

    records = []
    for t in tpts:
        try:
            r = annotate_region_cdna_transcript1(args, q, t, db)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue
        records.append(r)

    format_records(records, q.op, args)

    return records

def annotate_region_protein_transcript1(args, q, t, db):

    """Annotate region based on one transcript

    Args:
        args (argparse.Namespace): command line arguments
        q (record.QueryREG): query of region
        t (transcript.Transcript): transcript
        db (annodb.AnnoDB): annotation database

    Returns:
        r (record.Record): record

    Raises:
        IncompatibleTranscriptError: transcript is incompatible
    """

    # reference
    tnuc_beg = q.beg*3 - 2
    tnuc_end = q.end*3
    natrefseq = t.getseq(tnuc_beg, tnuc_end)
    refrefseq = reverse_complement(natrefseq) if t.strand == '-' else natrefseq
    taa_natrefseq = translate_seq(natrefseq)

    ## checks
    if q.tpt and t.name != q.tpt:
        raise IncompatibleTranscriptError('Transcript name unmatched')
    if q.end*3 > t.cdslen():
        raise IncompatibleTranscriptError('codon nonexistent')
    if q.beg_aa and q.beg_aa != taa_natrefseq[0]:
        raise IncompatibleTranscriptError('beginning reference amino acid unmatched')
    if q.end_aa and q.end_aa != taa_natrefseq[-1]:
        raise IncompatibleTranscriptError('ending reference amino acid unmatched')
    if q.refseq and not re.match(q.refseq.replace('x','[A-Z]'), taa_natrefseq):
        raise IncompatibleTranscriptError('reference sequence unmatched')

    # transcript info
    r = Record()
    r.chrm = t.chrm
    r.tname = t.format()
    r.gene = t.gene_name
    r.strand = t.strand

    # region
    r.reg = RegCDSAnno(t)
    r.reg.from_taa_range(q.beg, q.end)

    # g-syntax
    r.gnuc_beg, r.gnuc_end = t.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
    r.gnuc_range = '%d_%d' % (r.gnuc_beg, r.gnuc_end)

    # c-syntax
    r.tnuc_range = '%d_%d' % (tnuc_beg, tnuc_end)

    # p-syntax
    r.taa_range = '%s%d_%s%d' % (aaf(taa_natrefseq[0], args), q.beg, aaf(taa_natrefseq[-1], args), q.end) if q.beg != q.end else '%d%s' % (q.beg, aaf(taa_natrefseq[0], args))

    # info
    r.append_info('protein_sequence=%s;cDNA_sequence=%s;gDNA_sequence=%s' % (
        printseq(taa_natrefseq, args), printseq(natrefseq, args), printseq(refrefseq, args)))

    return r

def annotate_region_protein(args, q, tpts, db):

    """Annotation protein region on all transcripts

    Print or return records

    Args:
        args (argparse.Namespace): command line arguments
        q (record.QueryREG): query of region
        tpts (lists of transcript.Transcript): transcripts
        db (annodb.AnnoDB): annotation database

    Returns:
        records (list of record.Record): a list of records

    Raises:
        IncompatibleTranscriptError: transcript is incompatible
    """

    records = []
    for t in tpts:
        try:
            r = annotate_region_protein_transcript1(args, q, t, db)
        except IncompatibleTranscriptError:
            continue
        except SequenceRetrievalError:
            continue
        records.append(r)

    format_records(records, q.op, args)

    return records

def annotate_region_gdna_intergenic_point(args, q, reg):

    """annotate gDNA intergenic point

    Args:
        args (argparse.Namespace): command line arguments
        q (record.QueryREG): query of region
        reg (RegAnno): point region

    Returns:
        r (record.Record): record
    """

    r = Record()
    r.reg = reg
    r.chrm = q.tok
    r.set_promoter()

    r.gnuc_pos = q.pos if hasattr(q,'pos') else q.beg
    r.pos = r.gnuc_pos

    # # # annotate extra noncoding features
    # if 'GENCODE' in args.ffhs:
    #     iis = set()
    #     for entry in args.ffhs['GENCODE'].fetch(tok, beg, end+1):
    #         fields = entry.strip().split('\t')
    #         info = dict(re.findall(r'\s*([^"]*) "([^"]*)";', fields[8]))
    #         if 'gene_type' in info:
    #             ii = info['gene_type']
    #             if 'gene_name' in info: ii += '(%s)' % info['gene_name']
    #             iis.add(ii)
    #     r.info = 'gene_type=%s;' % ','.join(list(iis))
    return r

def annotate_region_gdna_genic_point(args, q, reg):

    """annotate gDNA genic point

    Args:
        args (argparse.Namespace): command line arguments
        q (record.QueryREG): query of region
        reg (RegAnno): point region

    Returns:
        r (record.Record): record
    """
    r = Record()
    r.reg = reg
    r.chrm = q.tok
    r.set_promoter()

    c, p = reg.t.gpos2codon(q.pos)
    r.append_info("is_gene_body")
    r.tname = reg.t.format()
    r.gene = reg.t.gene_name if reg.t.gene_name else '.'
    r.strand = reg.t.strand

    if p.tpos == 0 and reg.t.transcript_type == 'protein_coding':
        if c.seq in standard_codon_table:
            r.taa_ref = aaf(standard_codon_table[c.seq], args)
        r.taa_pos = c.index
        if args.aacontext>0 and r.taa_ref:
            aa1 = aaf(reg.t.taa_range2aa_seq(
                c.index-args.aacontext if c.index>=args.aacontext else 0, c.index-1), args)
            aa2 = aaf(reg.t.taa_range2aa_seq(c.index+1, c.index+args.aacontext), args)
            r.append_info('aacontext=%s[%s]%s' % (aa1, r.taa_ref, aa2))

    r.gnuc_pos = q.pos
    r.pos = q.pos
    r.gnuc_ref = faidx.refgenome.fetch_sequence(q.tok, q.pos, q.pos)
    r.tnuc_pos = p
    r.tnuc_ref = r.gnuc_ref if c.strand == '+' else complement(r.gnuc_ref)
    r.append_info('codon_pos=%s' % ('-'.join(map(str, c.locs)),))

    return r

def annotate_region_gdna_genic_span(args, q, reg):

    """annotate gDNA genic span

    Args:
        args (argparse.Namespace): command line arguments
        q (record.QueryREG): query of region
        reg (RegSpanAnno): span region

    Returns:
        r (record.Record): record
    """

    r = Record()
    r.tname = reg.t.format()
    # r.pos = '%d-%d' % (q.beg, q.end)
    r.chrm = q.tok
    r.reg = reg
    r.set_promoter()
    r.gene = reg.t.gene_name if reg.t.gene_name else '.'
    r.strand = reg.t.strand

    r.gnuc_range = '%d_%d' % (q.beg, q.end)

    r.set_splice('included')

    c1, p1 = reg.t.gpos2codon(q.beg)
    c2, p2 = reg.t.gpos2codon(q.end)

    if reg.t.strand == '+':
        r.tnuc_range = '%s_%s' % (p1,p2)
    else:
        r.tnuc_range = '%s_%s' % (p2,p1)

    # if protein coding transcript and not purely
    # intronic region, set amino acid
    if reg.t.transcript_type == 'protein_coding' and not same_intron(p1, p2):
        c1, p1 = reg.t.intronic_lean(p1, 'g_greater')
        c2, p2 = reg.t.intronic_lean(p2, 'g_smaller')

        if len(c1.seq) != 3 or len(c2.seq) != 3:
            r.append_info("truncated_refseq_at_boundary_(start_codon_seq_%s_and_end_codon_seq_%s)" % (c1.seq, c2.seq))
        elif c1.index == c2.index:
            # print 1
            r.taa_ref = aaf(c1.aa(), args)
            r.taa_pos = c1.index
            r.append_info('codon=%s' % c1.locformat())
        else:
            # print 2
            if reg.t.strand == '+':
                r.taa_range = '%s%d_%s%d' % (
                    aaf(c1.aa(), args), c1.index, c2.aa(), c2.index)
            else:
                r.taa_range = '%s%d_%s%d' % (
                    aaf(c2.aa(), args), c2.index, c1.aa(), c1.index)
            r.append_info('start_codon=%s;end_codon=%s' % (c1.locformat(), c2.locformat()))
    return r

def annotate_region_gdna_intergenic_span(args, q, reg):

    """annotate gDNA intergenic span
    Print or return records

    Args:
        args ((argparse.Namespace): command line arguments
        q (record.QueryREG): query of region
        reg (RegSpanAnno): span region

    Returns:
        r (record.Record): record
    """
    r = Record()
    r.reg = reg
    r.chrm = q.tok
    r.set_promoter()
    r.gnuc_range = '%d_%d' % (q.beg, q.end)
    # r.pos = '%d-%d' % (q.beg, q.end)
    r.set_splice('included')

    if hasattr(reg.b1, 't') and hasattr(reg.b2, 't'):
        if reg.b1.t == reg.b2.t:
            r.tname = reg.b1.t.format()
            r.gene = reg.b1.t.gene_name
            r.strand = reg.b1.t.strand
        else:
            r.tname = '%s,%s' % (reg.b1.t.format(), reg.b2.t.format())
            r.gene = '%s,%s' % (reg.b1.t.gene_name, reg.b2.t.gene_name)
            r.strand = '%s,%s' % (reg.b1.t.strand, reg.b2.t.strand)
    elif hasattr(reg.b1, 't'):
        r.tname = reg.b1.t.format()+','
        r.gene = reg.b1.t.gene_name+','
        r.strand = reg.b1.t.strand
    elif hasattr(reg.b2, 't'):
        r.tname = ','+reg.b2.t.format()
        r.gene = reg.b2.t.gene_name+','
        r.strand = reg.b2.t.strand

    return r

def annotate_region_gdna(args, q, db):

    """Annotate gDNA region
    Print or return records

    Args:
        args (argparse.Namespace): command line arguments
        q (record.QueryREG): query of region
        db (annodb.AnnoDB): annotation database

    Returns:
        records (list of record.Record): a list of records

    Raises:
        IncompatibleTranscriptError: transcript is incompatible
    """

    normalize_reg(q)
    records = []
    for reg in describe(args, q, db):
        # check transcript ID if given
        if q.tpt and hasattr(reg, 't') and reg.t.name != q.tpt:
            continue
        if isinstance(reg, RegAnno):       # point annotation
            if hasattr(reg, 'intergenic'): # intergenic
                r = annotate_region_gdna_intergenic_point(args, q, reg)
            elif hasattr(reg, 't'):        # genic
                r = annotate_region_gdna_genic_point(args, q, reg)
            else:
                raise Exception()          # shouldn't reach here

        elif isinstance(reg, RegSpanAnno): # span annotation
            if hasattr(reg, 't'):          # genic
                r = annotate_region_gdna_genic_span(args, q, reg)
            else:                          # intergenic
                r = annotate_region_gdna_intergenic_span(args, q, reg)
        else:
            raise Exception()              # shouldn't reach

        records.append(r)

    format_records(records, q.op, args)

    return records

def annotate_gene(args, q, tpts, db):

    """annotate gene
    """

    records = []
    for t in tpts:
        r = Record()
        r.chrm = t.chrm
        r.gene = q.gene.name
        r.gnuc_range = '%d_%d' % (t.beg, t.end)
        r.tnuc_range = '%d_%d' % (1, t.cdslen())
        if t.strand == '+':
            tss_tes = (t.chrm, t.beg-args.prombeg, t.beg+args.promend)
        else:
            tss_tes = (t.chrm, t.end-args.promend, t.end+args.prombeg)
        r.append_info('promoter=%s:%d_%d' % tss_tes)
        r.append_info('#exons=%d' % len(t.exons))
        if t.transcript_type == 'protein_coding':
            r.append_info('cds=%s:%d_%d' % (t.chrm, t.cds_beg, t.cds_end))
            r.taa_range = '%s%d_%s%d' % (aaf(t.taa2aa(1), args), 1, aaf(t.taa2aa(t.cdslen()//3), args), t.cdslen()//3)
            if args.pp or args.ppp:
                r.append_info('ref_protein_seq=%s' % aaf(t.get_proteinseq(), args))
            if t.cdslen() % 3 != 0:
                r.append_info('truncated_refseq_at_boundary')

        r.strand = t.strand
        r.tname = t.format()
        r.reg = 'whole_transcript'

        records.append(r)

    format_records(records, q.op, args)
    # format_records(r, rs, q.op, args)
    # format_all(rs, q.op, args)

    return records
