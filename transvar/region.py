""" annotate region """
from err import *
from record import *
from transcripts import *
from describe import *

def _annotate_region_cdna(args, q, t, db):

    if q.tpt and t.name != q.tpt:
        raise IncompatibleTranscriptError('Transcript name unmatched')
    t.ensure_seq()

    r = Record()
    r.chrm = t.chrm
    r.tname = t.name
    r.gene = t.gene.name
    r.strand = t.strand

    t.check_exon_boundary(q.beg)
    t.check_exon_boundary(q.end)

    gnuc_beg = t.tnuc2gnuc(q.beg)
    gnuc_end = t.tnuc2gnuc(q.end)
    r.gnuc_beg = min(gnuc_beg, gnuc_end)
    r.gnuc_end = max(gnuc_beg, gnuc_end)
    tnuc_coding_beg = q.beg.pos if q.beg.tpos <= 0 else q.beg.pos+1 # base after intron
    tnuc_coding_end = q.end.pos if q.end.tpos >= 0 else q.end.pos-1 # base before intron
    gnuc_coding_beg = t._tnuc2gnuc(tnuc_coding_beg)
    gnuc_coding_end = t._tnuc2gnuc(tnuc_coding_end)

    r.refrefseq = faidx.refgenome.fetch_sequence(t.chrm, r.gnuc_beg, r.gnuc_end)
    r.natrefseq = reverse_complement(r.refrefseq) if t.strand == '-' else r.refrefseq
    if q.ref and r.natrefseq != q.ref:
        raise IncompatibleTranscriptError()

    r.gnuc_range = '%d_%d%s' % (r.gnuc_beg, r.gnuc_end, r.refrefseq) if r.gnuc_beg != r.gnuc_end else '%d%s' % (r.gnuc_beg, r.refrefseq)
    r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
    r.tnuc_range = '%s_%s%s' % (q.beg, q.end, r.natrefseq) if q.beg != q.end else '%s%s' % (q.beg, r.natrefseq)

    r.reg = describe_genic(args, t.chrm, r.gnuc_beg, r.gnuc_end, t, db)
    expt = r.set_splice('included')
    if hasattr(r.reg, 'cover_cds') and r.reg.cover_cds: #region_in_exon(np, q.beg, q.end):
        c1, p1 = t.intronic_lean(q.beg, 'c_greater')
        c2, p2 = t.intronic_lean(q.end, 'c_smaller')
        
        if c1.index == c2.index:
            r.taa_pos = c1.index
            r.taa_ref = codon2aa(c1.seq)
        else:
            r.taa_ref = translate_seq(t.seq[c1.index*3-3:c2.index*3])
            r.taa_range = '%d_%d%s' % (c1.index, c2.index, r.taa_ref)

    return r

def annotate_region_cdna(args, q, tpts, db):

    found = False
    for t in tpts:
        try:
            r = _annotate_region_cdna(args, q, t, db)
        except IncompatibleTranscriptError:
            continue
        except UnknownChromosomeError:
            continue
        found = True
        r.format(q.op)

    if not found:
        r = Record()
        r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))
        r.format(q.op)

    return


def annotate_region_protein(args, q, tpts, db):

    found = False
    for t in tpts:
        try:
            if q.tpt and t.name != q.tpt:
                raise IncompatibleTranscriptError('Transcript name unmatched')
            t.ensure_seq()

            r = Record()
            r.chrm = t.chrm
            r.tname = t.name

            if q.beg*3 > len(t) or q.end*3 > len(t):
                raise IncompatibleTranscriptError('codon nonexistent')

            tnuc_beg = q.beg*3 - 2
            tnuc_end = q.end*3
            r.gnuc_beg, r.gnuc_end = t.tnuc_range2gnuc_range(tnuc_beg, tnuc_end)
            r.natrefseq = t.seq[tnuc_beg-1:tnuc_end]
            r.refrefseq = reverse_complement(r.natrefseq) if t.strand == '-' else r.natrefseq
            taa_natrefseq = translate_seq(r.natrefseq)
            if q.beg_aa and q.beg_aa != taa_natrefseq[0]:
                raise IncompatibleTranscriptError('beginning reference amino acid unmatched')
            if q.end_aa and q.end_aa != taa_natrefseq[-1]:
                raise IncompatibleTranscriptError('ending reference amino acid unmatched')
            if q.refseq and not re.match(q.refseq.replace('x','[A-Z]'), taa_natrefseq): # != q.refseq:
                raise IncompatibleTranscriptError('reference sequence unmatched')
            r.tnuc_range = '%d_%d' % (tnuc_beg, tnuc_end)
            r.gnuc_range = '%d_%d' % (r.gnuc_beg, r.gnuc_end)
            r.taa_range = '%s%d_%s%d' % (taa_natrefseq[0], q.beg, taa_natrefseq[-1], q.end) if q.beg != q.end else '%d%s' % (q.beg, taa_natrefseq[0])
            r.pos = '%d-%d' % (r.gnuc_beg, r.gnuc_end)
            r.gene = t.gene.name
            r.strand = t.strand

            r.reg = RegCDSAnno(t)
            r.reg.from_taa_range(q.beg, q.end)
            r.append_info('protein_sequence=%s;cDNA_sequence=%s;gDNA_sequence=%s' % (printseq(taa_natrefseq), printseq(r.natrefseq), printseq(r.refrefseq)))
        except IncompatibleTranscriptError:
            continue
        except UnknownChromosomeError:
            continue
        found = True
        r.format(q.op)

    if not found:
        r = Record()
        r.append_info('no_valid_transcript_found_(from_%s_candidates)' % len(tpts))

        r.format(q.op)

    return


def annotate_region_gdna(args, q, db):
    
    normalize_reg(q)
    for reg in describe(args, q, db):
        r = Record()
        r.reg = reg
        r.chrm = q.tok

        r.set_promoter()

        if isinstance(reg, RegAnno):

            if hasattr(reg, 'intergenic'):

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
                
            elif hasattr(reg, 't'):
                c, p = reg.t.gpos2codon(q.pos)

                r.append_info("is_gene_body")
                r.tname = reg.t.format()
                r.gene = reg.t.gene.name
                r.strand = reg.t.strand

                if p.tpos == 0 and reg.t.transcript_type == 'protein_coding':
                    if c.seq in standard_codon_table:
                        r.taa_ref = standard_codon_table[c.seq]
                    r.taa_pos = c.index

                r.gnuc_pos = q.pos
                r.pos = q.pos
                r.gnuc_ref = faidx.refgenome.fetch_sequence(q.tok, q.pos, q.pos)
                r.tnuc_pos = p
                r.tnuc_ref = r.gnuc_ref if c.strand == '+' else complement(r.gnuc_ref)
                r.append_info('codon_pos=%s' % ('-'.join(map(str, c.locs)),))
                
            else:
                raise Exception() # shouldn't reach here

        elif isinstance(reg, RegSpanAnno):

            r.gnuc_range = '%d_%d' % (q.beg, q.end)
            r.pos = '%d-%d' % (q.beg, q.end)
            r.set_splice('included')
            
            if hasattr(reg, 't'):
                r.tname = reg.t.format()
            else:
                if hasattr(reg.b1, 't') and hasattr(reg.b2, 't'):
                    if reg.b1.t == reg.b2.t:
                        r.tname = reg.b1.t.format()
                    else:
                        r.tname = '%s,%s' % (reg.b1.t.format(), reg.b2.t.format())
                elif hasattr(reg.b1, 't'):
                    r.tname = reg.b1.t.format()+','
                elif hasattr(reg.b2, 't'):
                    r.tname = ','+reg.b2.t.format()

            if hasattr(reg, 't'):
                r.gene = reg.t.gene.name if reg.t.gene.name else '.'
                r.strand = reg.t.strand
                c1, p1 = reg.t.gpos2codon(q.beg)
                c2, p2 = reg.t.gpos2codon(q.end)

                if reg.t.strand == '+':
                    r.tnuc_range = '%s_%s' % (p1,p2)
                else:
                    r.tnuc_range = '%s_%s' % (p2,p1)

                # if protein coding transcript and not purely intronic region, set amino acid
                if reg.t.transcript_type == 'protein_coding' and not same_intron(p1, p2):
                    c1, p1 = reg.t.intronic_lean(p1, 'g_greater')
                    c2, p2 = reg.t.intronic_lean(p2, 'g_smaller')

                    if len(c1.seq) != 3 or len(c2.seq) != 3:
                        r.append_info("truncated_refseq_at_boundary_(start_codon_seq_%s_and_end_codon_seq_%s)" % (c1.seq, c2.seq))
                    elif c1.index == c2.index:
                        r.taa_ref = c1.aa()
                        r.taa_pos = c1.index
                        r.append_info('codon=%s' % c1.locformat())
                    else:
                        if reg.t.strand == '+':
                            r.taa_range = '%s%d_%s%d' % (c1.aa(), c1.index, c2.aa(), c2.index)
                        else:
                            r.taa_range = '%s%d_%s%d' % (c2.aa(), c2.index, c1.aa(), c1.index)
                        r.append_info('start_codon=%s;end_codon=%s' % (c1.locformat(), c2.locformat()))

        else:
            raise Exception()   # shouldn't reach

        r.format(q.op)
                


