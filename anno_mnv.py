from transcripts import *
from record import *
from anno_reg import __annotate_reg_intergenic, _annotate_reg_single_gene

# def __annotate_mnv_gene(args, q, t):

#     c1, p1, reg1 = t.gpos2codon(q.beg)
#     c2, p2, reg2 = t.gpos2codon(q.end)

#     reg1f = reg1.format()
#     reg2f = reg2.format()
#     if reg1f == reg2f:
#         regf = reg1f
#     else:
#         regf = '%s - %s'

#     r = Record()
#     r.chrm = t.chrm
#     r.tname = t.name
#     r.reg = '%s (%s, %s)' % (t.gene.name, t.strand, reg.format())
#     r.pos = q.pos

#     # at the ends of retained intron transcripts from ENSEMBL,
#     # codon sequence is not always of length 3
#     if p.tpos == 0:
#         r.taa_ref = codon2aa(c.seq)
#         r.taa_pos = c.index

#         if q.alt:
#             if c.strand == "+":
#                 alt_seq = set_seq(c.seq, c.locs.index(q.pos), q.alt)
#             else:
#                 alt_seq = set_seq(c.seq, 2-c.locs.index(q.pos), complement(q.alt))
#             r.taa_alt = codon2aa(alt_seq)

#     r.gnuc_pos = q.pos
#     # r.gnuc_ref = c.refseq()[c.locs.index(q.pos)]
#     r.gnuc_ref = q.ref
#     r.gnuc_alt = q.alt
#     r.tnuc_pos = p
#     if c.strand == '+':
#         r.tnuc_ref = r.gnuc_ref
#         r.tnuc_alt = r.gnuc_alt
#     else:
#         r.tnuc_ref = complement(r.gnuc_ref)
#         r.tnuc_alt = complement(r.gnuc_alt) if r.gnuc_alt else ''

#     r.append_info('CodonPos=%s' % '-'.join(map(str, c.locs)))
#     r.append_info('NCodonSeq=%s' % c.seq)
#     if hasattr(r, 'taa_ref') and hasattr(r, 'taa_alt'):
#         r.append_info('Syn' if r.taa_ref == r.taa_alt else 'Nonsyn')

#     return r

def add_mnv_single_gene(r, q, t):

    """ add reference and alternative annotation
    r.cbeg, r.pbeg, r.regbeg
    r.cend, r.pend, r.regend are all preset in _annotate_reg_single_gene
    """
    r.gnuc_refseq = q.refseq
    r.gnuc_altseq = q.altseq
    r.gnuc_range = '%d_%d%s>%s' % (q.beg, q.end, r.gnuc_refseq, r.gnuc_altseq)
    if t.strand == '+':
        r.tnuc_refseq = r.gnuc_refseq
        r.tnuc_altseq = r.gnuc_altseq
        r.tnuc_beg = r.pbeg
        r.tnuc_end = r.pend
    else:
        r.tnuc_refseq = reverse_complement(r.gnuc_refseq)
        r.tnuc_altseq = reverse_complement(r.gnuc_altseq)
        r.tnuc_beg = r.pend
        r.tnuc_end = r.pbeg
    r.tnuc_range = '%s_%s%s>%s' % (r.tnuc_beg, r.tnuc_end,
                                   r.tnuc_refseq, r.tnuc_altseq)


    # r._reg_ is the raw RegSpanAnno object
    # predict protein change if entire region is in exon
    if r._reg_.in_exon():
        r.taa_range = t.tnuc_mnv_coding(r.tnuc_beg.pos, r.tnuc_end.pos, r.tnuc_altseq, r)
    return

def _annotate_mnv_gene(args, q, db):

    tpts = [t for t in db.get_transcripts(q.tok, q.beg, q.end)]
    genes = list(set([t.gene for t in tpts]))
    if len(genes) == 1:         # TODO should judge by the number of genes, there might be nested genes.
        for t in tpts:
            r = _annotate_reg_single_gene(args, q, t)
            add_mnv_single_gene(r, q, t)
            yield r

def _annotate_mnv(args, q, db):

    # check reference sequence
    gnuc_refseq = faidx.refgenome.fetch_sequence(q.tok, q.beg, q.end)
    if q.refseq and gnuc_refseq != q.refseq:
        r = Record()
        r.chrm = q.tok
        r.pos = '%d-%d (block substitution)' % (q.beg, q.end)
        r.info = "InvalidReference:%s(%s)" % (q.refseq, gnuc_refseq)
        r.format(q.op)
        err_print("Invalid reference base %s (%s), maybe wrong reference?" % (q.refseq, gnuc_refseq))
        return
    else:                       # make sure q.refseq exists
        q.refseq = gnuc_refseq

    # check if location in a gene
    gene_found = False
    for r in _annotate_mnv_gene(args, q, db):
        r.format(q.op)
        gene_found = True

    if not gene_found:
        r = __annotate_reg_intergenic(args, db, q.tok, q.beg, q.end)
        # r.gnuc_range = '%d_%d%s>%s'
        if hasattr(r, "gnuc_range"):
            r.gnuc_range += '%s>%s' % (q.refseq, q.altseq)
        # TODO : add to reflect reference and alternative
        r.format(q.op)
