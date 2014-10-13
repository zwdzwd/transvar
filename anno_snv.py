from transcripts import *
from record import *

def _annotate_snv_gene(args, q, thash):

    if args.longest:
        tc_iter = gpos2codon_longest(thash, q.tok, q.pos)
    else:
        tc_iter = gpos2codon(thash, q.tok, q.pos)

    found = False
    for t, c in tc_iter:
        if isinstance(c, Codon):
            found = True

            r = Record()
            r.chrm = t.chrm
            r.tname = t.name
            r.reg = '%s (%s, coding)' % (t.gene.name, t.strand)
            r.pos = '-'.join(map(str, c.locs))

            # at the ends of retained intron transcripts from ENSEMBL,
            # codon sequence is not always of length 3
            if c.seq in standard_codon_table:
                r.taa_ref = standard_codon_table[c.seq]
            r.taa_pos = c.index
            if q.alt:
                if c.strand == "+":
                    alt_seq = set_seq(c.seq, c.locs.index(q.pos), q.alt)
                else:
                    alt_seq = set_seq(c.seq, 2-c.locs.index(q.pos), complement(q.alt))
                r.taa_alt = standard_codon_table[alt_seq]

            r.gnuc_pos = q.pos
            r.gnuc_ref = c.refseq()[c.locs.index(q.pos)]
            r.gnuc_alt = q.alt

            if c.strand == '+':
                r.tnuc_ref = r.gnuc_ref
                r.tnuc_alt = r.gnuc_alt
                r.tnuc_pos = (c.index-1)*3 + c.locs.index(q.pos) + 1
            else:
                r.tnuc_ref = complement(r.gnuc_ref)
                r.tnuc_alt = complement(r.gnuc_alt) if r.gnuc_alt else ''
                r.tnuc_pos = c.index*3 - c.locs.index(q.pos)

            yield r


def _annotate_snv(args, q, thash):

    # check if location in a gene
    gene_found = False
    for r in _annotate_snv_gene(args, q, thash):
        r.format(q.op)
        gene_found = True

    if not gene_found:
        # annotate noncoding
        pass

    #     elif isinstance(c, NonCoding):
    #         found = True

    #         r = Record()
    #         r.chrm = t.chrm
    #         r.gnuc_pos = q.pos
    #         r.tname = t.name
    #         r.reg = '%s (%s noncoding)' % (t.gene.name, t.strand)
    #         r.info = c.format()
    #         r.format(q.op)

    # if not found:
    #     r = Record()
    #     r.gnuc_ref = q.ref
    #     r.gnuc_alt = q.alt
    #     r.gnuc_pos = q.pos
    #     r.info = 'status=NoValidTranscriptFound'
    #     r.format(q.op)
