from transcripts import *
from record import *
from anno_reg import __annotate_reg_intergenic
import locale
locale.setlocale(locale.LC_ALL, 'en_US')

def __annotate_snv_gene(args, q, t):

    c, p, reg = t.gpos2codon(q.tok, q.pos)

    r = Record()
    r.chrm = t.chrm
    r.tname = t.name
    r.reg = '%s (%s, %s)' % (t.gene.name, t.strand, reg)
    r.pos = q.pos

    # at the ends of retained intron transcripts from ENSEMBL,
    # codon sequence is not always of length 3
    if p.tpos == 0:
        r.taa_ref = codon2aa(c.seq)
        r.taa_pos = c.index

        if q.alt:
            if c.strand == "+":
                alt_seq = set_seq(c.seq, c.locs.index(q.pos), q.alt)
            else:
                alt_seq = set_seq(c.seq, 2-c.locs.index(q.pos), complement(q.alt))
            r.taa_alt = codon2aa(alt_seq)

    r.gnuc_pos = q.pos
    # r.gnuc_ref = c.refseq()[c.locs.index(q.pos)]
    r.gnuc_ref = faidx.refgenome.fetch_sequence(q.tok, q.pos, q.pos)
    r.gnuc_alt = q.alt
    r.tnuc_pos = p
    if c.strand == '+':
        r.tnuc_ref = r.gnuc_ref
        r.tnuc_alt = r.gnuc_alt
    else:
        r.tnuc_ref = complement(r.gnuc_ref)
        r.tnuc_alt = complement(r.gnuc_alt) if r.gnuc_alt else ''

    r.info = 'CodonPos=%s;NCodonSeq=%s' % ('-'.join(map(str, c.locs)), c.seq)

    return r

def _annotate_snv_gene(args, q, db):

    tpts = [t for t in db.get_transcripts(q.tok, q.pos, q.pos)]
    if tpts:
        if args.longest:
            tpts.sort(key=lambda t: len(t), reverse=True)
            tpts = tpts[:1]

        for t in tpts:
            try:
                r = __annotate_snv_gene(args, q, t)
            except IncompatibleTranscriptError:
                continue
            yield r

    # if args.longest:
    #     tc_iter = gpos2codon_longest(thash, q.tok, q.pos)
    # else:
    #     tc_iter = gpos2codon(thash, q.tok, q.pos)

    # found = False
    # for t, c in tc_iter:
    #     if isinstance(c, Codon):
    #         found = True


def _annotate_snv(args, q, db):
    # check if location in a gene
    gene_found = False
    for r in _annotate_snv_gene(args, q, db):
        r.format(q.op)
        gene_found = True

    if not gene_found:
        r = __annotate_reg_intergenic(args, db, q.tok, q.pos, q.pos)
        r.format(q.op)
        # # annotate noncoding
        # r = Record()
        # tu, td  = thash.get_closest_transcripts(q.tok, q.pos, q.pos)
        # if tu:
        #     up = 'up: %s bp to %s' % (
        #         locale.format('%d', q.pos - tu.end, grouping=True), tu.gene.name)
        # else:
        #     up = 'up: %s bp to 5-telomere' % (
        #         locale.format('%d', q.pos, grouping=True), )
        # if td:
        #     down = 'down: %s bp to %s' % (
        #         locale.format('%d', td.beg - q.pos, grouping=True), td.gene.name)
        # else:
        #     down = 'down: %s bp to 3-telomere' % (
        #         locale.format('%d', reflen(q.tok)-q.pos, grouping=True), )
        # r.reg = 'Intergenic (%s, %s)' % (up, down)
        # r.format(q.op)

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
