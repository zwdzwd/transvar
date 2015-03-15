from transcripts import *
from record import *
import copy
import locale
locale.setlocale(locale.LC_ALL, '')

def _annotate_reg_gene_point(args, q, t):
    
    c, p, reg = t.gpos2codon(q.pos)
    r = Record()
    r.chrm = t.chrm
    r.tname = t.name
    r.reg = '%s (%s, %s)' % (t.gene.name, t.strand, reg.format())

    # at the ends of retained intron transcripts from ENSEMBL,
    # codon sequence is not always of length 3
    if p.tpos == 0:
        if c.seq in standard_codon_table:
            r.taa_ref = standard_codon_table[c.seq]
        r.taa_pos = c.index

    r.gnuc_pos = q.pos
    r.pos = q.pos
    r.gnuc_ref = faidx.refgenome.fetch_sequence(q.tok, q.pos, q.pos)
    r.tnuc_pos = p
    r.tnuc_ref = r.gnuc_ref if c.strand == '+' else complement(r.gnuc_ref)
    r.info = 'CodonPos=%s' % ('-'.join(map(str, c.locs)),)
    return r

def _annotate_reg_single_gene(args, q, t):

    r = Record()
    r.chrm = t.chrm
    r.tname = t.name
    # TODO: t.overlap_region change to RegSpanAnno framework
    # r.reg = '%s (%s, %s)' % (t.gene.name, t.strand, t.overlap_region(q.beg, q.end))
    r._reg_ = t.describe_span(q.beg, q.end)
    r.reg = '%s (%s, %s)' % (t.gene.name, t.strand, r._reg_.format())
    r.pos = '%d-%d' % (q.beg, q.end)

    r.gnuc_range = '%d_%d' % (q.beg, q.end)
    r.cbeg, r.pbeg, r.regbeg = t.gpos2codon(q.beg)
    r.cend, r.pend, r.regend = t.gpos2codon(q.end)
    if t.strand == '+':
        r.tnuc_range = '%s_%s' % (r.pbeg, r.pend)
    else:
        r.tnuc_range = '%s_%s' % (r.pend, r.pbeg)

    # detect cross boundary case
    if q.beg < t.cds_beg and q.end >= t.cds_beg:
        r.append_info('%s_loss' % ('start' if t.strand == '+' else 'stop'))
    if q.end > t.cds_end and q.beg <= t.cds_end:
        r.append_info('%s_loss' % ('stop' if t.strand == '+' else 'start'))

    if r.pbeg.tpos == 0 and r.pend.tpos == 0:
        if t.strand == '+':
            if not same_intron(r.pbeg, r.pend):
                if r.cbeg.index == r.cend.index:
                    r.taa_ref = r.cbeg.aa()
                    r.taa_pos = r.cbeg.index
                else:
                    r.taa_range = '%s%d_%s%d' % (r.cbeg.aa(), r.cbeg.index,
                                                 r.cend.aa(), r.cend.index)
        else:
            if not same_intron(r.pbeg, r.pend):
                if r.cbeg.index == r.cend.index:
                    r.taa_ref = r.cbeg.aa()
                    r.taa_pos = r.cbeg.index
                else:
                    r.taa_range = '%s%d_%s%d' % (r.cend.aa(), r.cend.index,
                                                 r.cbeg.aa(), r.cbeg.index)
        r.append_info('BEGCodon=%s;ENDCodon=%s' % (
            '-'.join(map(str, r.cbeg.locs)), '-'.join(map(str, r.cend.locs))))

    return r

def _annotate_reg_gene_long_range(args, q, tpts, genes, db):

    r = Record()
    r.chrm = tpts[0].chrm
    r.pos = '%d-%d' % (q.beg, q.end)
    r.reg = '%s bp covering %d genes' % (
        locale.format("%d", q.end-q.beg+1, grouping=True), len(genes))
    if len(genes) <= 5:
        r.reg += (' (%s)' % ';'.join([g.name for g in genes]))
    r.gnuc_range = '%d_%d' % (q.beg, q.end)

    qbeg = copy.copy(q)
    qbeg.end = qbeg.beg
    qend = copy.copy(q)
    qend.beg = qend.end
    for rbeg in __annotate_reg(args, qbeg, db):
        for rend in __annotate_reg(args, qend, db):
            r.tname = 'BEG=%s,END=%s' % (rbeg.tname, rend.tname)
            infocols = []
            infocols.append('BEGreg=%s' % rbeg.reg)
            infocols.append('BEGid=%s' % rbeg.format_id())
            infocols.append('ENDreg=%s' % rend.reg)
            infocols.append('ENDid=%s' % rend.format_id())
            r.info = ';'.join(infocols)
            yield r

def _annotate_reg_gene(args, q, db):

    tpts = [t for t in db.get_transcripts(q.tok, q.beg, q.end)]
    if tpts:
        if q.beg == q.end:
            if args.longest:
                tpts.sort(key=lambda t: len(t), reverse=True)
                tpts = tpts[:1]

            q.pos = q.beg
            for t in tpts:
                yield _annotate_reg_gene_point(args, q, t)
        else:
            genes = list(set([t.gene for t in tpts]))
            if len(genes) == 1:
                if args.longest:
                    tpts.sort(key=lambda t: len(t), reverse=True)
                    tpts = tpts[:1]

                for t in tpts:
                    yield _annotate_reg_single_gene(args, q, t)
            else:
                for r in _annotate_reg_gene_long_range(args, q, tpts, genes, db):
                    yield r

def __annotate_reg_intergenic(args, db, tok, beg, end):

    """ the function is also used in anno_snv.py """

    # annotate noncoding
    r = Record()
    r.chrm = tok
    tu, td  = db.get_closest_transcripts(tok, beg, end)
    if tu:
        up = 'up: %s bp to %s' % (
            locale.format('%d', beg - tu.end, grouping=True), tu.gene.name)
    else:
        up = 'up: %s bp to 5-telomere' % (
            locale.format('%d', beg, grouping=True), )
    if td:
        down = 'down: %s bp to %s' % (
            locale.format('%d', td.beg - end, grouping=True), td.gene.name)
    else:
        down = 'down: %s bp to 3-telomere' % (
            locale.format('%d', reflen(tok)-end, grouping=True), )
    r.reg = 'Intergenic (%s, %s)' % (up, down)
    if beg == end:
        r.gnuc_pos = beg
        r.pos = beg
    else:
        r.gnuc_range = '%d_%d' % (beg, end)
        r.pos = '%d-%d' % (beg, end)

    # # annotate extra noncoding features
    if 'GENCODE' in args.ffhs:
        iis = set()
        for entry in args.ffhs['GENCODE'].fetch(tok, beg, end+1):
            fields = entry.strip().split('\t')
            info = dict(re.findall(r'\s*([^"]*) "([^"]*)";', fields[8]))
            if 'gene_type' in info:
                ii = info['gene_type']
                if 'gene_name' in info: ii += '(%s)' % info['gene_name']
                iis.add(ii)
        r.info = 'gene_type=%s;' % ','.join(list(iis))

    return r
    
                    
def __annotate_reg(args, q, db):

    # check if location in a gene
    gene_found = False
    for r in _annotate_reg_gene(args, q, db):
        yield r
        gene_found = True

    if not gene_found:
        yield __annotate_reg_intergenic(args, db, q.tok, q.beg, q.end)

def _annotate_reg(args, q, db):
    normalize_reg(q)
    for r in __annotate_reg(args, q, db):
        r.format(q.op)
