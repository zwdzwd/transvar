from transcripts import *
from record import *
import copy
import locale
locale.setlocale(locale.LC_ALL, 'en_US')


def _annotate_reg_gene_point(args, q, t):
    
    c, p, reg = t.gpos2codon(q.tok, q.pos)
    r = Record()
    r.chrm = t.chrm
    r.tname = t.name
    r.reg = '%s (%s, %s)' % (t.gene.name, t.strand, reg)

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

def _annotate_reg_gene_short_range(args, q, t):

    r = Record()
    r.chrm = t.chrm
    r.tname = t.name
    r.reg = '%s (%s, %s)' % (t.gene.name, t.strand, t.overlap_region(q.beg, q.end))
    r.pos = '%d-%d' % (q.beg, q.end)
    cbeg, pbeg, regbeg = t.gpos2codon(q.tok, q.beg)
    cend, pend, regend = t.gpos2codon(q.tok, q.end)
    r.gnuc_range = '%d_%d' % (q.beg, q.end)
    if t.strand == '+':
        r.tnuc_range = '%s_%s' % (pbeg, pend)
        if cbeg.index == cend.index:
            r.taa_ref = cbeg.aa()
            r.taa_pos = cbeg.index
        else:
            r.taa_range = '%s%d_%s%d' % (cbeg.aa(), cbeg.index, cend.aa(), cend.index)
    else:
        r.tnuc_range = '%s_%s' % (pend, pbeg)
        if cbeg.index == cend.index:
            r.taa_ref = cbeg.aa()
            r.taa_pos = cbeg.index
        else:
            r.taa_range = '%s%d_%s%d' % (cend.aa(), cend.index, cbeg.aa(), cbeg.index)
    r.info = 'BEGCodon=%s;ENDCodon=%s' % (
        '-'.join(map(str, cbeg.locs)), '-'.join(map(str, cend.locs)))
    
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
        if args.longest:
            tpts.sort(key=lambda t: len(t), reverse=True)
            tpts = tpts[:1]

        if q.beg == q.end:
            q.pos = q.beg
            for t in tpts:
                yield _annotate_reg_gene_point(args, q, t)
        else:
            genes = list(set([t.gene for t in tpts]))
            if len(genes) == 1:
                for t in tpts:
                    yield _annotate_reg_gene_short_range(args, q, t)
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
