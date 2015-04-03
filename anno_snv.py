from transcripts import *
from record import *
# from anno_reg import __annotate_reg_intergenic
from describe import *
# import locale
# locale.setlocale(locale.LC_ALL, '')

def _annotate_snv(args, q, db):

    # check reference base
    gnuc_ref = faidx.refgenome.fetch_sequence(q.tok, q.pos, q.pos)
    if q.ref and gnuc_ref != q.ref:
        
        r = Record()
        r.chrm = q.tok
        r.pos = q.pos
        r.info = "invalid_reference_base_%s_(expect_%s)" % (q.ref, gnuc_ref)
        r.format(q.op)
        err_print("invalid reference base %s (expect %s), maybe wrong reference?" % (q.ref, gnuc_ref))
        return
    
    else:
        q.ref = gnuc_ref

    for reg in describe(args, q, db):

        r = Record()
        r.reg = reg
        r.chrm = q.tok
        r.gnuc_pos = q.pos
        r.pos = r.gnuc_pos
        r.gnuc_ref = gnuc_ref
        r.gnuc_alt = q.alt if q.alt else ''
        
        if hasattr(reg, 't'):

            c,p = reg.t.gpos2codon(q.pos, args)

            r.tname = reg.t.format()
            r.gene = reg.t.gene.name
            r.strand = reg.t.strand
            r.tnuc_pos = p

            if c.strand == '+':
                r.tnuc_ref = r.gnuc_ref
                r.tnuc_alt = r.gnuc_alt
            else:
                r.tnuc_ref = complement(r.gnuc_ref)
                r.tnuc_alt = complement(r.gnuc_alt) if r.gnuc_alt else ''

            if p.tpos == 0:
                if c.seq in standard_codon_table:
                    if c.seq in standard_codon_table:
                        r.taa_ref = standard_codon_table[c.seq]
                    r.taa_pos = c.index

                if q.alt:
                    if c.strand == '+':
                        alt_seq = set_seq(c.seq, c.locs.index(q.pos), q.alt)
                    else:
                        alt_seq = set_seq(c.seq, 2-c.locs.index(q.pos), complement(q.alt))
                    r.taa_alt = codon2aa(alt_seq)
                    if r.taa_alt == r.taa_ref:
                        r.append_info('synonymous')
                    else:
                        r.append_info('nonsynonymous')

                r.append_info('codon_pos=%s' % (c.locformat(),))
                r.append_info('ref_codon_seq=%s' % c.seq)
                
        r.format(q.op)
