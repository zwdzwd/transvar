from transcripts import *
from record import *
import copy
from describe import *

def _annotate_reg(args, q, db):
    
    normalize_reg(q)
    for reg in describe(args, q, db):
        r = Record()
        r.reg = reg
        r.chrm = q.tok
        if hasattr(reg, 'promoter') and reg.promoter > 0:
            r.append_info('promoter_overlap_%d_bp(%1.2f%%)' % (reg.promoter, float(reg.promoter)/(q.end-q.beg)*100))

        if isinstance(reg, RegAnno):

            if hasattr(reg, 'intergenic'):

                r.gnuc_pos = q.pos if hasattr(q,'pos') else q.beg
                r.pos = r.gnus_pos

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
                c, p = reg.t.gpos2codon(q.pos, args)

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

            if hasattr(reg, 't'):
                r.tname = reg.t.format()
                r.gene = reg.t.gene.name if reg.t.gene.name else '.'
                r.strand = reg.t.strand
                c1, p1 = reg.t.gpos2codon(q.beg, args)
                c2, p2 = reg.t.gpos2codon(q.end, args)

                if reg.t.strand == '+':
                    r.tnuc_range = '%s_%s' % (p1,p2)
                else:
                    r.tnuc_range = '%s_%s' % (p2,p1)

                # if protein coding transcript and not purely intronic region, set amino acid
                if reg.t.transcript_type == 'protein_coding' and not same_intron(p1, p2):
                    c1, p1 = reg.t.intronic_lean(c1, p1, 'g_greater')
                    c2, p2 = reg.t.intronic_lean(c2, p2, 'g_smaller')

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
                


