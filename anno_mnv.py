from transcripts import *
from record import *
# from anno_reg import __annotate_reg_intergenic, _annotate_reg_single_gene
from describe import *

def _annotate_mnv(args, q, db):

    # check reference sequence
    gnuc_refseq = faidx.refgenome.fetch_sequence(q.tok, q.beg, q.end)
    if q.refseq and gnuc_refseq != q.refseq:
        
        r = Record()
        r.chrm = q.tok
        r.pos = '%d-%d (block substitution)' % (q.beg, q.end)
        r.info = "invalid_reference_seq_%s_(expect_%s)" % (q.refseq, gnuc_refseq)
        r.format(q.op)
        err_print("Warning: %s invalid reference %s (expect %s), maybe wrong reference?" % (q.op, q.refseq, gnuc_refseq))
        return
    
    else:                       # make sure q.refseq exists
        q.refseq = gnuc_refseq

    for reg in describe(args, q, db):

        r = Record()
        r.reg = reg
        r.chrm = q.tok
        r.gnuc_refseq = q.refseq
        r.gnuc_altseq = q.altseq
        r.gnuc_range = '%d_%d%s>%s' % (q.beg, q.end, r.gnuc_refseq, r.gnuc_altseq)

        if hasattr(reg, 't'):

            t = reg.t
            r.tname = t.format()
            r.gene = t.gene.name
            r.strand = t.strand

            c1, p1 = t.gpos2codon(q.beg, intronic_policy="g_greater")
            c2, p2 = t.gpos2codon(q.end, intronic_policy="g_smaller")

            if t.strand == '+':
                tnuc_beg = p1
                tnuc_end = p2
                tnuc_refseq = q.refseq
                tnuc_altseq = q.altseq
            else:
                tnuc_beg = p2
                tnuc_end = p1
                tnuc_refseq = reverse_complement(q.refseq)
                tnuc_altseq = reverse_complement(q.altseq)
            r.tnuc_range = '%s_%s%s>%s' % (tnuc_beg, tnuc_end, tnuc_refseq, tnuc_altseq)

            if r.reg.entirely_in_cds():
                try:
                    t.tnuc_mnv_coding(tnuc_beg.pos, tnuc_end.pos, tnuc_altseq, r)
                except IncompatibleTranscriptError as inst:
                    _beg, _end, _seqlen = inst
                    r.append_info('mnv_(%s-%s)_at_truncated_refseq_of_length_%d' % (_beg, _end, _seqlen))

        elif isinstance(reg, RegSpanAnno):

            tnames = []
            strands = []
            genes = []
            if hasattr(reg.b1, 't'):
                if reg.b1.t.name not in tnames:
                    tnames.append(reg.b1.t.name)
                    strands.append(reg.b1.t.strand)
                    genes.append(reg.b1.t.gene.name)
                    
            if hasattr(reg.b2, 't'):
                if reg.b2.t.name not in tnames:
                    tnames.append(reg.b2.t.name)
                    strands.append(reg.b2.t.strand)
                    genes.append(reg.b2.t.gene.name)

            r.tname = ','.join(tnames)
            r.strand = ','.join(strands)
            r.gene = ','.join(genes)

        r.format(q.op)
