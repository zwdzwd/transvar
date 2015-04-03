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
        err_print("invalid reference base %s (expect %s), maybe wrong reference?" % (q.refseq, gnuc_refseq))
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

            r.tname = reg.t.format()
            r.gene = reg.t.gene.name
            r.strand = reg.t.strand

            c1, p1 = reg.t.gpos2codon(q.beg, args, intronic_policy="g_greater")
            c2, p2 = reg.t.gpos2codon(q.end, args, intronic_policy="g_smaller")

            if t.strand == '+':
                k = (p1, p2, q.refseq, q.altseq)
            else:
                k = (p2, p1, reverse_complement(q.refseq), reverse_complement(q.altseq))
            r.tnuc_range = '%s_%s%s>%s' % k
            
            if r.reg.entire_in_cds():
                try:
                    r.taa_range = reg.t.tnuc_mnv_coding(r.tnuc_beg.pos, r.tnuc_end.pos, r.tnuc_altseq, r)
                except IncompatibleTranscriptError:
                    raise Exception()
            
        r.format(q.op)
