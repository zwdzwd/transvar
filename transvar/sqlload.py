from sqlmodel import *
from utils import normalize_chrm
from parser import *
Session = sessionmaker(bind=engine, autoflush=False)
session = Session()

Base.metadata.create_all(engine)
def ensure_elem(eclass, ename):
    elems = session.query(eclass).filter_by(name=ename).all()
    if elems:
        return elems[0]
    else:
        elem = eclass(name=ename)
        session.add(elem)
        session.commit()
        return elem

def load_elem(rv_s, src_s, fn, func):

    rv = ensure_elem(DRefVersion, rv_s)
    src = ensure_elem(DSource, src_s)
    
    name2gene = {}
    func(fn, name2gene)
    for name, gene in name2gene.iteritems():
        gns = session.query(DGene).filter_by(name=name).all()   # query and find
        if gns:
            g = gns[0]
        else:
            g = DGene(name=name)
            session.add(g)
            session.flush()

        for t in gene.tpts:
            chrm = ensure_elem(DChromosome, normalize_chrm(t.chrm))
            ttype = ensure_elem(DTranscriptType, t.transcript_type)

            # infer cds begin and cds end
            if t.cds:
                t.cds.sort()
                cds_beg = t.cds[0][0]
                cds_end = t.cds[-1][1]
            elif hasattr(t, 'cds_beg'):
                cds_beg = t.cds_beg
                cds_end = t.cds_end
            else:
                cds_beg = -1
                cds_end = -1

            dt = DTranscript(
                chrm_id=chrm.id,
                transcript_type_id=ttype.id,
                beg = t.beg,
                end = t.end,
                source_id = src.id,
                refversion_id = rv.id,
                name = t.name,
                cds_beg = cds_beg,
                cds_end = cds_end,
                gene_id = g.id,
                strand = 1 if t.strand == '-' else 0,
            )
            session.add(dt)
            session.flush()
            for exbeg, exend in t.exons:
                e = DExon(tid=dt.id, beg=exbeg, end=exend)
                session.add(e)

    session.commit()

if __name__ == '__main__':

    # load_elem('hg19', 'UCSC', 'transvar.download/hg19.ucsc.txt.gz', parse_ucsc_refgene)
    # load_elem('hg19', 'Ensembl', 'transvar.download/hg19.ensembl.gtf.gz', parse_ensembl_gtf)
    # load_elem('hg19', 'GENCODE', 'transvar.download/hg19.gencode.gtf.gz', parse_gencode_gtf)
    # load_elem('hg19', 'CCDS', 'transvar.download/hg19.ccds.txt', parse_ccds_table)
    # load_elem('hg19', 'RefSeq', 'transvar.download/hg19.refseq.gff.gz', parse_refseq_gff)
    # # load_elem('hg19', 'AceView', 'transvar.download/hg19.aceview.gff.gz', parse_aceview_transcripts)
    # load_elem('hg38', 'Ensembl', 'transvar.download/hg38.ensembl.gtf.gz', parse_ensembl_gtf)
    # load_elem('hg38', 'UCSC', 'transvar.download/hg38.ucsc.txt.gz', parse_ucsc_refgene)
    # load_elem('hg38', 'GENCODE', 'transvar.download/hg38.gencode.gtf.gz', parse_gencode_gtf)
    # load_elem('hg38', 'RefSeq', 'transvar.download/hg38.refseq.gff.gz', parse_refseq_gff)
    # load_elem('hg18', 'UCSC', 'transvar.download/hg18.ucsc.txt.gz', parse_ucsc_refgene)
    # load_elem('hg18', 'GENCODE', 'transvar.download/hg18.gencode.gtf.gz', parse_gencode_gtf)
    # # load_elem('hg18', 'Ensembl', 'transvar.download/hg18.ensembl.gtf.gz', parse_ensembl_gtf)
    # load_elem('hg18', 'CCDS', 'transvar.download/hg18.ccds.txt', parse_ccds_table)
    load_elem('hg18', 'RefSeq', 'transvar.download/hg18.refseq.gff.gz', parse_refseq_gff)
    # # load_elem('hg18', 'AceView', 'transvar.download/hg18.aceview.gff.gz', parse_aceview_transcripts)
    load_elem('mm10', 'Ensembl' ,'transvar.download/mm10.ensembl.gtf.gz', parse_ensembl_gtf)
    load_elem('mm10', 'RefSeq' ,'transvar.download/mm10.refseq.gff.gz', parse_refseq_gff)
    load_elem('mm10', 'CCDS', 'transvar.download/mm10.ccds.txt', parse_ccds_table)
