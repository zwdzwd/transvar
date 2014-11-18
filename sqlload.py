from sqlmodel import *
import transcripts as trs
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

    rv = ensure_elem(RefVersion, rv_s)
    src = ensure_elem(Source, src_s)
    ft_cds = ensure_elem(FeatureType, 'protein_coding')
    
    name2gene = {}
    func(fn, name2gene)
    for name, gene in name2gene.iteritems():
        gns = session.query(Gene).filter_by(name=name).all()   # query and find
        if gns:
            g = gns[0]
        else:
            g = Gene(name=name)
            session.add(g)
        for transcript in gene.tpts:
            chrm = normalize_chrm(transcript.chrm)
            chms = session.query(Chromosome).filter_by(name=chrm).all()
            if chms:
                chm = chms[0]
            else:
                chm = Chromosome(name=chrm)
                session.add(chm)
                session.flush()

            f = Feature(ftype=ft_cds.id,
                        chrm_id=chm.id,
                        beg=transcript.beg,
                        end=transcript.end,
                        source_id=src.id,
                        refversion_id=rv.id)
            session.add(f)
            session.flush()

            if transcript.cds:
                transcript.cds.sort()
                cds_beg = transcript.cds[0][0]
                cds_end = transcript.cds[-1][1]
            elif hasattr(transcript, 'cds_beg'):
                cds_beg = transcript.cds_beg
                cds_end = transcript.cds_end
            else:
                continue

            t = Transcript(
                id = f.id,
                name = transcript.name,
                cds_beg = cds_beg,
                cds_end = cds_end,
                gene_id = g.id,
                strand = 1 if transcript.strand == '-' else 0,
            )
            session.add(t)
            session.flush()
            for exbeg, exend in transcript.exons:
                e = Exon(tid=t.id, beg=exbeg, end=exend)
                session.add(e)

    session.commit()

if __name__ == '__main__':
    load_elem('hg19', 'UCSC', 'transvar.download/hg19.ucsc.txt.gz', trs.parse_ucsc_refgene)
    load_elem('hg19', 'GENCODE', 'transvar.download/hg19.gencode.gtf.gz', trs.parse_gencode_gtf)
    load_elem('hg19', 'Ensembl', 'transvar.download/hg19.ensembl.gtf.gz', trs.parse_ensembl_gtf)
    load_elem('hg19', 'CCDS', 'transvar.download/hg19.ccds.txt', trs.parse_ccds_table)
    load_elem('hg19', 'RefSeq', 'transvar.download/hg19.refseq.gff.gz', trs.parse_refseq_gff)
    load_elem('hg19', 'AceView', 'transvar.download/hg19.aceview.gff.gz', trs.parse_aceview_transcripts)
    load_elem('hg38', 'Ensembl', 'transvar.download/hg38.ensembl.gtf.gz', trs.parse_ensembl_gtf)
    load_elem('hg38', 'UCSC', 'transvar.download/hg38.ucsc.txt.gz', trs.parse_ucsc_refgene)
    load_elem('hg38', 'GENCODE', 'transvar.download/hg38.gencode.gtf.gz', trs.parse_gencode_gtf)
    load_elem('hg38', 'RefSeq', 'transvar.download/hg38.refseq.gff.gz', trs.parse_refseq_gff)
    load_elem('hg18', 'UCSC', 'transvar.download/hg18.ucsc.txt.gz', trs.parse_ucsc_refgene)
    load_elem('hg18', 'GENCODE', 'transvar.download/hg18.gencode.gtf.gz', trs.parse_gencode_gtf)
    # load_elem('hg18', 'Ensembl', 'transvar.download/hg18.ensembl.gtf.gz', trs.parse_ensembl_gtf)
    load_elem('hg18', 'CCDS', 'transvar.download/hg18.ccds.txt', trs.parse_ccds_table)
    load_elem('hg18', 'RefSeq', 'transvar.download/hg18.refseq.gff.gz', trs.parse_refseq_gff)
    load_elem('hg18', 'AceView', 'transvar.download/hg18.aceview.gff.gz', trs.parse_aceview_transcripts)
