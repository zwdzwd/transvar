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

def load_hg19_gencode():
    
    rv = ensure_elem(RefVersion, 'hg19')
    src = ensure_elem(Source, 'GENCODE')
    ft_cds = ensure_elem(FeatureType, 'protein_coding')
    
    name2gene = {}
    trs.parse_gencode_gtf('transvar.download/hg19.gencode.gtf.gz', name2gene)
    for name, gene in name2gene.iteritems():
        gns = session.query(Gene).filter_by(name=name).all()   # query and find
        if gns:
            g = gns[0]
        else:
            g = Gene(name=name)
            session.add(g)
        for transcript in gene.tpts:
            chms = session.query(Chromosome).filter_by(name=transcript.chrm).all()
            if chms:
                chm = chms[0]
            else:
                chm = Chromosome(name=transcript.chrm)
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

            if not transcript.cds: continue
            transcript.cds.sort()
            t = Transcript(
                id = f.id,
                name = transcript.name,
                cds_beg = transcript.cds[0][0],
                cds_end = transcript.cds[-1][1],
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
    load_hg19_gencode()
