import sqlalchemy
from sqlalchemy import create_engine
engine=create_engine('mysql://transvar_user:transvar_user@localhost/transvar_anno', echo=True)
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()
from sqlalchemy import Column, Integer, String, ForeignKey, Boolean
from sqlalchemy.orm import relationship

class FeatureType(Base):
    __tablename__ = 'feature_type'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100))
    features = relationship('Feature', backref='type')

class Source(Base):
    __tablename__ = 'source'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100))
    features = relationship('Feature', backref='source')
    
class RefVersion(Base):
    __tablename__ = 'refversion'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100))
    features = relationship('Feature', backref='refversion')

class Chromosome(Base):
    __tablename__ = 'chromosome'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100))
    features = relationship('Feature', backref='chrm')

class Feature(Base):
    __tablename__ = 'feature'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key = True, autoincrement=True)
    ftype = Column(Integer, ForeignKey('feature_type.id'))
    chrm_id = Column(Integer, ForeignKey('chromosome.id'))
    beg = Column(Integer, index=True)
    end = Column(Integer, index=True)
    source_id = Column(Integer, ForeignKey('source.id'))    # primary key?
    refversion_id = Column(Integer, ForeignKey('refversion.id'))
    transcripts = relationship('Transcript', backref='feature')

class Gene(Base):
    __tablename__ = 'gene'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key = True, autoincrement=True)
    name = Column(String(100), primary_key=True)
    transcripts = relationship('Transcript', backref='gene')

class Transcript(Base):
    __tablename__ = 'transcript'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, ForeignKey('feature.id'), primary_key = True)
    cds_beg = Column(Integer)
    cds_end = Column(Integer)
    gene_id = Column(Integer, ForeignKey('gene.id'))
    version = Column(Integer, primary_key = True)
    strand = Column(Boolean)
    exons = relationship('Exon', backref='transcript')

class Exon(Base):
    __tablename__ = 'exon'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key=True)
    tid = Column(Integer, ForeignKey('transcript.id'))
    beg = Column(Integer),
    end = Column(Integer),

# class User(Base):
#     __tablename__ = 'users'
#     id = Column(Integer, primary_key=True)
#     name = Column(String(100))
#     fullname = Column(String(100))
#     password = Column(String(100))
#     def __repr__(self):
#         return "<User(name='%s', fullname='%s', password='%s')>" % (self.name, self.fullname, self.password)

Base.metadata.create_all(engine)
from sqlalchemy.orm import sessionmaker
Session = sessionmaker(bind=engine, autoflush=False)
session = Session()

rv = RefVersion(name='hg19')
session.add(rv)
src = Source(name='GENCODE')
session.add(src)
session.commit()

fts = session.query(FeatureType).filter_by(name='protein_coding').all()
if fts:
    ft_cds = fts[0]
else:
    ft_cds = FeatureType(name='protein_coding')
    session.add(ft_cds)
    session.commit()

fts = session.query(FeatureType).filter_by(name='exon').all()
if fts:
    ft_ex = fts[0]
else:
    ft_ex = FeatureType(name='exon')
    session.add(ft_ex)
    session.commit()
    
import transcripts as trs
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
            cds_beg=transcript.cds[0][0],
            cds_end=transcript.cds[-1][1],
            gene_id = g.id,
            strand = 1 if transcript.strand == '-' else 0,
        )
        session.add(t)
        session.flush()
        for exbeg, exend in transcript.exons:
            e = Exon(tid=t.id, beg=exbeg, end=exend)
            session.add(e)

session.commit()

# transcript = Transcript(cds_beg=100, cds_end=200)
