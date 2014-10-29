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
    
class Feature(Base):
    __tablename__ = 'feature'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key = True, autoincrement=True)
    ftype = Column(Integer, ForeignKey('feature_type.id'))
    beg = Column(Integer)
    end = Column(Integer)
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
    id = Column(Integer, ForeignKey('feature.id'), primary_key=True)
    tid = Column(Integer, ForeignKey('transcript.id'))
    exonno = Column(Integer)

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
Session = sessionmaker(bind=engine)
session = Session()

rv = RefVersion(name='hg19')
session.add(rv)
src = Source(name='GENCODE')
session.add(src)
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
        fts = session.query(FeatureType).filter_by(name='protein_coding').all()
        if fts:
            ft = fts[0]
        else:
            ft = FeatureType(name='protein_coding')
            session.add(ft)

        f = Feature(ftype=ft.id,
                    beg=transcript.beg,
                    end=transcript.end,
                    source_id=src.id,
                    refversion_id=rv.id)
        session.add(f)
    #     transcript = Transcript(cds_beg, cds_end, gene.id)
    #     session.add(transcript)
    #     for exbeg, exend in transcript.exons:
    #         exon = Exon(beg, end)
    #         session.add(exon)

session.commit()

# transcript = Transcript(cds_beg=100, cds_end=200)
