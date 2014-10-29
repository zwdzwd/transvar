import sqlalchemy
from sqlalchemy import create_engine
engine=create_engine('mysql://transvar_user:transvar_user@localhost/transvar_anno', echo=True)
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()
from sqlalchemy import Column, Integer, String, ForeignKey, Boolean
from sqlalchemy.orm import relationship, sessionmaker

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
    name = Column(String(30)),
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
    beg = Column(Integer)
    end = Column(Integer)

# class User(Base):
#     __tablename__ = 'users'
#     id = Column(Integer, primary_key=True)
#     name = Column(String(100))
#     fullname = Column(String(100))
#     password = Column(String(100))
#     def __repr__(self):
#         return "<User(name='%s', fullname='%s', password='%s')>" % (self.name, self.fullname, self.password)

# transcript = Transcript(cds_beg=100, cds_end=200)
