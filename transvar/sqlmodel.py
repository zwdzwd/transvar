import sqlalchemy
from sqlalchemy import create_engine
engine=create_engine('mysql://transvar_user:transvar_user@localhost/transvar_anno', echo=False)
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()
from sqlalchemy import Column, Integer, String, ForeignKey, Boolean
from sqlalchemy.orm import relationship, sessionmaker

# class FeatureType(Base):
#     __tablename__ = 'feature_type'
#     __table_args__ = {'mysql_engine':'InnoDB'}
#     id = Column(Integer, primary_key=True, autoincrement=True)
#     name = Column(String(100))
#     features = relationship('Feature', backref='type')

class DSource(Base):
    
    __tablename__ = 'source'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100))
    transcripts = relationship('DTranscript', backref='source')
    
class DRefVersion(Base):
    
    __tablename__ = 'refversion'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100))
    transcripts = relationship('DTranscript', backref='refversion')

class DChromosome(Base):

    __tablename__ = 'chromosome'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key=True, autoincrement=True)
    name = Column(String(100))
    transcripts = relationship('DTranscript', backref='chrm')
    # source_id = Column(Integer, ForeignKey('source.id'))    # primary key?
    # refversion_id = Column(Integer, ForeignKey('refversion.id'))
    # beg = Column(Integer, index=True)
    # end = Column(Integer, index=True)

# class Feature(Base):

#     __tablename__ = 'feature'
#     __table_args__ = {'mysql_engine':'InnoDB'}
#     id = Column(Integer, primary_key = True, autoincrement=True)
#     ftype = Column(Integer, ForeignKey('feature_type.id'))

class DGene(Base):

    __tablename__ = 'gene'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key = True, autoincrement=True)
    name = Column(String(100), primary_key=True)
    transcripts = relationship('DTranscript', backref='gene')

class DTranscriptType(Base):

    __tablename__ = 'transcript_type'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    transcripts = relationship('DTranscript', backref='transcript_type')

class DTranscript(Base):

    __tablename__ = 'transcript'
    __table_args__ = {'mysql_engine':'InnoDB'}
    id = Column(Integer, autoincrement=True, primary_key = True)
    name = Column(String(30))
    transcript_type_id = Column(Integer, ForeignKey('transcript_type.id'))
    chrm_id = Column(Integer, ForeignKey('chromosome.id'))
    source_id = Column(Integer, ForeignKey('source.id'))    # primary key?
    refversion_id = Column(Integer, ForeignKey('refversion.id'))
    beg = Column(Integer, index=True)
    end = Column(Integer, index=True)
    cds_beg = Column(Integer)
    cds_end = Column(Integer)
    gene_id = Column(Integer, ForeignKey('gene.id'))
    version = Column(Integer, default=0)
    strand = Column(Boolean)
    exons = relationship('DExon', backref='transcript')

class DExon(Base):
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
