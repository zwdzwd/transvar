"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Wanding Zhou, Tenghui Chen, Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

from .transcripts import *
from .localdb import TransVarDB
from . import parser
from pickle import load

class AnnoDB():

    """ AnnoDB keeps a collection of TransVarDB """
    
    def __init__(self, args, config):

        if args.refversion:
            self.rv = args.refversion
        elif 'refversion' in config.defaults():
            self.rv = config.get('DEFAULT', 'refversion')
        else:
            err_warn('please specify reference version, either in transvar.cfg, or as argument --refversion')
            sys.exit(1)
        
        replace_defaults(args, config)
        
        faidx.init_refgenome(args.reference if args.reference else None)
        self.session = None
        self.name2gene = {}
        self.name2trnx = {}
        self.thash = None

        self.dbs = []
        if args.ensembl:
            self.dbs.append(TransVarDB(args.ensembl, source='Ensembl'))
        if args.gencode:
            self.dbs.append(TransVarDB(args.gencode, source='GENCODE'))
        if args.kg:
            self.dbs.append(TransVarDB(args.kg, source='KnownGene'))
        if args.ucsc:
            self.dbs.append(TransVarDB(args.ucsc, source='UCSCRefGene'))
        if args.refseq:
            self.dbs.append(TransVarDB(args.refseq, source='RefSeq'))
        if args.ccds:
            self.dbs.append(TransVarDB(args.ccds, source='CCDS'))
        if args.aceview:
            self.dbs.append(TransVarDB(args.aceview, source='AceView'))
        if args.kg:
            self.dbs.append(TransVarDB(args.kg, source='KnownGene'))

        if args.uniprot:
            idmap = load(open(args.uniprot))
            for db in self.dbs:
                db.idmap = idmap
        self.config = config
        self.args = args
        self.resources = {}
        self.init_resource()

        # The following in-memory processing is unfinished
        # One needs to handle integratively transcript ID
        # input and gene ID input. NOT necessary.
        # The speed-up may not worth the effort.
        if args.mem:
            for db in self.dbs:
                db.parse_all(self.name2gene, self.name2trnx)

            self.thash = THash()
            genes = set(self.name2gene.values())
            for g in genes:
                for t in g.tpts:
                    self.thash.insert(t)

    def init_resource(self):
        """ init features and other annotation resources """
        for rname in ['dbsnp']:
            if self.config.has_option(self.rv, 'dbsnp'):
                import tabix
                self.resources['dbsnp'] = tabix.open(self.config.get(self.rv, 'dbsnp'))

        self.features = []
        for rname in self.config.options(self.rv):
            featdb =  self.config.get(self.rv, rname)
            if featdb.endswith('.featuredb'):
                self.features.append((rname,tabix.open(featdb)))

    def query_feature(self, r, chrm, beg, end):
        """ find all the dbsnp in a range """
        for rname, feat in self.features:
            for fields in tabix_query(feat, chrm, int(beg), int(end)):
                r.append_info('[feature:%s]=%s|%s:%s_%s' %
                              (rname, fields[3], fields[0], fields[1], fields[2]))

    def _query_dbsnp_(self, chrm, beg, end, ref=None, alt=None):

        dbsnps = []
        if 'dbsnp' in self.resources:
            if beg == end and (alt is None or len(alt)==1): # SNV
                ret = tabix_query(
                    self.resources['dbsnp'], normalize_chrm_dbsnp(chrm), int(beg), int(end))
                for fields in ret:
                    if int(fields[1]) != int(beg):
                        continue
                    alts = fields[4].split(',')
                    if ref is not None and ref != fields[3]:
                        continue
                    if alt is not None and alt not in alts:
                        continue
                    if alt is None:
                        for alt in alts:
                            dbsnps.append('%s(%s:%s%s>%s)' % (fields[2], chrm, fields[1], fields[3], alt))
                    else:
                        dbsnps.append('%s(%s:%s%s>%s)' % (fields[2], chrm, fields[1], fields[3], alt))
            else:               # indels and mnv
                ret = tabix_query(
                    self.resources['dbsnp'], normalize_chrm_dbsnp(chrm), int(beg)-1, int(end))
                for fields in ret:
                    if int(fields[1]) != int(beg)-1:
                        continue
                    if len(fields[3]) == 1 and len(fields[4]) == 1: # SNP from dbSNP
                        continue
                    alts = [_[1:] for _ in fields[4].split(',')]
                    if ref is not None and ref != fields[3][1:]:
                        continue
                    if alt is not None and alt not in alts:
                        continue
                    if alt is None:
                        for alt in alts:
                            dbsnps.append('%s(%s:%s%s>%s)' % (fields[2], chrm, fields[1], fields[3], fields[3][0]+alt))
                    else:
                        dbsnps.append('%s(%s:%s%s>%s)' % (fields[2], chrm, fields[1], fields[3], fields[3][0]+alt))
        return dbsnps

    def query_dbsnp_range(self, r, beg, end, alt):

        """ find all the dbsnp in a range """
        dbsnps = self._query_dbsnp_(r.chrm, beg, end, alt=alt)
        if dbsnps:
            r.append_info('dbsnp='+','.join(dbsnps))
        
    def query_dbsnp_codon(self, r, codon, taa_alt):

        """ find all the dbsnp in a codon """
        dbsnps = []
        for tnuc_altseq in reverse_codon_table[taa_alt]:
            subs = []
            for i, (tnuc_refbase, tnuc_altbase) in enumerate(zip(codon.seq, tnuc_altseq)):
                if tnuc_refbase != tnuc_altbase:
                    subs.append((i, tnuc_altbase))
                    
            if len(subs) == 1:
                i, tnuc_altbase = subs[0]
                gnuc_pos = codon.tloc(i)
                gnuc_alt = tnuc_altbase if codon.strand == '+' else complement(tnuc_altbase)
                dbsnps.extend(self._query_dbsnp_(r.chrm, gnuc_pos, gnuc_pos, alt=gnuc_alt))

            if len(subs) == 2:
                i1 = subs[0][0]
                i2 = subs[1][0]
                if i1+1 == i2 and abs(codon.tloc(i1)-codon.tloc(i2)) == 1:
                    gnuc_beg = min(codon.tloc(i1), codon.tloc(i2))
                    gnuc_end = max(codon.tloc(i1), codon.tloc(i2))
                    tnuc_alt = subs[0][1]+subs[1][1]
                    gnuc_alt = tnuc_alt if codon.strand == '+' else reverse_complement(tnuc_alt)
                    dbsnps.extend(self._query_dbsnp_(r.chrm, gnuc_beg, gnuc_end, gnuc_alt))

            if len(subs) == 3:
                i1 = subs[0][0]
                i2 = subs[2][0]
                if i1+2 == i2 and abs(codon.tloc(i1)-codon.tloc(i2)) == 2:
                    gnuc_beg = min(codon.tloc(i1), codon.tloc(i2))
                    gnuc_end = max(codon.tloc(i1), codon.tloc(i2))
                    tnuc_alt = codon.seq
                    gnuc_alt = codon.seq if codon.strand == '+' else reverse_complement(codon.seq)
                    dbsnps.extend(self._query_dbsnp_(r.chrm, gnuc_beg, gnuc_end, gnuc_alt))

        if dbsnps:
            r.append_info('dbsnp='+','.join(dbsnps))
        
    def query_dbsnp(self, r, pos, ref=None, alt=None):

        dbsnps = self._query_dbsnp_(r.chrm, pos, pos, ref, alt)
        if dbsnps:
            r.append_info('dbsnp='+','.join(dbsnps))

    def get_gene(self, name):

        for db in self.dbs:
            for g in db.get(name):
                yield g

    def get_transcripts(self, chrm, beg, end=None, flanking=0):

        for db in self.dbs:
            for t in db.get_by_loc(chrm, beg, end, flanking):
                yield t

    def get_closest_transcripts_upstream(self, chrm, pos):

        max_end = -1
        max_t = None
        for db in self.dbs:
            t = db.get_closest_upstream(chrm, pos)
            if t is not None and t.end > max_end:
                max_t = t

        return max_t

    def get_closest_transcripts_downstream(self, chrm, pos):

        min_beg = 30000000000
        min_t = None
        for db in self.dbs:
            t = db.get_closest_downstream(chrm, pos)
            if t is not None and t.beg < min_beg:
                min_t = t

        return min_t

    def get_closest_transcripts(self, chrm, beg, end):
        """ closest transcripts upstream and downstream """
        return (self.get_closest_transcripts_upstream(chrm, beg),
                self.get_closest_transcripts_downstream(chrm, end))
