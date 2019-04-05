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

from . import argparse
from . import parser
import sys
import re, os
from .utils import *
from .transcripts import *
from pickle import load, dump
from . import faidx
from . import tabix
import subprocess

## tabix and bgzip is now a dependency, but only when featuredb and indexing is used
# tabix_path = '%s/tabix' % os.path.abspath(os.path.dirname(__file__))
tabix_path = 'tabix'
# bgzip_path = '%s/bgzip' % os.path.abspath(os.path.dirname(__file__))
bgzip_path = 'bgzip'

p_trxn_version=re.compile(r'(.*)\.(\d+)$')
class TransVarDB():

    """ hold transcripts and genes
    In TransVar, transcripts are indexed in two ways to allow both access from name and from coordinates.
    *.transvardb is ordered by gene name and index contains location to the first item
    *.transvardb.loc_idx is a bed-like file ordered by coordinates
    Different from TransVarDB, FeatureDB is only indexed by coordinates.
    """

    def __init__(self, dbfn=None, source=None):

        self.name2gene = {}
        self.idmap = {}
        if dbfn is None: return

        self.dbfn = dbfn
        self.dbfh = open(self.dbfn, 'rt')

        idxfn = dbfn+'.gene_idx'
        self.gene_idx = load(open(idxfn, 'rb'))

        idxfn = dbfn+'.trxn_idx'
        self.trnx_idx = load(open(idxfn, 'rb'))

        self.alias_idx = None
        self.loc_idx = None
        self.source = source

    ##########################
    # parsers for transvardb #
    ##########################

    def parse_trnx(self, gname=None):

        """ parse name-indexed transcript file
        .transvardb file

        if gname is None, read one transcript
        if gname is not None, read the number of transcripts with
        the given gname

        parsing starts from dbfh current location
        """
        for line in self.dbfh:
            fields = line.strip('\n').split('\t')
            if gname is not None and fields[0] != gname:
                break
            t = Transcript()
            t.gene_name = fields[0] # important when name == None
            t.name = fields[1]
            t.version = int(fields[2])
            t.transcript_type = fields[3]
            t.beg = int(fields[4])
            t.end = int(fields[5])
            t.chrm = fields[6]
            t.strand = fields[7]
            t.cds_beg = int(fields[8])
            t.cds_end = int(fields[9])
            t.exons = eval(fields[10])
            if fields[11]:
                t.aliases = fields[11].split(';')
            t.gene_dbxref = fields[12]
            t.source = self.source
            yield t
            if gname is None:
                break

    def parse_trnx_loc(self, fields):

        """ parse location-indexed transcript file
        .transvardb.loc_idx file
        return only 1 line
        """
        t = Transcript()
        t.chrm = fields[0]
        t.beg = int(fields[1])
        t.end = int(fields[2])
        t.gene_name = fields[3] # important when name == None
        t.name = fields[4]
        t.version = int(fields[5])
        t.transcript_type = fields[6]
        t.strand = fields[7]
        t.cds_beg = int(fields[8])
        t.cds_end = int(fields[9])
        t.exons = eval(fields[10])
        if fields[11]:
            t.aliases = fields[11].split(';')
        t.gene_dbxref = fields[12]
        t.source = self.source
        return t

    def parse_all(self, name2gene, name2trnx):

        """ parse the whole location-indexed transcript file
        .transvardb.loc_idx file
        this is for in-memory processing
        """

        self.dbfh.seek(0)
        for line in self.dbfh:
            fields = line.strip('\n').split('\t')
            t = parse_trnx_loc(fields)
            if t.gene_name in name2gene:
                g = name2gene[t.gene_name]
            else:
                g = Gene(t.gene_name)
                name2gene[t.gene_name] = g
            g.link_t(t)
            if t.name not in name2trnx:
                name2trnx[t.name] = []
            name2trnx[t.name].append(t)


    ########################################################
    # retrieve transcripts by gene name or transcript name #
    ########################################################

    def get(self, name, search_alias = True, strictversion = False):

        """ get by either gene name or transcript name
        use the following order:
        1) gene name
        2) trnx name
        3) trxn name with version info
        4) alias from idmap
        """

        no_hit = True
        for g in self.get_by_gene(name):
            no_hit = False
            yield g

        if no_hit:
            g = self.get_by_trnx(name, strictversion)
            if g is not None:
                no_hit = False
                yield g

        ## deprecates alias_idx and uses idmap
        # # try alias_idx
        # if no_hit and self.alias_idx is not None and search_alias:
        #     _name = name
        #     m = p_trxn_version.match(name)
        #     if m:
        #         _name = m.group(1)

        #     for g in self.get_by_alias(_name):
        #         no_hit = False
        #         yield g

        # try idmap
        if no_hit and len(self.idmap) > 0 and search_alias:
            if name in self.idmap:
                for name2 in self.idmap[name]:
                    for g in self.get(name2, search_alias = False):
                        yield g

    def get_by_gene(self, name):

        """ get by gene name """
        if name in self.gene_idx:
            pos = self.gene_idx[name]
            self.dbfh.seek(pos)
            g = Gene(name)
            for t in self.parse_trnx(gname=name):
                g.link_t(t)
            yield g

    def get_by_trnx(self, name, strictversion = False):

        """ read in a gene by transcript name 
        the gene has only the transcript(s) with that transcript name
        """
        m = p_trxn_version.match(name)
        if m:
            name = m.group(1)
            version = int(m.group(2))
        else:
            version = None # no version info

        if name not in self.trnx_idx:
            return None

        poses = self.trnx_idx[name] # transcript ID might not be unique
        g = None
        for pos in poses:
            self.dbfh.seek(pos)
            t = next(self.parse_trnx(), None)

            if t is None:
                return None
            elif g is None:
                g = Gene(t.gene_name)
            t.gene = g
            
            if (strictversion and 
                version is not None and 
                t.version != version):
                continue
            g.tpts.append(t)

        return g

    def get_by_alias(self, alias):

        """ read a gene by alias of transcripts """
        if self.alias_idx is None and os.path.exists(self.dbfn+'.alias_idx'):
            self.alias_idx = load(open(self.dbfn+'.alias_idx', 'rb'))

        if self.alias_idx is not None and alias in self.alias_idx:
            poses = self.alias_idx[alias]
            name2gene = {}
            for pos in poses:
                self.dbfh.seek(pos)
                t = next(self.parse_trnx(), None)
                if t is None:
                    continue
                elif t.gene_name in name2gene:
                    g = name2gene[t.gene_name]
                else:
                    g = Gene(t.gene_name)
                    name2gene[g.name] = g
                g.link_t(t)

            for g in name2gene.values():
                yield g

    ##############################################
    # retrieve transcript by genomic coordinates #
    ##############################################

    def _ensure_loc_idx(self):
        if self.loc_idx is None:
            idx_fn = self.dbfn+'.loc_idx'
            if not os.path.exists(idx_fn):
                err_die("Missing location index. Consider rerunning the transvar index command")
            self.loc_idx = tabix.open(idx_fn)

    def _iloc_query(self, chrm, beg, end):

        self._ensure_loc_idx()
        return tabix_query(self.loc_idx, chrm, beg, end)

    def get_by_loc(self, chrm, beg, end=None, flanking=0):

        """ get transcript if between begin and end """
        self._ensure_loc_idx()
        if not end: end = beg
        chrm = normalize_chrm(chrm)
        for fields in self._iloc_query(chrm,beg-flanking,end+flanking):
            yield self.parse_trnx_loc(fields)

    def get_closest_upstream(self, chrm, pos):
        pos = int(pos)
        chrm = normalize_chrm(chrm)
        s = 50000
        for p in range(pos, -1, -s):
            fs = [f for f in self._iloc_query(chrm, p-s, p) if int(f[2])<pos]
            if fs:
                fs.sort(key=lambda f:int(f[2]), reverse=True)
                return self.parse_trnx_loc(fs[0])
        return None

    def get_closest_downstream(self, chrm, pos):
        pos = int(pos)
        chrm = normalize_chrm(chrm)
        s = 50000
        chrmlen = faidx.refgenome.chrm2len(chrm)
        for p in range(pos, chrmlen, s):
            fs = [f for f in self._iloc_query(chrm, p, p+s) if int(f[1])>pos]
            if fs:
                fs.sort(key=lambda f:int(f[1]))
                return self.parse_trnx_loc(fs[0])
        return None

    def get_closest(self, chrm, beg, end):
        """ closest transcripts upstream and downstream """
        return (self.get_closest_upstream(chrm, beg), self.get_closest_downstream(chrm, end))


    #####################################
    # index transcripts from raw files ##
    #####################################

    def index(self, raw_fns):

        # each class that subclassed TransVarDB should have parse_raw
        self.parse_raw(*raw_fns)
        # results are stored in self.name2gene and self.idmap
        
        # set cds_beg and cds_end
        set_cds_boundary(self.name2gene)

        ## .transvardb
        names = sorted(self.name2gene.keys())
        dbfn = raw_fns[0]+'.transvardb'
        dbfh = open(dbfn, 'wt')
        gene_idx = {}
        trnx_idx = {}
        # alias_idx = {}          # hold transcript aliases
        tpts = []
        for name in names:
            g = self.name2gene[name]
            for t in g.tpts:
                t.gene_name = g.name
                tpts.append((t.chrm, t.beg, t.end, t))
                pos = dbfh.tell()
                if g.name not in gene_idx: # first location, each gene record one position
                    gene_idx[g.name] = pos

                if t.name in trnx_idx: # precise location
                    trnx_idx[t.name].append(pos)
                else:
                    trnx_idx[t.name] = [pos]

                # for alias in t.aliases:
                #     if alias in alias_idx:
                #         alias_idx[alias].append(pos)
                #     else:
                #         alias_idx[alias] = [pos]

                dbfh.write('%s\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n' %
                           (g.name, t.name, t.version, t.transcript_type, t.beg, t.end, t.chrm,
                            t.strand, t.cds_beg, t.cds_end, t.exons, ';'.join(t.aliases), g.dbxref))

        ############################################
        ## .gene_idx - index gene name
        ############################################
        idxfn = dbfn+'.gene_idx'
        dump(gene_idx, open(idxfn, 'wb'), 2)

        ############################################
        ## .trxn_idx - index transcript name
        ############################################
        idxfn = dbfn+'.trxn_idx'
        dump(trnx_idx, open(idxfn, 'wb'), 2)

        ############################################
        ## .?.idmap_idx - mappings to gene names
        ## or transcript IDs
        ############################################
        for map_name, mapping in self.idmap.items():
            idxfn = dbfn + '.' + map_name + '.idmap_idx'
            dump(mapping, open(idxfn, 'wb'), 2)
        
        # if len(alias_idx) > 0:
        #     idxfn = dbfn+'.alias_idx'
        #     dump(alias_idx, open(idxfn, 'wb'), 2)

        ############################################
        ## .loc_idx - tab-index genomic locations
        ############################################
        idxfn = dbfn+'.loc_idx'
        tpts.sort()
        s = ''
        for chrm, beg, end, t in tpts:
            s += '%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n' % (
                t.chrm, t.beg, t.end, t.gene_name, t.name, t.version, t.transcript_type,
                t.strand, t.cds_beg, t.cds_end, t.exons, ';'.join(t.aliases), t.gene.dbxref)

        ## call external tabix
        with open(idxfn, 'wb') as fh:
            p = subprocess.Popen(bgzip_path, stdout=fh, stdin=subprocess.PIPE)
            p.communicate(input=s.encode('utf-8'))

        subprocess.check_call([tabix_path, '-p', 'bed', idxfn])

class FeatureDB():

    def parse_bed(self, bed_fn, db_fn):

        """ bed format indexing takes only the first four columns,
        the annotation is the fourth column
        """
        bed_fh = opengz(bed_fn)
        with open(db_fn+'.presort', 'wt') as outfile:
            for line in bed_fh:
                fields = line.strip('\n').split('\t')
                if len(fields) < 4:
                    continue
                outfile.write('%s\t%s\t%s\t%s\n' % (
                    normalize_chrm(fields[0]), fields[1], fields[2], fields[3]))

    def parse_gff(self, gff_fn, db_fn):

        """ GFF: seqname, source, feature, start, end, score, strand, frame, attribute
        indexing made a bed file with seqname, start, end, feature
        """
        gff_fh = opengz(gff_fn)
        with open(db_fn+'.presort', 'wt') as outfile:
            for line in gff_fh:
                fields = line.strip('\n').split('\t')
                if len(fields) < 4:
                    continue
                outfile.write('%s\t%s\t%s\t%s\n' % (
                    normalize_chrm(fields[0]), fields[3], fields[4], fields[2]))

    def parse_vcf(self, vcf_fn, db_fn):

        """ VCF: #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INF,
        indexing made a bed file with CHROM, POS, POS+len(REF), ID|REF|ALT """
        vcf_fh = opengz(vcf_fn)
        with open(db_fn+'.presort', 'wt') as outfile:
            for line in vcf_fh:
                if line.startswith('#'):
                    continue
                fields = line.strip('\n').split('\t')
                if len(fields) < 7:
                    continue
                outfile.write('%s\t%s\t%d\t%s|%s|%s\n' % (
                    normalize_chrm(fields[0]), fields[1], int(fields[1])+len(fields[3]),
                    fields[2],fields[3],fields[4]))

    def index(self, fn, raw_format, is_sorted):
        db_fn = fn+'.featuredb'
        if raw_format == 'bed':
            self.parse_bed(fn, db_fn)
        elif raw_format == 'vcf':
            self.parse_vcf(fn, db_fn)
        elif raw_format == 'gff':
            self.parse_gff(fn, db_fn)
        else:
            raise Exception('Unknown format, must be a bug.\n')

        if is_sorted:
            subprocess.check_call(['mv', db_fn+'.presort', db_fn+'.sort'])
        else:
            with open(db_fn+'.sort', 'wt') as outfile:
                subprocess.check_call(['sort', '-T', os.path.dirname(os.path.abspath(db_fn)),
                                       '-k', '1,1', '-k', '2,2n', db_fn+'.presort'], stdout=outfile)
        with open(db_fn, 'wt') as outfile:
            subprocess.check_call([bgzip_path, '-c', db_fn+'.sort'], stdout=outfile)
        subprocess.check_call([tabix_path, '-p', 'bed', db_fn])

        # clean
        subprocess.check_call(['rm', '-f', db_fn+'.presort'])
        subprocess.check_call(['rm', '-f', db_fn+'.sort'])

def set_cds_boundary(name2gene):

    for g in name2gene.values():
        for t in g.tpts:
            t.exons.sort()
            if not (hasattr(t, 'cds_beg') and hasattr(t, 'cds_end')):
                if t.cds:
                    t.cds.sort()
                    t.cds_beg = t.cds[0][0]
                    t.cds_end = t.cds[-1][1]
                elif hasattr(t,'beg') and hasattr(t,'end'):
                    t.cds_beg = t.beg
                    t.cds_end = t.end
                else:
                    t.cds_beg = t.exons[0][0]
                    t.cds_end = t.exons[-1][1]

            if len(t.exons) == 0:
                t.exons = t.cds[:]

def mapping_append(m, k, v):

    if k in m:
        if v not in m[k]:
            m[k].append(v)
    else:
        m[k] = {v}

def index_idmap(idmap_fn, idx_fh):

    fh = opengz(idmap_fn)
    idmap = {}
    for i, line in enumerate(fh):
        if line.startswith('#'): continue
        fields = line.strip('\n').split('\t')
        mapping_append(idmap, fields[0], fields[1])

    dump(idmap, idx_fh, 2)

    return

class EnsemblDB(TransVarDB):

    def __init__(self, dbfn=None):
        # super(self.__class__, self).__init__(dbfn)
        TransVarDB.__init__(self, dbfn)

    def parse_raw(self, gtf_fn):
        """
        This parses the new GTF after or equal to hg19.
        The function does not handle hg18.
        gtf file is gffv2
        parser does not assume order in the GTF file
        
        The result is stored to self.name2gene
        """

        gtf_fh = opengz(gtf_fn)

        ## the following is a workaround
        if (gtf_fn.startswith('Mus_musculus.NCBIM37')):
            return self.parse_raw0(gtf_fn)
        
        id2ent = {}
        cnt = 0
        
        ## id mappings
        geneID_to_geneName = {}
        proteinID_to_transcriptID = {}
        
        for i, line in enumerate(gtf_fh):
            if line.startswith('#'): continue
            fields = line.strip('\n').split('\t')
            info = dict(re.findall(r'\s*([^";]*) "([^"]*)";', fields[8]))
            # info = dict([_.split('=') for _ in fields[8].split(';')])
            if fields[2] == 'gene':
                gene_id = info['gene_id']
                if gene_id not in id2ent:
                    id2ent[gene_id] = Gene(gene_type=info['gene_biotype'])
                g = id2ent[gene_id]
                if 'gene_name' in info:
                    g.name = info['gene_name'].upper()
                else:
                    g.name = gene_id
                if g.name not in self.name2gene:
                    self.name2gene[g.name] = g

                ## mapping from geneID to geneName
                mapping_append(geneID_to_geneName, gene_id, g.name)
                        
                g.beg = int(fields[3])
                g.end = int(fields[4])

            elif fields[2] == 'transcript':

                # there exits two transcript format in ensembl gtf
                # the old version has no 'transcript_biotype'
                # the equivalent transcript_biotype is fields[1]
                tid = info['transcript_id']
                if tid not in id2ent:
                    transcript_type = info['transcript_biotype'] if 'transcript_biotype' in info else fields[1]
                    id2ent[tid] = Transcript(transcript_type=transcript_type)
                t = id2ent[tid]
                t.chrm = normalize_chrm(fields[0])
                t.strand = fields[6]
                t.beg = int(fields[3])
                t.end = int(fields[4])
                t.name = info['transcript_id']
                gene_id = info['gene_id']
                if gene_id not in id2ent:
                    id2ent[gene_id] = Gene(gene_type=info['gene_biotype'])
                t.gene = id2ent[gene_id]
                t.gene.tpts.append(t)
                t.source = 'Ensembl'
                cnt += 1
            elif fields[2] == 'exon':
                tid = info['transcript_id']
                if tid not in id2ent:
                    transcript_type = info['transcript_biotype'] if 'transcript_biotype' in info else fields[1]
                    id2ent[tid] = Transcript(transcript_type=transcript_type)
                t = id2ent[tid]
                t.exons.append((int(fields[3]), int(fields[4])))
            elif fields[2] == 'CDS':
                tid = info['transcript_id']
                if tid not in id2ent:
                    transcript_type = info['transcript_biotype'] if 'transcript_biotype' in info else fields[1]
                    id2ent[tid] = Transcript(transcript_type=transcript_type)
                t = id2ent[tid]
                t.cds.append((int(fields[3]), int(fields[4])))
                if 'protein_id' in info:
                    ## mapping from geneID to geneName
                    mapping_append(proteinID_to_transcriptID, info['protein_id'], tid)
                    if info['protein_id'] not in t.aliases:
                        t.aliases.append(info['protein_id'])

        self.idmap = {
            'protein_id': proteinID_to_transcriptID,
            'gene_id': geneID_to_geneName}

        err_print("loaded %d transcripts from Ensembl GTF file." % cnt)

        return

    #################################################################
    # TransVar stops supporting hg18 raw format by default any more
    #################################################################
    def parse_raw0(self, gtf_fn):
        """
        This parses the old ensembl GTF before or equal to hg18.
        The function does not handle hg19 or later.
        """

        gtf_fh = opengz(gtf_fn)
        tid2transcript = {}
        cnt = 0
        proteinID_to_transcriptID = {}
        for line in gtf_fh:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            info = dict(re.findall(r'\s*([^"]*) "([^"]*)";', fields[8]))
            if fields[2] == "exon":
                if info['transcript_id'] in tid2transcript:
                    t = tid2transcript[info['transcript_id']]
                else:
                    t = Transcript(transcript_type=fields[1])
                    t.chrm = normalize_chrm(fields[0])
                    t.strand = fields[6]
                    t.name = info['transcript_id']
                    tid2transcript[t.name] = t
                    if info['gene_name'] in self.name2gene:
                        g = self.name2gene[info['gene_name']]
                    else:
                        g = Gene()
                        g.name = info['gene_name']
                        self.name2gene[g.name] = g
                    t.gene = g
                    g.tpts.append(t)
                    t.source = 'Ensembl'
                    cnt += 1
                t.exons.append((int(fields[3]), int(fields[4])))
            elif fields[2] == 'CDS':
                tid = info['transcript_id']
                if tid in tid2transcript:
                    t = tid2transcript[tid]
                else:
                    t = Transcript(transcript_type=fields[1])
                    t.chrm = normalize_chrm(fields[0])
                    t.strand = fields[6]
                    t.name = tid
                    tid2transcript[t.name] = t
                    if info['gene_name'] in self.name2gene:
                        g = self.name2gene[info['gene_name']]
                    else:
                        g = Gene()
                        g.name = info['gene_name']
                        self.name2gene[g.name] = g
                    t.gene = g
                    g.tpts.append(t)
                    t.source = 'Ensembl'
                    cnt += 1
                t.cds.append((int(fields[3]), int(fields[4])))
                if 'protein_id' in info:
                    ## mapping from geneID to geneName
                    mapping_append(proteinID_to_transcriptID, info['protein_id'], tid)
                    if info['protein_id'] not in t.aliases:
                        t.aliases.append(info['protein_id'])

        self.idmap = {
            'protein_id': proteinID_to_transcriptID}

        for t in list(tid2transcript.values()):
            t.exons.sort()
            t.beg = t.exons[0][0]
            t.end = t.exons[-1][1]

        err_print("loaded %d transcripts from Ensembl GTF file." % cnt)

class CCDSDB(TransVarDB):

    def __init__(self, dbfn=None):
        TransVarDB.__init__(self, dbfn)

    def parse_raw(self, ccds_fn):

        ccds_fh = open(ccds_fn, 'rt')
        ccds_fh.readline()
        cnt = 0
        for line in ccds_fh:
            fields = line.strip('\n').split('\t')
            if fields[5] != 'Public':
                continue
            gene_name = fields[2].upper()
            if gene_name not in self.name2gene:
                self.name2gene[gene_name] = Gene(name=gene_name)

            g = self.name2gene[gene_name]
            t = Transcript()
            t.chrm = normalize_chrm(fields[0])
            t.strand = fields[6]
            t.cds_beg = int(fields[7])+1
            t.cds_end = int(fields[8])+1

            # without UTR information, take CDS boundary as the exon boundary
            t.beg = t.cds_beg
            t.end = t.cds_end

            t.name = fields[4]
            m = p_trxn_version.match(t.name)
            if m:
                t.name = m.group(1)
                t.version = int(m.group(2))

            # note that CCDS do not annotate UTR, so all the exons are equivalently cds
            t.exons = [(int(b)+1, int(e)+1) for b,e in re.findall(r"[\s\[](\d+)-(\d+)[,\]]", fields[9])]
            t.source = 'CDDS'
            t.gene = g
            g.tpts.append(t)
            cnt += 1

        err_print("loaded %d transcripts from CCDS table." % cnt)


class RefSeqDB(TransVarDB):

    def __init__(self, dbfn=None):
        TransVarDB.__init__(self, dbfn)

    def parse_raw(self, gff_fn):

        id2ent = {}
        gff_fh = opengz(gff_fn)
        reg = None
        cnt = 0

        ## id mappings
        geneID_to_geneName = {}
        HGNC_to_geneName = {}
        proteinID_to_transcriptID = {}
        # TODO: HRPD and MIM
        
        for line in gff_fh:
            if line.startswith('#'): continue
            fields = line.strip('\n').split('\t')
            if (len(fields) < 9):
                continue
            # print line.strip()
            info = dict([_.split('=') for _ in fields[8].split(';')])
            if fields[2] == 'region':
                if 'chromosome' in info:
                    reg = Region(info['chromosome'], int(fields[3]), int(fields[4]))

                if 'map' in info and info['map']=='unlocalized':
                    reg.unlocalized = True

                # else:
                # reg = None
            elif (reg and fields[2] == 'gene' and
                  ('pseudo' not in info or info['pseudo'] != 'true')):

                if reg.unlocalized:
                    continue

                gene_name = info['Name'].upper()
                if gene_name in self.name2gene:
                    g = self.name2gene[gene_name]
                    if hasattr(g, '_gene_id') and g._gene_id != info['ID']:
                        # if a gene_name appears twice, then all
                        # the subsequent occurrences are all ignored.
                        continue
                else:
                    g = Gene(name=gene_name)
                    self.name2gene[gene_name] = g
                g._gene_id = info['ID']
                g.beg = int(fields[3])
                g.end = int(fields[4])
                id2ent[info['ID']] = g
                if 'Dbxref' in info:
                    g.dbxref = info['Dbxref']
                    _dbxref = dict(re.findall(r'([^:,\s]*):([^;,\s]*)', g.dbxref))
                    if 'GeneID' in _dbxref:
                        mapping_append(geneID_to_geneName, _dbxref['GeneID'], gene_name)
                    if 'HGNC' in _dbxref:
                        # import pdb; pdb.set_trace()
                        mapping_append(HGNC_to_geneName, _dbxref['HGNC'], gene_name)

            elif (fields[2] in ['mRNA', 'ncRNA', 'rRNA', 'tRNA']
                  and 'Parent' in info and info['Parent'] in id2ent):

                if reg.unlocalized:
                    continue

                if fields[2] == 'mRNA':
                    fields[2] = 'protein_coding'
                if fields[2] == 'ncRNA' and 'ncrna_class' in info:
                    fields[2] = info['ncrna_class']
                t = Transcript(transcript_type=fields[2])
                t.chrm = normalize_chrm(reg.name)
                t.strand = fields[6]
                t.beg = int(fields[3])
                t.end = int(fields[4])
                t.name = info['Name'] if 'Name' in info else info['product']
                m = p_trxn_version.match(t.name)
                if m:
                    t.name = m.group(1)
                    t.version = int(m.group(2))
                t.gene = id2ent[info['Parent']]
                t.gene.tpts.append(t)
                t.source = 'RefSeq'
                id2ent[info['ID']] = t
                cnt += 1

            elif fields[2] == 'exon' and info['Parent'] in id2ent:

                if reg.unlocalized:
                    continue

                t = id2ent[info['Parent']]
                if (isinstance(t, Gene)):
                    g = t
                    if not hasattr(g, 'gene_t'):
                        g.gene_t = Transcript()
                        g.tpts.append(g.gene_t)
                        g.gene_t.chrm = normalize_chrm(reg.name)
                        g.gene_t.strand = fields[6]
                        g.gene_t.gene = g
                        g.gene_t.beg = g.beg
                        g.gene_t.end = g.end
                        g.gene_t.source = 'RefSeq'
                        cnt += 1
                    t = g.gene_t
                t.exons.append((int(fields[3]), int(fields[4])))

            elif fields[2] == 'CDS' and info['Parent'] in id2ent:

                if reg.unlocalized:
                    continue

                t = id2ent[info['Parent']]
                if (isinstance(t, Gene)):
                    g = t
                    if not hasattr(g, 'gene_t'):
                        g.gene_t = Transcript()
                        g.tpts.append(g.gene_t)
                        g.gene_t.chrm = normalize_chrm(reg.name)
                        g.gene_t.strand = fields[6]
                        g.gene_t.gene = g
                        g.gene_t.beg = g.beg
                        g.gene_t.end = g.end
                        g.gene_t.source = 'RefSeq'
                        cnt += 1
                    t = g.gene_t
                t.cds.append((int(fields[3]), int(fields[4])))
                if 'protein_id' in info:
                    protein_id = info['protein_id']
                    m = p_trxn_version.match(protein_id)
                    if m:
                        protein_id = m.group(1)
                        protein_version = int(m.group(2)) # TODO: don't discard this info
                    if protein_id not in t.aliases:
                        t.aliases.append(protein_id)
                        mapping_append(proteinID_to_transcriptID, protein_id, t.name)

        self.idmap = {
            'GeneID': geneID_to_geneName,
            'HGNC': HGNC_to_geneName,
            'protein_id': proteinID_to_transcriptID}

        err_print("loaded %d transcripts from RefSeq GFF3 file." % cnt)

class AceViewDB(TransVarDB):

    def __init__(self, dbfn=None):
        TransVarDB.__init__(self, dbfn)

    def parse_raw(self, aceview_gff_fn):

        id2tpt = {}
        aceview_fh = opengz(aceview_gff_fn)
        for line in aceview_fh:
            if line.startswith('#'): continue
            if line.startswith('AceView'): continue
            fields = line.strip('\n').split('\t')
            if len(fields) != 9:
                continue
            if len(fields) < 9: continue # the old transcript definition (hg18) from AceView is a bit corrupted.
            info = dict(re.findall(r'\s*(\S+) (\S+);', fields[8]))
            if fields[2] == 'CDS':
                gene_name = info['gene_id'].upper()
                if gene_name in self.name2gene:
                    g = self.name2gene[gene_name]
                else:
                    g = Gene(name=gene_name)
                    self.name2gene[gene_name] = g

                if info['transcript_id'] in id2tpt:
                    t = id2tpt[info['transcript_id']]
                else:
                    t = Transcript()
                    t.chrm = normalize_chrm(fields[0])
                    t.strand = fields[6]
                    t.name = info['transcript_id']
                    id2tpt[t.name] = t
                    t.gene = g
                    g.tpts.append(t)
                    t.source = 'AceView'

                t.cds.append((int(fields[3]), int(fields[4])))

            elif fields[2] == 'exon':
                gene_name = info['gene_id'].upper()
                if gene_name in self.name2gene:
                    g = self.name2gene[gene_name]
                else:
                    g = Gene(name=gene_name)
                    self.name2gene[gene_name] = g

                if info['transcript_id'] in id2tpt:
                    t = id2tpt[info['transcript_id']]
                else:
                    t = Transcript()
                    t.chrm = normalize_chrm(fields[0])
                    t.strand = fields[6]
                    t.name = info['transcript_id']
                    id2tpt[t.name] = t
                    t.gene = g
                    g.tpts.append(t)
                    t.source = 'AceView'

                t.exons.append((int(fields[3]), int(fields[4])))

        # skip transcripts without CDS, e.g., LOC391566.aAug10-unspliced
        for tid, t in id2tpt.items():
            if t.cds and t.exons:
                t.exons.sort()
                t.beg = t.exons[0][0]
                t.end = t.exons[-1][1]
            else:
                t.gene.tpts.remove(t)

        err_print("loaded %d transcripts from AceView GFF file." % len(id2tpt))

class GENCODEDB(TransVarDB):

    def __init__(self, dbfn=None):
        TransVarDB.__init__(self, dbfn)

    def parse_raw(self, gencode_fn):

        gtf_fh = opengz(gencode_fn)
        id2ent = {}
        cnt = 0
        
        for line in gtf_fh:
            if line.startswith('#'): continue
            fields = line.strip('\n').split('\t')
            info = dict(re.findall(r'\s*([^";]*) "([^"]*)";', fields[8]))
            if fields[2] == 'gene':
                gene_name = info['gene_name'].upper()
                gid = info['gene_id']
                if gene_name in self.name2gene:
                    g = self.name2gene[gene_name]
                    id2ent[gid] = g
                else:
                    if gid not in id2ent:
                        id2ent[gid] = Gene(name=gene_name, gene_type=info['gene_type'])
                    g = id2ent[gid]
                    self.name2gene[gene_name] = g
                g.beg = int(fields[3])
                g.end = int(fields[4])
                # if info['gene_type'] == 'pseudogene':
                #     g.pseudo = True
            elif fields[2] == 'transcript':
                tid = info['transcript_id']
                if tid not in id2ent:
                    id2ent[tid] = Transcript(transcript_type=info['transcript_type'])

                t = id2ent[tid]
                t.chrm = normalize_chrm(fields[0])
                t.strand = fields[6]
                t.beg = int(fields[3])
                t.end = int(fields[4])
                t.name = tid
                m = p_trxn_version.match(t.name)
                if m:
                    t.name = m.group(1)
                    t.version = int(m.group(2))
                gid = info['gene_id']
                if gid not in id2ent:
                    id2ent[gid] = Gene(gene_type=info['gene_type'])
                t.gene = id2ent[gid]
                t.gene.tpts.append(t)
                t.source = 'GENCODE'
                id2ent[t.name] = t
                cnt += 1
                if 'protein_id' in info:
                    if info['protein_id'] not in t.aliases:
                        t.aliases.append(info['protein_id'])
            elif fields[2] == 'exon':
                tid = info['transcript_id']
                if tid not in id2ent:
                    id2ent[tid] = Transcript(transcript_type=info['transcript_type'])
                t = id2ent[tid]
                t.exons.append((int(fields[3]), int(fields[4])))
            elif fields[2] == 'CDS':
                tid = info['transcript_id']
                if tid not in id2ent:
                    id2ent[tid] = Transcript(transcript_type=info['transcript_type'])
                if tid not in id2ent: id2ent[tid] = Transcript()
                t = id2ent[tid]
                t.cds.append((int(fields[3]), int(fields[4])))
                if 'protein_id' in info:
                    if info['protein_id'] not in t.aliases:
                        t.aliases.append(info['protein_id'])

        err_print("loaded %d transcripts from GENCODE GTF file." % cnt)

class UCSCRefGeneDB(TransVarDB):

    def __init__(self, dbfn=None):
        TransVarDB.__init__(self, dbfn)

    def parse_raw(self, map_file):
        """ start 1-based, end 1-based """
        cnt = 0
        for line in opengz(map_file):
            if line.startswith('#'): continue
            fields = line.strip('\n').split('\t')
            if fields[13] != 'cmpl' or fields[14] != 'cmpl':
                continue
            gene_name = fields[12].upper()
            if gene_name in self.name2gene:
                gene = self.name2gene[gene_name]
            else:
                gene = Gene(name=gene_name)
                self.name2gene[gene_name] = gene
            t = Transcript()
            t.name = fields[1]
            t.chrm = normalize_chrm(fields[2])
            t.strand = fields[3]
            t.beg    = int(fields[4])+1
            t.end    = int(fields[5])
            t.cds_beg = int(fields[6])+1
            t.cds_end = int(fields[7])
            t.source = 'UCSC_refGene'
            ex_begs, ex_ends = fields[9], fields[10]

            for ex_beg, ex_end in zip([int(x)+1 for x in ex_begs.strip(',').split(',')],
                                      list(map(int, ex_ends.strip(',').split(',')))):
                t.exons.append((ex_beg, ex_end))

            t.exons.sort() # keep exons sorted
            gene.tpts.append(t)
            t.gene = gene
            cnt += 1

        err_print('loaded %d transcripts from UCSC refgene.' % cnt)

class UCSCKnownGeneDB(TransVarDB):

    """ NOTE!!! knowngene gene name is considered an alias to the transcript ID.
    the recorded gene name is just the transcript ID """
    def __init__(self, dbfn=None):
        TransVarDB.__init__(self, dbfn)

    def parse_raw(self, kg_fn, alias_fn):

        kg_fh = opengz(kg_fn)
        id2aliases = {}
        if alias_fn:
            alias_fh = opengz(alias_fn)
            for line in alias_fh:
                if line.startswith('#'): continue
                fields = line.strip('\n').split('\t')
                if fields[0] == fields[1]:
                    continue
                if fields[0] in id2aliases:
                    id2aliases[fields[0]].append(fields[1])
                else:
                    id2aliases[fields[0]] = [fields[1]]
        else:
            err_warn("No alias file provided")
            
        self.idmap = {
            'gene_symbol': id2aliases}

        cnt = 0
        for line in kg_fh:
            if line.startswith('#'): continue
            fields = line.strip('\n').split('\t')

            t = Transcript()
            t.name = fields[0]
            t.aliases = id2aliases[t.name] # alias is the full name
            m = p_trxn_version.match(t.name)
            if m:
                t.name = m.group(1)
                t.version = int(m.group(2))

            if t.name in self.name2gene:
                g = self.name2gene[t.name]
            else:
                g = Gene(t.name)
                self.name2gene[t.name] = g

            t.chrm = normalize_chrm(fields[1])
            t.strand = fields[2]
            t.beg = int(fields[3])
            t.end = int(fields[4])
            t.cds_beg = int(fields[5])
            t.cds_end = int(fields[6])
            t.source = 'UCSC_knownGene'
            ex_begs, ex_ends = fields[8], fields[9]
            for ex_beg, ex_end in zip(list(map(int, ex_begs.strip(',').split(','))),
                                      list(map(int, ex_ends.strip(',').split(',')))):
                t.exons.append((ex_beg, ex_end))
            t.exons.sort()
            g.tpts.append(t)
            t.gene = g
            cnt += 1

        err_print("loaded %d transcripts from UCSC knownGene table." % cnt)


def main_index(args):

    """ this takes care of indexing
    1) gene/transcripts;
    2) other general features (TFBS, histone etc);
    3) reference;
    4) alias to gene/transcripts
    """
    # gene / transcripts
    if args.ensembl:
        db = EnsemblDB()
        db.index([args.ensembl])

    if args.ccds:
        db = CCDSDB()
        db.index([args.ccds])

    if args.refseq:
        db = RefSeqDB()
        db.index([args.refseq])

    if args.aceview:
        db = AceViewDB()
        db.index([args.aceview])

    if args.gencode:
        db = GENCODEDB()
        db.index([args.gencode])

    if args.kg:
        db = UCSCKnownGeneDB()
        db.index([args.kg, args.alias])

    if args.ucsc:
        db = UCSCRefGeneDB()
        db.index([args.ucsc])

    # features
    if args.gff:
        db = FeatureDB()
        db.index(args.gff, 'gff', args.sorted)

    if args.bed:
        db = FeatureDB()
        db.index(args.bed, 'bed', args.sorted)

    if args.vcf:
        db = FeatureDB()
        db.index(args.vcf, 'vcf', args.sorted)

    if args.idmap:
        if args.output is None:
            err_die("Please provide output through -o.")
            
        index_idmap(args.idmap, args.output)

    # TODO!!! let users define their own id mapping
    # if args.mapid:              # custom id mapping
    #     idxfn = args.ouptut

    # make aliases through idmaps
    if args.uniprot:
        # mapping from uniprot ID to multiple IDs
        uniprot2multi_ids = parser.parse_uniprot_mapping(args.uniprot)
        dump(uniprot2multi_ids, open(args.uniprot+'.idx','wb'), 2)

    # references
    if args.reference and args.reference != "_DEF_":
        from . import config
        config.samtools_faidx(args.reference)

def main():

    import itertools
    db = TransVarDB(sys.argv[1])
    print(list(itertools.chain(*[g.tpts for g in db.get('ZNF418')])))
    print(list(itertools.chain(*[g.tpts for g in db.get('TTK')])))
    print(list(itertools.chain(*[g.tpts for g in db.get('ENST00000574474')])))
    print(list(itertools.chain(*[g.tpts for g in db.get('TP53')])))
    print(list(itertools.chain(*[g.tpts for g in db.get('EGFR')])))
    print(list(itertools.chain(*[g.tpts for g in db.get('A2M')])))
    print(list(itertools.chain(*[g.tpts for g in db.get('NM_005228')])))
    print(db.get_by_loc('chr1',1022955, 1023955))

if __name__ == "__main__":
    main()
