import transcripts as trs
import sys
import faidx

class AnnoDB():
    
    def __init__(self, args):
        faidx.init_refgenome(args.reference if args.reference else None)
        self.session = None
        self.name2gene = None
        self.thash = None
        args.ffhs = {}
        if args.dbsnp:
            import pysam
            args.dbsnp_fh = pysam.Tabixfile(args.dbsnp)
        else:
            args.dbsnp_fh = None
        
        if args.sql:
            import sqlmodel
            self.sqlmodel = sqlmodel
            self.session = sqlmodel.sessionmaker(bind=sqlmodel.engine,
                                                 autoflush=False)()
        else:
            self.name2gene, self.thash = parse_annotation(args)


    def get_gene(self, name):

        if self.session:
            gs = self.session.query(self.sqlmodel.Gene).filter_by(name=name).all()
            if not gs: return None
            g = gs[0]
            gene = trs.Gene(name)
            for t in g.transcripts:
                f = t.feature
                tpt = trs.Transcript()
                tpt.chrm = f.chrm.name
                tpt.strand = '-' if t.strand == 1 else '+'
                tpt.name = t.name
                tpt.beg = f.beg
                tpt.end = f.end
                tpt.cds_beg = t.cds_beg
                tpt.cds_end = t.cds_end
                tpt.gene = gene
                tpt.source = f.source.name
                gene.tpts.append(tpt)
                for ex in t.exons:
                    tpt.exons.append((ex.beg, ex.end))
                t.exons.sort()
            return gene
        elif self.name2gene:
            if name in self.name2gene:
                return self.name2gene[name]
            else:
                return None

MAXCHRMLEN=300000000
def normalize_chrm(chrm):

    if chrm == '23': chrm = 'X'
    if chrm == '24': chrm = 'Y'
    if not chrm.startswith('chr'):
        chrm = 'chr'+chrm

    return chrm

def reflen(chrm):
    return faidx.refgenome.faidx[normalize_chrm(chrm)][0]

def printseq(seq):

    if len(seq) > 10:
        return '%s..%s' % (seq[:3], seq[-3:])
    else:
        return seq

class THash():

    def __init__(self):

        self.key2transcripts = {}
        self.binsize = 100000

    def add_transcript_by_key(self, k, t):
        if k in self.key2transcripts:
            self.key2transcripts[k].append(t)
        else:
            self.key2transcripts[k] = [t]

    def insert(self, t):
        chrm = normalize_chrm(t.chrm)
        for ki in xrange(t.cds_beg/self.binsize, t.cds_end/self.binsize+1):
            k = (chrm, ki)
            self.add_transcript_by_key(k, t)
        # k1 = (chrm, t.cds_beg / self.binsize)
        # k2 = (chrm, t.cds_end / self.binsize)
        # if k1 == k2:
        #     self.add_transcript_by_key(k1, t)
        # else:
        #     self.add_transcript_by_key(k1, t)
        #     self.add_transcript_by_key(k2, t)

    def get_transcripts_cds(self, chrm, beg, end=None, flanking=0):
        
        """ get transcript if between CDS beginning and end """
        if not end: end = beg
        chrm = normalize_chrm(chrm)
        kbeg = int(beg) / self.binsize
        kend = int(end) / self.binsize
        for ki in xrange(kbeg, kend+1):
            k = (chrm, ki)
            if k in self.key2transcripts:
                for t in self.key2transcripts[k]:
                    if t.cds_beg-flanking <= end and t.cds_end+flanking >= beg:
                        yield t

    def get_transcripts(self, chrm, beg, end=None, flanking=0):

        """ get transcript if between beginning and end """

        if not end: end = beg
        chrm = normalize_chrm(chrm)
        kbeg = int(beg) / self.binsize
        kend = int(end) / self.binsize
        for ki in xrange(kbeg, kend+1):
            k = (chrm, ki)
            # print ki, kbeg, kend, k
            if k in self.key2transcripts:
                for t in self.key2transcripts[k]:
                    if t.beg-flanking <= end and t.end+flanking >= beg:
                        yield t

    def get_closest_transcripts_upstream(self, chrm, pos):
        pos = int(pos)
        chrm = normalize_chrm(chrm)
        for ki in xrange(pos/self.binsize, -1, -1):
            k = (chrm, ki)
            if k in self.key2transcripts:
                tpts = [t for t in self.key2transcripts[k] if t.end < pos]
                if tpts:
                    tpts.sort(key=lambda t: t.end, reverse=True)
                    return tpts[0]
        return None

    def get_closest_transcripts_downstream(self, chrm, pos):
        pos = int(pos)
        chrm = normalize_chrm(chrm)
        for ki in xrange(pos/self.binsize, MAXCHRMLEN/self.binsize, 1):
            k = (chrm, ki)
            if k in self.key2transcripts:
                tpts = [t for t in self.key2transcripts[k] if t.beg > pos]
                if tpts:
                    tpts.sort(key=lambda t: t.beg)
                    return tpts[0]
        return None
    
    def get_closest_transcripts(self, chrm, beg, end):

        """ closest transcripts upstream and downstream """
        return (self.get_closest_transcripts_upstream(chrm, beg),
                self.get_closest_transcripts_downstream(chrm, end))

    # def __init__(self):
    #     # chrm => bin => list of transcripts
    #     self.chrm2bin = {}
    #     self.bins = 3000
    #     self.binsize = 100000

    # def insert(self, tpt):

    #     if tpt.chrm not in self.chrm2bin:
    #         self.chrm2bin[tpt.chrm] = [[]]*self.bins
    #     else:
    #         bin1 = self.chrm2bin[tpt.chrm][tpt.cds_beg/self.binsize]
    #         bin2 = self.chrm2bin[tpt.chrm][tpt.cds_end/self.binsize]
    #         if bin1 == bin2:
    #             bin1.append(tpt)
    #         else:
    #             bin1.append(tpt)
    #             bin2.append(tpt)

    # def get_transcripts(self,chrm, pos, std):
    #     pos = int(pos)
    #     if chrm not in self.chrm2bin:
    #         return []

    #     tpts = []
    #     for tpt in self.chrm2bin[chrm][pos / self.binsize]:
    #         if tpt.cds_beg <= pos and tpt.cds_end >= pos:
    #             if std and tpt != tpt.gene.std_tpt:
    #                 continue
    #             tpts.append(tpt)

    #     return tpts


# class THash():

#     def __init__(self):
#         self.chrm2transciprts = {}
#         self.chrm2locs = {}

#     def index(self):
#         self.chrm2transcriptsloc = {}
#         for chrm, l in self.chrm2transcripts.iteritems():
#             self.chrm2locs[chrm] = [t.beg for t in l]

#     def insert(self, t):
#         if t.chrm not in self.chrm2transciprts:
#             self.chrm2transciprts[t.chrm] = []
#         l = self.chrm2transcripts[t.chrm]
#         i = bisect_left(l, t)
#         l.insert(t, i)

#     def get_transcript(self, chrm, pos, std_only=False):
#         pos = int(pos)
#         i = bisect_left(self.chrm2locs[chrm], pos)
#         l = []
#         while True:
#             t = chrm2transcripts[i]
#             if t.beg > t: break
#             if t.end 
#             if t.end >= pos:
#                 l.append(t)
#             i += 1
#             l.append(t)

def get_config(config, option, rv=None):

    if not config.has_section(rv):
        err_print('[Warning] %s (%s) has no default value, please specify' % (option, rv))
        return None
    
    if config.has_option(rv, option):
        return config.get(rv, option)
    else:
        err_print('[Warning] %s (%s) has no default value, please specify' % (option, rv))
        return None

def replace_defaults(args, config):
    if args.refversion:
        rv = args.refversion
    elif 'refversion' in config.defaults():
        rv = config.get('DEFAULT', 'refversion')
    else:
        rv = 'hg19'

    def _set_arg_(argname, rv):
        if getattr(args, argname) == '_DEF_':
            setattr(args, argname, get_config(config, argname, rv))

    argnames = ['ensembl', 'reference', 'refseq', 'ccds',
                'gencode', 'ucsc', 'custom', 'kg', 'aceview', 'dbsnp']
    for argname in argnames:
        _set_arg_(argname, rv)

    _set_arg_('uniprot', 'idmap')

def parse_annotation(args):

    name2gene = {}
    if args.ensembl:
        trs.parse_ensembl_gtf(args.ensembl, name2gene)
    
    if args.refseq:
        trs.parse_refseq_gff(args.refseq, name2gene)

    if args.ccds:
        trs.parse_ccds_table(args.ccds, name2gene)

    if args.gencode:
        trs.parse_gencode_gtf(args.gencode, name2gene)
        try:
            import pysam
            args.ffhs['GENCODE'] = pysam.Tabixfile(args.gencode)
        except:
            err_print("Cannot import non-coding features (may need pysam).")

    if args.ucsc:
        trs.parse_ucsc_refgene(args.ucsc, name2gene)

    if args.custom:
        trs.parse_ucsc_refgene_customized(args.custom, name2gene)

    if args.kg:
        trs.parse_ucsc_kg_table(args.kg, args.alias, name2gene)

    if args.aceview:
        trs.parse_aceview_transcripts(args.aceview, name2gene)

    # remove genes without transcripts
    names_no_tpts = []
    for name, gene in name2gene.iteritems():
        # print gene, len(gene.tpts)
        if not gene.tpts:
            names_no_tpts.append(name)
    for name in names_no_tpts:
        del name2gene[name]
    sys.stderr.write('[%s] Loaded %d genes.\n' % (__name__, len(name2gene)))

    # index transcripts in a gene
    thash = THash()
    genes = set(name2gene.values())
    for g in genes:
        invalid_tpts = []
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

            thash.insert(t)
            if len(t) == 0:     # if no exon, use cds
                t.exons = t.cds[:]

        g.std_tpt = g.longest_tpt()

    if args.uniprot:
        tid2uniprot = trs.parse_uniprot_mapping(args.uniprot)
        name2protein = {}
        for name, gene in name2gene.iteritems():
            for tpt in gene.tpts:
                if tpt.name in tid2uniprot:
                    uniprot = tid2uniprot[tpt.name]
                    if uniprot not in name2protein:
                        name2protein[uniprot] = trs.Gene(uniprot)
                    prot = name2protein[uniprot]
                    prot.tpts.append(tpt)
        return name2protein, thash
    else:
        return name2gene, thash

def parser_add_annotation(parser):

    parser.add_argument('--refversion', nargs='?', default=None,
                        help='reference version (hg18, hg19, hg38 etc) (config key: refversion)')
    parser.add_argument('--reference', nargs='?', default='_DEF_',
                        help='indexed reference fasta (with .fai) (config key: reference)')
    parser.add_argument('--ensembl', nargs='?', default=None, const='_DEF_',
                        help='Ensembl GTF transcript annotation (config key: ensembl)')
    parser.add_argument('--gencode', nargs='?', default=None, const='_DEF_',
                        help='GENCODE GTF transcript annotation (config key: gencode)')
    parser.add_argument('--kg', nargs='?', default=None, const='_DEF_',
                        help='UCSC knownGene transcript annotation (config key: known_gene)')
    parser.add_argument('--alias', nargs='?', default=None, const='_DEF_',
                        help='UCSC knownGene aliases (without providing aliases, only the knownGene id can be searched (config key: known_gene_alias)')
    parser.add_argument('--ucsc', nargs='?', default=None, const='_DEF_',
                        help='UCSC transcript annotation table (config key: ucsc')
    parser.add_argument('--refseq', nargs='?', default=None, const='_DEF_',
                        help='RefSeq transcript annotation (config key: refseq)')
    parser.add_argument('--ccds', nargs='?', default=None, const='_DEF_',
                        help='CCDS transcript annotation table (config key: ccds)')
    parser.add_argument('--aceview', nargs='?', default=None, const='_DEF_',
                        help='AceView GFF transcript annotation (config key: aceview)')
    parser.add_argument('--custom', nargs='?', default=None, const='_DEF_',
                        help='A customized transcript table with sequence (config key: custom)')
    parser.add_argument('--uniprot', nargs='?', default=None, const='_DEF_',
                        help='use uniprot ID rather than gene id (config key: uniprot)')
    parser.add_argument('--dbsnp', nargs='?', default=None, const='_DEF_',
                        help='dbSNP information in annotation')
    parser.add_argument('--sql', action='store_true', help='SQL mode')


    return

class Indices:

    def __init__(self):
        self.spans = []

    def extend(self, start, end):
        self.spans.append((start, end))

    def extract(self, lst):
        result = []
        for start, end in self.spans:
            if not end:
                end = len(lst)
            result.extend([lst[_] for _ in xrange(start, end)])

        return result

def parse_indices(indstr):
    indices = Indices()
    if not indstr: return indices
    rgs = indstr.split(',')
    for rg in rgs:
        if rg.find('-') >= 0:
            pair = rg.split('-')
            if not pair[0]:
                pair[0] = 0
            if not pair[1]:
                pair[1] = None
            indices.extend(int(pair[0])-1 if pair[0] else 0,
                           int(pair[1]) if pair[1] else None)
        else:
            indices.extend(int(rg)-1, int(rg))

    return indices


def opengz(fn):
    
    if fn.endswith('.gz'):
        import gzip
        fh = gzip.open(fn)
    else:
        fh = open(fn)

    return fh

def err_die(msg, fn):

    sys.stderr.write('[%s] %s\n' % (fn, msg))
    sys.stderr.write('[%s] Abort\n' % fn)
    sys.exit(1)

def err_warn(msg, fn):
    sys.stderr.write('[%s] Warning: %s\n' % (fn, msg))

def err_raise(cls, msg, fn):
    raise cls('[%s] Exception: %s' % (fn, msg))

def err_print(msg):
    sys.stderr.write('%s\n' % str(msg))


def double_trim(seq1, seq2):

    # trim head
    head_trim = 0
    while seq1 and seq2:
        if seq1[0] == seq2[0]:
            head_trim += 1
            seq1 = seq1[1:]
            seq2 = seq2[1:]
        else:
            break

    # trim tail
    tail_trim = 0
    while seq1 and seq2:
        if seq1[-1] == seq2[-1]:
            tail_trim += 1
            seq1 = seq1[:-1]
            seq2 = seq2[:-1]
        else:
            break

    return (seq1, seq2, head_trim, tail_trim)
