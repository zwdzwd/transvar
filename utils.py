import transcripts as trs
import sys
import faidx

def normalize_chrm(chrm):

    if not chrm.startswith('chr'):
        chrm = 'chr'+chrm

    return chrm

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
        k1 = (chrm, t.cds_beg / self.binsize)
        k2 = (chrm, t.cds_end / self.binsize)
        if k1 == k2:
            self.add_transcript_by_key(k1, t)
        else:
            self.add_transcript_by_key(k1, t)
            self.add_transcript_by_key(k2, t)

    def get_transcripts(self, chrm, pos, std_only=False):
        chrm = normalize_chrm(chrm)
        k = (chrm, int(pos) / self.binsize)
        if k in self.key2transcripts:
            for t in self.key2transcripts[k]:
                if t.cds_beg <= pos and t.cds_end >= pos:
                    if std_only and t != t.gene.std_tpt:
                        continue
                    yield t

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

def parse_annotation(args):

    faidx.init_refgenome(args.ref if args.ref else None)

    name2gene = {}
    if args.ensembl:
        trs.parse_ensembl_gtf(args.ensembl, name2gene)
    
    if args.refseq:
        trs.parse_refseq_gff(args.refseq, name2gene)

    if args.ccds:
        trs.parse_ccds_table(args.ccds, name2gene)

    if args.gencode:
        trs.parse_gencode_gtf(args.gencode, name2gene)

    if args.ucsc:
        trs.parse_ucsc_refgene_customized(args.ucsc, name2gene)

    if args.kg:
        trs.parse_ucsc_kg_table(args.kg, args.alias, name2gene)

    if args.aceview:
        trs.parse_aceview_transcripts(args.aceview, name2gene)

    if args.dbsnp:
        import pysam
        args.dbsnp_fh = pysam.Tabixfile(args.dbsnp)
    else:
        args.dbsnp_fh = None

    # remove genes without transcripts
    names_no_tpts = []
    for name, gene in name2gene.iteritems():
        # print gene, len(gene.tpts)
        if not gene.tpts:
            names_no_tpts.append(name)
    for name in names_no_tpts:
        del name2gene[name]
    sys.stderr.write('[%s] Loaded %d genes.\n' % (__name__, len(name2gene)))

    # index transcripts in a genee
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

def parser_add_annotation(parser, d):

    parser.add_argument('--ref', nargs='?',
                        default=d['reference'] if 'reference' in d else None,
                        help='indexed reference fasta (with .fai) (config key: reference)')
    parser.add_argument('--ensembl', nargs='?', default=None,
                        const=d['ensembl'] if 'ensembl' in d else None,
                        help='Ensembl GTF transcript annotation (config key: ensembl)')
    parser.add_argument('--gencode', nargs='?', default=None,
                        const=d['gencode'] if 'gencode' in d else None,
                        help='GENCODE GTF transcript annotation (config key: gencode)')
    parser.add_argument('--kg', nargs='?', default=None,
                        const=d['known_gene'] if 'known_gene' in d else None,
                        help='UCSC knownGene transcript annotation (config key: known_gene)')
    parser.add_argument('--alias', nargs='?', default=None,
                        const=d['known_gene_alias'] if 'known_gene_alias' in d else None,
                        help='UCSC knownGene aliases (without providing aliases, only the knownGene id can be searched (config key: known_gene_alias)')
    parser.add_argument('--ucsc', nargs='?', default=None,
                        const=d['ucsc'] if 'ucsc' in d else None,
                        help='customized UCSC transcript annotation table (config key: ucsc')
    parser.add_argument('--refseq', nargs='?', default=None,
                        const=d['refseq'] if 'refseq' in d else None,
                        help='RefSeq transcript annotation (config key: refseq)')
    parser.add_argument('--ccds', nargs='?', default=None,
                        const=d['ccds'] if 'ccds' in d else None,
                        help='CCDS transcript annotation table (config key: ccds)')
    parser.add_argument('--aceview', nargs='?', default=None,
                        const=d['aceview'] if 'aceview' in d else None,
                        help='AceView GFF transcript annotation (config key: aceview)')
    parser.add_argument('--uniprot', nargs='?', default=None,
                        const=d['uniprot'] if 'uniprot' in d else None,
                        help='use uniprot ID rather than gene id (config key: uniprot)')
    parser.add_argument('--dbsnp', nargs='?', default=None,
                        const=d['dbsnp'] if 'dbsnp' in d else None,
                        help='dbSNP information in annotation')

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
