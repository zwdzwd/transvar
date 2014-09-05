import transcripts as trs
import sys
from faidx import RefGenome

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
        k1 = (t.chrm, t.cds_beg / self.binsize)
        k2 = (t.chrm, t.cds_end / self.binsize)
        if k1 == k2:
            self.add_transcript_by_key(k1, t)
        else:
            self.add_transcript_by_key(k1, t)
            self.add_transcript_by_key(k2, t)

    def get_transcripts(self, chrm, pos, std_only):
        k = int(pos) / self.binsize
        if k not in key2transcripts: return []
        l = []
        for t in self.key2transcripts[k]:
            if t.cds_beg <= pos and t.cds_end >= pos:
                if std_only and t != t.gene.std_tpt:
                    continue
                l.append(t)
        return l

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

    trs.Transcript.refseq = RefGenome(args.ref)

    name2gene = {}
    if args.ensembl:
        trs.parse_ensembl_gtf(args.ensembl, name2gene)
    
    if args.refseq:
        trs.parse_refseq_gff(args.refseq, name2gene)

    if args.ccds:
        trs.parse_ccds_table(args.ccds, name2gene)

    if args.gencode:
        trs.parse_gencode_gtf(args.gencode, name2gene)

    if args.ucsc2:
        trs.parse_ucsc_refgene_customized(args.ucsc2, name2gene)

    if args.kg:
        trs.parse_ucsc_kg_table(args.kg, args.alias, name2gene)

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
    # cnt = 0
    genes = set(name2gene.values())
    for g in genes:
        # cnt += 1
        # if cnt % 10000 == 0:
        #     sys.stderr.write("\r%d genes processed\033[K" % cnt)

        for t in g.tpts:
            if not (hasattr(t, 'cds_beg') and hasattr(t, 'cds_end')):
                if t.cds:
                    t.cds.sort()
                    t.cds_beg = t.cds[0][0]
                    t.cds_end = t.cds[-1][1]
                else:
                    t.cds_beg = t.beg
                    t.cds_end = t.end

            thash.insert(t)
            if len(t) == 0:     # if no exon, use cds
                t.exons = t.cds[:]

        g.std_tpt = g.longest_tpt()

    return name2gene, thash

def parser_add_annotation(parser):

    parser.add_argument('--ref', required=True, help='indexed reference fasta (with .fai)')
    parser.add_argument('--ensembl', default=None, help='Ensembl GTF transcript annotation')
    parser.add_argument('--gencode', default=None, help='GENCODE GTF transcript annotation')
    parser.add_argument('--kg', default=None, help='UCSC knownGene transcript annotation')
    parser.add_argument('--alias', default=None, help='UCSC knownGene aliases (without providing aliases, only the knownGene id can be searched')
    parser.add_argument('--ucsc2', default=None,  help='customized UCSC transcript annotation')
    parser.add_argument('--refseq', default=None, help='RefSeq transcript annotation')
    parser.add_argument('--ccds', default=None, help='CCDS transcript annotation table')

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
