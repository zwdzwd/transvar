from transcripts import *

class THash():

    def __init__(self):
        # chrm => bin => list of transcripts
        self.chrm2bin = {}
        self.bins = 300
        self.binsize = 1000000

    def insert(self, tpt):

        if tpt.chrm not in self.chrm2bin:
            self.chrm2bin[tpt.chrm] = [[]]*self.bins
        else:
            bin1 = self.chrm2bin[tpt.chrm][tpt.cds_beg/self.binsize]
            bin2 = self.chrm2bin[tpt.chrm][tpt.cds_end/self.binsize]
            if bin1 == bin2:
                bin1.append(tpt)
            else:
                bin1.append(tpt)
                bin2.append(tpt)

    def get_transcripts(self,chrm, pos, std):
        pos = int(pos)
        if chrm not in self.chrm2bin:
            return []

        tpts = []
        for tpt in self.chrm2bin[chrm][pos / self.binsize]:
            if tpt.cds_beg <= pos and tpt.cds_end >= pos:
                if std and tpt != tpt.gene.std_tpt:
                    continue
                tpts.append(tpt)

        return tpts

def parse_annotation(args):

    name2gene = {}
    if args.ensembl:
        parse_ensembl_gtf(args.ensembl, name2gene, thash)
    
    if args.refseq:
        parse_refseq_gff(args.refseq, name2gene, thash)

    if args.ccds:
        parse_ccds_table(args.ccds, name2gene, thash)

    if args.gencode:
        parse_gencode_gtf(args.gencode, name2gene, thash)

    if args.ucsc2:
        parse_ucsc_refgene_customized(args.ucsc2, name2gene, thash)

    if args.kg:
        parse_ucsc_kg_table(args.kg, args.alias, name2gene, thash)

    # establish reverse link
    # and the primary transcript
    thash = THash()
    names_no_tpts = []
    for name, gene in name2gene.iteritems():
        if not gene.tpts:
            names_no_tpts.append(name)
            continue
        for t in gene.tpts:
            if not (t.cds_beg and t.cds_end):
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

        gene.std_tpt = gene.longest_tpt()

    # remove genes without transcripts
    for name in names_no_tpts:
        del name2gene[name]

    sys.stderr.write("Loaded %d genes with transcripts.\n" % len(name2gene))

    return name2gene, thash

def parser_add_annotation(parser):

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
