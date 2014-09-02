

def complement(base):

    return {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N'
    }[base]

def reverse_complement(seq):
    
    return ''.join([complement(base) for base in reversed(seq)])


def set_seq(seq, pos, base):

    return ''.join([base if p == pos else b for p, b in enumerate(seq)])

standard_codon_table = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', }
stop_codons  = [ 'TAA', 'TAG', 'TGA', ]
start_codons = [ 'TTG', 'CTG', 'ATG', ]

reverse_codon_table = {
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'C': ['TGT', 'TGC'],
    'E': ['GAG', 'GAA'],
    'D': ['GAC', 'GAT'],
    'G': ['GGT', 'GGG', 'GGA', 'GGC'],
    'F': ['TTT', 'TTC'],
    'I': ['ATC', 'ATA', 'ATT'],
    'H': ['CAT', 'CAC'],
    'K': ['AAG', 'AAA'],
    'M': ['ATG'],
    'L': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
    'N': ['AAC', 'AAT'],
    'Q': ['CAA', 'CAG'],
    'P': ['CCT', 'CCG', 'CCA', 'CCC'],
    'S': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'],
    'R': ['AGG', 'AGA', 'CGA', 'CGC', 'CGG', 'CGT'],
    'T': ['ACA', 'ACG', 'ACT', 'ACC'],
    'W': ['TGG'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'Y': ['TAT', 'TAC'],
    '*': ['TAA', 'TAG', 'TGA']
}

# site in codon follow the genomic order.
# no matter the strand the positive or negative, first site has
# the smallest genomic coordinate
class Codon():
    
    # chrm, locs, strand
    def __init__(self):

        self.gene   = None
        self.pos    = -1        # codon position
        self.chrm   = "NA"
        self.locs   = (-1,-1,-1)
        self.strand = "NA"
        self.seq    = "NA"      # natural sequence, not the actural sequence, can be directly mapped to amino acids
        self.index  = -1

    def __repr__(self):
        if self.locs:
            return '<Codon %s:%d,%d,%d (%s)>' % (self.chrm, self.locs[0], self.locs[1], self.locs[2], self.strand)
        else:
            return '<Codon unknown>'

    def __eq__(self, other):

        return ((self.chrm, self.strand, self.sites[0].loc) == (other.chrm, other.strand, other.sites[0].loc))

    def __hash__(self):

        return hash((self.chrm, self.sites[0].loc))

    def format(self):
        if self.locs:
            return "%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s" % (self.gene.name, self.index, self.chrm, self.locs[0], self.locs[1], self.locs[2], self.seq, self.strand)
        else:
            return "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"

class NonCoding():

    def __init__(self):

        self.gene = None
        self.region = ''
        self.closest_coding_pos = -1
        self.relative_coding_pos = 0

    def format(self):

        return "%s\t%s\t%d\t%d" % (self.gene.name, self.region, self.closest_coding_pos, self.relative_coding_pos)
        
class Transcript():

    def __init__(self, chrm, strand, start, end, seq):

        self.chrm   = chrm
        self.start  = start
        self.end    = end
        self.exons  = []
        self.seq    = seq
        self.strand = strand
        self.gene   = None
        
    def __repr__(self):
        if self.gene:
            return "<Transcript for %s: %s(%s):%d-%d>" % (self.gene.name, self.chrm, self.strand, self.start, self.end)
        else:
            return "<Empty Transcript>"

    def is_standard(self):
        return self == self.gene.std_tpt

    def cpos2codon(self, cpos):

        """ all coordinates, exons, cds are 1-based
        i.e., (200,300) means the first base is 200
        the last base is 300
        """
        cpos = int(cpos)
        if self.strand == "+":
            np = []
            for beg, end in self.exons:
                np += range(max(beg, self.cds_beg),
                            min(self.cds_end, end)+1)
            assert len(np) == len(self.seq)

            ni = cpos*3
            if ni <= len(np):
                codon        = Codon()
                codon.index  = cpos
                codon.locs   = tuple(np[ni-3:ni])
                codon.gene   = self.gene
                codon.chrm   = self.chrm
                codon.strand = self.strand
                codon.seq    = self.seq[ni-3:ni]
                return codon
            else:
                return Codon()

        else:
            np = []
            for beg, end in reversed(self.exons):
                np += range(min(self.cds_end, end),
                            max(beg, self.cds_beg)-1,-1)
            assert len(np) == len(self.seq)

            ni = cpos*3
            if ni <= len(np):
                codon        = Codon()
                codon.index  = cpos
                codon.locs   = tuple(reversed(np[ni-3:ni]))
                codon.gene   = self.gene
                codon.chrm   = self.chrm
                codon.strand = self.strand
                codon.seq    = self.seq[ni-3:ni]
                return codon
            else:
                return Codon()

    def npos2codon(self, chrm, npos):
        npos = int(npos)
        if chrm != self.chrm:
            raise Exception("Wrong chromosome.\n")
        
        if self.cds_beg > npos:
            nc = NonCoding()
            nc.region = "noncoding 5'end"
            nc.gene = self.gene
            nc.closest_coding_pos = 1
            nc.relative_coding_pos = npos - self.cds_beg
            return nc

        if self.cds_end < npos:
            nc = NonCoding()
            nc.region = "noncoding 3'end"
            nc.gene = self.gene
            nc.closest_coding_pos = len(self.seq)
            nc.relative_coding_pos = npos - self.cds_end
            return nc

        if self.strand == "+":
            np = []
            for beg, end in self.exons:
                np += range(max(beg, self.cds_beg),
                            min(self.cds_end, end)+1)
            assert len(np) == len(self.seq)

            for i, pos in enumerate(np):
                if npos == pos:
                    codon        = Codon()
                    codon.chrm   = self.chrm
                    codon.gene   = self.gene
                    codon.strand = self.strand
                    codon.index  = i/3 + 1
                    codon.seq    = self.seq[i-i%3:i-i%3+3]
                    codon.locs   = np[i-i%3:i-i%3+3]
                    codon.region = 'coding'
                    return codon
                if npos < pos:
                    nc = NonCoding()
                    nc.gene = self.gene
                    nc.region = 'intronic'
                    # print np[i-1], npos, pos
                    if npos - np[i-1] < pos - npos:
                        nc.relative_coding_pos = npos - np[i-1]
                        nc.closest_coding_pos = i # 1-based
                    else:
                        nc.relative_coding_pos = npos - pos
                        nc.closest_coding_pos = i+1 # 1-based
                    return nc
        else:
            np = []
            for beg, end in reversed(self.exons):
                np += range(min(self.cds_end, end),
                            max(beg, self.cds_beg)-1,-1)
            assert len(np) == len(self.seq)

            for i, pos in enumerate(np):
                if npos == pos:
                    codon        = Codon()
                    codon.chrm   = self.chrm
                    codon.gene   = self.gene
                    codon.strand = self.strand
                    codon.index  = i/3 + 1
                    codon.seq    = self.seq[i-i%3:i-i%3+3]
                    codon.locs   = tuple(reversed(np[i-i%3:i-i%3+3]))
                    codon.region = 'coding'
                    return codon
                if npos > pos:
                    nc = NonCoding()
                    nc.gene = self.gene
                    nc.region = 'intronic'
                    # print pos, npos, np[i-1]
                    if np[i-1] - npos < npos - pos:
                        nc.relative_coding_pos = np[i-1] - npos
                        nc.closest_coding_pos = i # 1-based
                    else:
                        nc.relative_coding_pos = pos - npos
                        nc.closest_coding_pos = i+1 # 1-based
                    return nc


class Gene():

    def __init__(self, name):

        self.name    = name
        self.tpts    = []
        self.std_tpt = None

    def __repr__(self):
        return "<Gene: %s>" % self.name

    def longest_tpt(self):

        return max(self.tpts, key=lambda x: len(x.seq))

    def chrm(self):
        
        return self.std_tpt.chrm

    def strand(self):
        
        return self.std_tpt.strand

    def cpos2codon(self, cpos):
        """ based on the longest transcript """

        return self.std_tpt.cpos2codon(cpos)

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

def parse_annotation(map_file):

    name2gene = {}
    thash = THash()
    for line in open(map_file):
        fields = line.strip().split()
        gene_name = fields[0]
        if gene_name in name2gene:
            gene = name2gene[gene_name]
        else:
            gene = Gene(gene_name)
            name2gene[gene_name] = gene

        tpt = Transcript(fields[1], fields[2],
                         int(fields[3]), int(fields[4]), fields[-1])
        tpt.cds_beg = int(fields[5])
        tpt.cds_end = int(fields[6])
        ex_begs, ex_ends = fields[8], fields[9]
        
        for ex_beg, ex_end in zip(map(int, ex_begs.split(',')),
                                  map(int, ex_ends.split(','))):
            tpt.exons.append((ex_beg, ex_end))
            
        tpt.exons.sort() # keep exons sorted
        gene.tpts.append(tpt)
        tpt.gene = gene
        thash.insert(tpt)

    for gene in name2gene.itervalues():
        gene.std_tpt = gene.longest_tpt()

    return name2gene, thash


def codondiff(c1, c2):

    diff = []
    for i in xrange(3):
        if c1[i] != c2[i]:
            diff.append(i)

    return diff

