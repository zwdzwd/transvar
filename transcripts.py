import sys, re
from utils import *
import faidx

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
    'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*'}
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
        self.chrm   = "NA"
        self.locs   = (-1,-1,-1)
        self.strand = "NA"
        self.seq    = '' # natural sequence, not the actural sequence, can be directly mapped to amino acids
        self.index  = -1

    def refseq(self):

        if self.strand == '+': return self.seq
        else: return reverse_complement(self.seq)

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

def reverse_tnuc_pos(codon, tnuc_pos):
    
    if codon.strand == '+':
        return codon.locs[(tnuc_pos-1)%3]
    else:
        return codon.locs[2-(tnuc_pos-1)%3]

def codondiff(c1, c2):

    diff = []
    for i in xrange(3):
        if c1[i] != c2[i]:
            diff.append(i)

    return diff

class NonCoding():

    def __init__(self):

        self.gene = None
        self.region = ''
        self.closest_coding_pos = -1
        self.relative_coding_pos = 0

    def format(self):

        return "%s\t%s\t%d\t%d" % (self.gene.name, self.region, self.closest_coding_pos, self.relative_coding_pos)
        
class Transcript():

    def __init__(self):

        """ chrm, strand, start, end, seq (optional), cds_beg, cds_end """

        self.gene   = None
        self.seq    = None
        self.name   = '.'
        self.exons  = []
        self.cds    = []

    def __len__(self):
        if self.seq:
            return len(self.seq)
        else:
            return reduce(lambda x,y: x+y,
                          [end-beg+1 for beg, end in self.exons], 0)

    def ensure_seq(self):
        """ return True when successful """
        if self.seq: return True
        if not faidx.refgenome:
            err_die("Please provide reference through --ref.\n")
        seq = faidx.refgenome.fetch_sequence(self.chrm, self.beg, self.end)

        if not seq: return False
        segs = []
        for ex_beg, ex_end in self.exons:
            beg = max(ex_beg, self.cds_beg)
            end = min(ex_end, self.cds_end)
            if beg <= end:
                segs.append(seq[beg-self.beg:end+1-self.beg])
        self.seq = ''.join(segs)
        if self.strand == '-':
            self.seq = reverse_complement(self.seq)
        return True

    def __repr__(self):
        if self.gene:
            return "<Transcript %s %s: %s(%s):%d-%d>" % (self.name, self.gene.name, self.chrm, self.strand, self.beg, self.end)
        else:
            return "<Empty Transcript>"

    def is_standard(self):
        return self == self.gene.std_tpt

    def tnuc_range2gnuc_range(self, tbeg, tend):

        """ convert transcript range to genomic range
        tbeg and tend are 1-based
        """
        if not self.ensure_seq(): return None
        if self.strand == "+":
            np = []
            for beg, end in self.exons:
                np += range(max(beg, self.cds_beg),
                            min(self.cds_end, end)+1)
            assert len(np) == len(self.seq)
            return np[tbeg-1], np[tend-1]
        else:
            np = []
            for beg, end in reversed(self.exons):
                np += range(min(self.cds_end, end),
                            max(beg, self.cds_beg)-1,-1)
            assert len(np) == len(self.seq)
            return np[tend-1], np[tbeg-1]

    def cpos2codon(self, cpos):

        """ all coordinates, exons, cds are 1-based
        i.e., (200,300) means the first base is 200
        the last base is 300
        """
        if not self.ensure_seq(): return None
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
                return None
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
                return None

    def npos2codon(self, chrm, npos):
        if not self.ensure_seq(): return None
        npos = int(npos)

        # no check chrm == self.chrm, due to differential
        # naming convention: chr12 vs 12.
        
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
        self.pseudo  = False

    def __repr__(self):
        return "<Gene: %s>" % self.name

    def longest_tpt(self):

        return max(self.tpts, key=lambda x: len(x))

    def chrm(self):
        
        return self.std_tpt.chrm

    def strand(self):
        
        return self.std_tpt.strand

    def cpos2codon(self, cpos):
        """ based on the longest transcript """

        return self.std_tpt.cpos2codon(cpos)

def parse_ucsc_refgene_customized(map_file, name2gene):

    cnt = 0
    for line in open(map_file):
        fields = line.strip().split()
        gene_name = fields[0]
        if gene_name in name2gene:
            gene = name2gene[gene_name]
        else:
            gene = Gene(gene_name)
            name2gene[gene_name] = gene

        t = Transcript()
        t.chrm = fields[1]
        t.strand = fields[2]
        t.beg    = int(fields[3])
        t.end    = int(fields[4])
        t.seq    = fields[-1]
        t.cds_beg = int(fields[5])
        t.cds_end = int(fields[6])
        t.source = 'UCSC_refGene'
        t.name = '.'
        ex_begs, ex_ends = fields[8], fields[9]

        for ex_beg, ex_end in zip(map(int, ex_begs.split(',')),
                                  map(int, ex_ends.split(','))):
            t.exons.append((ex_beg, ex_end))
            
        t.exons.sort() # keep exons sorted
        gene.tpts.append(t)
        t.gene = gene
        cnt += 1

    sys.stderr.write('[%s] Loaded %d transcripts from UCSC refgene (customized).\n' % (__name__, cnt))

    return

class Region():

    def __init__(self, name, beg, end):

        self.name = name
        self.beg = beg
        self.end = end

def parse_refseq_gff(gff_fn, name2gene):

    id2ent = {}
    gff_fh = opengz(gff_fn)
    reg = None
    cnt = 0
    for line in gff_fh:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        # print line.strip()
        info = dict([_.split('=') for _ in fields[8].split(';')])
        if fields[2] == 'region':
            if 'Name' in info:
                reg = Region(info['Name'], int(fields[3]), int(fields[4]))
            else:
                reg = None
        elif (reg and fields[2] == 'gene' and
              ('pseudo' not in info or info['pseudo'] != 'true')):
            gene_name = info['Name']
            if gene_name in name2gene:
                g = name2gene[gene_name]
            else:
                g = Gene(gene_name)
                name2gene[gene_name] = g
            g.beg = int(fields[3])
            g.end = int(fields[4])
            id2ent[info['ID']] = g
        elif fields[2] == 'mRNA' and info['Parent'] in id2ent:
            t = Transcript()
            t.chrm = reg.name
            t.strand = fields[6]
            t.beg = int(fields[3])
            t.end = int(fields[4])
            t.name = info['Name']
            t.gene = id2ent[info['Parent']]
            t.gene.tpts.append(t)
            t.source = 'RefSeq'
            id2ent[info['ID']] = t
            cnt += 1
        elif fields[2] == 'exon' and info['Parent'] in id2ent:
            t = id2ent[info['Parent']]
            if (isinstance(t, Gene)):
                g = t
                if not hasattr(g, 'gene_t'):
                    g.gene_t = Transcript()
                    g.tpts.append(g.gene_t)
                    g.gene_t.chrm = reg.name
                    g.gene_t.strand = fields[6]
                    g.gene_t.gene = g
                    g.gene_t.beg = g.beg
                    g.gene_t.end = g.end
                    g.gene_t.source = 'RefSeq'
                    cnt += 1
                t = g.gene_t
            t.exons.append((int(fields[3]), int(fields[4])))
        elif fields[2] == 'CDS' and info['Parent'] in id2ent:
            t = id2ent[info['Parent']]
            if (isinstance(t, Gene)):
                g = t
                if not hasattr(g, 'gene_t'):
                    g.gene_t = Transcript()
                    g.tpts.append(g.gene_t)
                    g.gene_t.chrm = reg.name
                    g.gene_t.strand = fields[6]
                    g.gene_t.gene = g
                    g.gene_t.beg = g.beg
                    g.gene_t.end = g.end
                    g.gene_t.source = 'RefSeq'
                    cnt += 1
                t = g.gene_t
            t.cds.append((int(fields[3]), int(fields[4])))

    sys.stderr.write("[%s] Loaded %d transcripts from RefSeq GFF3 file.\n" % (__name__, cnt))

def parse_ensembl_gtf(gtf_fn, name2gene):
    """ gtf file is gffv2 """

    gtf_fh = opengz(gtf_fn)
    id2ent = {}
    cnt = 0
    for line in gtf_fh:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        info = dict(re.findall(r'\s*([^"]*) "([^"]*)";', fields[8]))
        # info = dict([_.split('=') for _ in fields[8].split(';')])
        if fields[2] == 'gene' and info['gene_biotype'] == 'protein_coding':
            gene_name = info['gene_name']
            if gene_name in name2gene:
                g = name2gene[gene_name]
            else:
                g = Gene(gene_name)
                name2gene[gene_name] = g
            g.beg = int(fields[3])
            g.end = int(fields[4])
            id2ent[info['gene_id']] = g
        elif fields[2] == 'transcript' and info['gene_biotype'] == 'protein_coding':
            t = Transcript()
            t.chrm = fields[0]
            t.strand = fields[6]
            t.beg = int(fields[3])
            t.end = int(fields[4])
            t.name = info['transcript_id']
            t.gene = id2ent[info['gene_id']]
            t.gene.tpts.append(t)
            t.source = 'Ensembl'
            id2ent[t.name] = t
            cnt += 1
        elif fields[2] == 'exon' and info['gene_biotype'] == 'protein_coding':
            t = id2ent[info['transcript_id']]
            t.exons.append((int(fields[3]), int(fields[4])))
        elif fields[2] == 'CDS' and info['gene_biotype'] == 'protein_coding':
            t = id2ent[info['transcript_id']]
            t.cds.append((int(fields[3]), int(fields[4])))

    sys.stderr.write("[%s] Loaded %d transcripts from Ensembl GTF file.\n" % (__name__, cnt))

def parse_ccds_table(ccds_fn, name2gene):

    ccds_fh = open(ccds_fn)
    ccds_fh.readline()
    cnt = 0
    for line in ccds_fh:
        fields = line.strip().split('\t')
        if fields[5] != 'Public':
            continue
        gene_name = fields[2]
        if gene_name not in name2gene:
            name2gene[gene_name] = Gene(fields[2])

        g = name2gene[gene_name]
        t = Transcript()
        t.chrm = fields[0]
        t.strand = fields[6]
        t.cds_beg = int(fields[7])+1
        t.cds_end = int(fields[8])+1

        # without UTR information, take CDS boundary as the exon boundary
        t.beg = t.cds_beg
        t.end = t.cds_end

        t.name = fields[4]
        t.cds = [(int(b)+1, int(e)+1) for b,e in re.findall(r"[\s\[](\d+)-(\d+)[,\]]", fields[9])]
        t.source = 'CDDS'
        t.gene = g
        g.tpts.append(t)
        cnt += 1

    sys.stderr.write("[%s] Loaded %d transcripts from CCDS table.\n" % (__name__, cnt))

def parse_ucsc_kg_table(kg_fn, alias_fn, name2gene):

    kg_fh = opengz(kg_fn)
    id2aliases = {}
    if alias_fn:
        alias_fh = opengz(alias_fn)
        for line in alias_fh:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if fields[0] in id2aliases:
                id2aliases[fields[0]].append(fields[1])
            else:
                id2aliases[fields[0]] = [fields[1]]

    cnt = 0
    for line in kg_fh:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        g = None
        if fields[0] in id2aliases:
            for alias in id2aliases[fields[0]]:
                if alias in name2gene:
                    g = name2gene[alias]
            if not g:
                g = Gene(fields[0])
            for alias in id2aliases[fields[0]]:
                name2gene[alias] = g
        else:
            if fields[0] in name2gene:
                g = name2gene[fields[0]]
            else:
                g = Gene(fields[0])
            name2gene[fields[0]] = g

        t = Transcript()
        t.name = fields[0]
        t.chrm = fields[1]
        t.strand = fields[2]
        t.beg = int(fields[3])
        t.end = int(fields[4])
        t.cds_beg = int(fields[5])
        t.cds_end = int(fields[6])
        t.source = 'UCSC_knownGene'
        ex_begs, ex_ends = fields[8], fields[9]
        for ex_beg, ex_end in zip(map(int, ex_begs.strip(',').split(',')),
                                  map(int, ex_ends.strip(',').split(','))):
            t.exons.append((ex_beg, ex_end))
        t.exons.sort()
        g.tpts.append(t)
        t.gene = g
        cnt += 1

    sys.stderr.write("[%s] Loaded %d transcripts from UCSC knownGene table.\n" % (__name__, cnt))

def parse_gencode_gtf(gencode_fn, name2gene):

    id2ent = {}
    gencode_fh = opengz(gencode_fn)
    cnt = 0
    for line in gencode_fh:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        info = dict(re.findall(r'\s*([^"]*) "([^"]*)";', fields[8]))
        if fields[2] == 'gene':
            gene_name = info['gene_name']
            if gene_name in name2gene:
                g = name2gene[gene_name]
            else:
                g = Gene(gene_name)
                name2gene[gene_name] = g
            g.beg = int(fields[3])
            g.end = int(fields[4])
            id2ent[info['gene_id']] = g
            if info['gene_type'] == 'pseudogene':
                g.pseudo = True
        elif fields[2] == 'transcript' and info['gene_type'] == 'protein_coding':
            t = Transcript()
            t.chrm = fields[0]
            t.strand = fields[6]
            t.beg = int(fields[3])
            t.end = int(fields[4])
            t.name = info['transcript_id']
            t.gene = id2ent[info['gene_id']]
            t.gene.tpts.append(t)
            t.source = 'GENCODE'
            id2ent[t.name] = t
            cnt += 1
        elif fields[2] == 'exon' and info['gene_type'] == 'protein_coding':
            t = id2ent[info['transcript_id']]
            t.exons.append((int(fields[3]), int(fields[4])))
        elif fields[2] == 'CDS' and info['gene_type'] == 'protein_coding':
            t = id2ent[info['transcript_id']]
            t.cds.append((int(fields[3]), int(fields[4])))

    sys.stderr.write("[%s] Loaded %d transcripts from GENCODE GTF file.\n" % (__name__, cnt))

def parse_aceview_transcripts(aceview_gff_fn, name2gene):

    id2tpt = {}
    aceview_fh = opengz(aceview_gff_fn)
    for line in aceview_fh:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        info = dict(re.findall(r'\s*(\S+) (\S+);', fields[8]))
        if fields[2] == 'CDS':
            gene_name = info['gene_id']
            if gene_name in name2gene:
                g = name2gene[gene_name]
            else:
                g = Gene(gene_name)
                name2gene[gene_name] = g

            if info['transcript_id'] in id2tpt:
                t = id2tpt[info['transcript_id']]
            else:
                t = Transcript()
                t.chrm = fields[0]
                t.strand = fields[6]
                t.name = info['transcript_id']
                id2tpt[t.name] = t
                t.gene = g
                g.tpts.append(t)
                t.source = 'AceView'

            t.cds.append((int(fields[3]), int(fields[4])))

        elif fields[2] == 'exon':
            gene_name = info['gene_id']
            if gene_name in name2gene:
                g = name2gene[gene_name]
            else:
                g = Gene(gene_name)
                name2gene[gene_name] = g

            if info['transcript_id'] in id2tpt:
                t = id2tpt[info['transcript_id']]
            else:
                t = Transcript()
                t.chrm = fields[0]
                t.strand = fields[6]
                t.name = info['transcript_id']
                id2tpt[t.name] = t
                t.gene = g
                g.tpts.append(t)
                t.source = 'AceView'

            t.exons.append((int(fields[3]), int(fields[4])))

    # skip transcripts without CDS, e.g., LOC391566.aAug10-unspliced
    for tid, t in id2tpt.iteritems():
        if t.cds and t.exons:
            t.exons.sort()
            t.beg = t.exons[0][0]
            t.end = t.exons[-1][1]
        else:
            t.gene.tpts.remove(t)

    sys.stderr.write("[%s] Loaded %d transcripts from AceView GFF file.\n" % (__name__, len(id2tpt)))


def extend_taa_seq(taa_pos_base, old_seq, new_seq):

    taa_pos = None
    termlen = None
    for i in xrange(len(new_seq)/3):
        taa_ref_run = standard_codon_table[old_seq[3*i:3*i+3]]
        taa_alt_run = standard_codon_table[new_seq[3*i:3*i+3]]
        # print i, old_seq[3*i:3*i+3], new_seq[3*i:3*i+3], taa_ref_run, taa_alt_run, taa_pos
        if taa_pos == None and taa_ref_run != taa_alt_run:
            taa_pos = i
            taa_ref = taa_ref_run
            taa_alt = taa_alt_run
        if taa_alt_run == '*':
            if taa_pos == None:
                err_die('Terminating codon encountered before difference.', __name__)
                return None
            termlen = i + 1 - taa_pos
            break
    if termlen == None:
        err_die('No terminating codon before the end of the new transcript.', __name__)
        return None
    taa_pos += taa_pos_base

    return taa_pos, taa_ref, taa_alt, termlen

