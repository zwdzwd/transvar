import sys, re
from utils import *
from err import *
import faidx
from record import *

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

def codon2aa(codonseq):
    if codonseq not in standard_codon_table:
        raise IncompatibleTranscriptError('Invalid codon sequence')
    return standard_codon_table[codonseq]

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

def aa2codon(aa):
    if aa not in reverse_codon_table:
        raise IncompatibleTranscriptError('Invalid amino acid')
    return reverse_codon_table[aa]

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

    def aa(self):
        return codon2aa(self.seq)

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

def tnuc2gnuc(np, tnuc_pos):
    """ np is the position array
    take integer as input
    """
    if tnuc_pos >= len(np):
        raise IncompatibleTranscriptError()
    return np[tnuc_pos-1]

def tnuc2gnuc2(np, tnuc_pos, tpt):
    """ take Pos as input """
    if tpt.strand == '-':
        return tnuc2gnuc(np, tnuc_pos.pos) - tnuc_pos.tpos
    else:
        return tnuc2gnuc(np, tnuc_pos.pos) + tnuc_pos.tpos

def check_exon_boundary(np, pos):
    """ check consistency with exon boundary """

    if pos.tpos > 0:
        if abs(tnuc2gnuc(np, pos.pos) - tnuc2gnuc(np, pos.pos+1)) == 1:
            raise IncompatibleTranscriptError()
    elif pos.tpos < 0:
        if abs(tnuc2gnuc(np, pos.pos) - tnuc2gnuc(np, pos.pos-1)) == 1:
            raise IncompatibleTranscriptError()

def region_in_exon(np, beg, end):
    """ region in tnuc positions """

    if beg.tpos != 0: return False
    if end.tpos != 0: return False
    for i in xrange(beg.pos, end.pos-1):
        if abs(np[i] - np[i+1]) != 1:
            return False
    return True

def region_in_intron(np, beg, end):
    """ region in tnuc positions """

    if beg.tpos == 0 or end.tpos == 0: return False
    if beg.pos == end.pos and beg.tpos*end.tpos > 0:
        return True
    if beg.pos+1 == end.pos and beg.tpos>0 and end.tpos<0:
        return True
    if end.pos+1 == beg.pos and beg.tpos<0 and end.tpos>0:
        return True
    return False

def tnuc_range2gnuc_range_(np, tbeg, tend):

    """ convert transcript range to genomic range
    tbeg and tend are 1-based
    """
    return min(np[tbeg-1], np[tend-1]), max(np[tbeg-1], np[tend-1])

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

    def region(self, gnuc_beg, gnuc_end):
        """ annotate genomic region with respect to this transcript """
        # check if gnuc_beg and gnuc_end are inside the genomic region
        pexon = None
        overlapping_exons = []
        for exon in self.exons:
            if (exon[0] <= gnuc_beg and exon[1] >= gnuc_end):
                _cds_beg = min(self.cds_beg, self.cds_end)
                _cds_end = max(self.cds_beg, self.cds_end)
                
                if gnuc_beg > _cds_beg and gnuc_end < _cds_end:
                    return 'Coding'
                elif gnuc_beg < _cds_beg and gnuc_end < _cds_beg:
                    return "5'UTR" if self.strand == '+' else "3'UTR"
                elif gnuc_beg > _cds_end and gnuc_end > _cds_end:
                    return "3'UTR" if self.strand == '+' else "5'UTR"
                elif gnuc_beg < _cds_beg:
                    return "5'UTR;coding" if self.strand == '+' else "3'UTR;coding"
                elif gnuc_end > _cds_end:
                    return "coding;3'UTR" if self.strand == '+' else "coding;5'UTR"
                else:
                    return "Unknown"
            if exon[0] >= gnuc_beg and exon[0] <= gnuc_end:
                overlapping_exons.append(exon)
            if pexon and gnuc_beg > pexon[1] and gnuc_end < exon[0]:
                return 'Intronic'
            pexon = exon

        if overlapping_exons:
            return 'Intronic;Exonic'
        else:
            return 'Unknown'

    def ensure_seq(self):
        """ return True when successful,
        potential reason include patch chromosomes
        """
        if self.seq: return
        if not faidx.refgenome:
            err_die("Please provide reference through --ref [reference fasta].", __name__)
        seq = faidx.refgenome.fetch_sequence(self.chrm, self.beg, self.end)

        if not seq: raise SequenceRetrievalError()
        segs = []
        for ex_beg, ex_end in self.exons:
            beg = max(ex_beg, self.cds_beg)
            end = min(ex_end, self.cds_end)
            if beg <= end:
                segs.append(seq[beg-self.beg:end+1-self.beg])
        self.seq = ''.join(segs)
        if self.strand == '-':
            self.seq = reverse_complement(self.seq)
        return

    def __repr__(self):
        if self.gene:
            return "<Transcript %s %s: %s(%s):%d-%d>" % (self.name, self.gene.name, self.chrm, self.strand, self.beg, self.end)
        else:
            return "<Empty Transcript>"

    def is_standard(self):
        return self == self.gene.std_tpt

    def position_array(self):
        if self.strand == "+":
            np = []
            for beg, end in self.exons:
                np += range(max(beg, self.cds_beg),
                            min(self.cds_end, end)+1)
        else:
            np = []
            for beg, end in reversed(self.exons):
                np += range(min(self.cds_end, end),
                            max(beg, self.cds_beg)-1,-1)

        return np

    def tnuc_range2gnuc_range(self, tbeg, tend):

        """ convert transcript range to genomic range
        tbeg and tend are 1-based
        """
        np = self.position_array()
        return tnuc_range2gnuc_range_(np, tbeg, tend)

    def taa2aa(self, taa):
        self.ensure_seq()
        if taa*3 > len(self.seq):
            raise IncompatibleTranscriptError('Incompatible reference amino acid')
        return codon2aa(self.seq[taa*3-3:taa*3])

    def taa_range2aa_seq(self, taa_beg, taa_end):

        if taa_beg*3 > len(self) or taa_end*3 > len(self):
            raise IncompatibleTranscriptError('codon nonexistent')

        return translate_seq(self.seq[taa_beg*3-3:taa_end*3])

    def tnuc2codon_(self, tnuc_pos):
        taa_pos = (tnuc_pos + 2) / 3
        codon = self.cpos2codon(taa_pos)
        pos_r = (tnuc_pos-1) % 3 # 0,1,2 for first, second and third base
        if self.strand == '-': pos_r = 2 - pos_r
        return codon, pos_r, codon.locs[pos_r]

    def _tnuc_range2exon_inds(self, tnuc_beg, tnuc_end):

        exoninds = []
        if self.strand == '+':
            for i, (beg, end) in enumerate(self.exons):
                exoninds.extend([i+1]*(min(self.cds_end, end)-max(beg, self.cds_beg)+1))
        else:
            for i, (beg, end) in enumerate(reversed(self.exons)):
                exoninds.extend([i+1]*(min(self.cds_end, end)-max(beg, self.cds_beg)+1))

        return sorted(list(set(exoninds[tnuc_beg-1:tnuc_end])))


    def tnuc_range2exon_inds(self, tnuc_beg, tnuc_end):

        return ';'.join(map(str, self._tnuc_range2exon_inds(tnuc_beg, tnuc_end)))

    def cpos2codon(self, cpos):

        """ all coordinates, exons, cds are 1-based
        i.e., (200,300) means the first base is 200
        the last base is 300
        cpos is taa_pos
        """
        self.ensure_seq()
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
                raise IncompatibleTranscriptError()
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
                raise IncompatibleTranscriptError()

    def _init_codon_(self, index):
        c = Codon()
        c.chrm = self.chrm
        c.gene = self.gene
        c.strand = self.strand
        c.index = index
        return c

    # def _gpos2codon_UTR(self, gpos, np):
    #     """ UTR region """
    #     if self.cds_beg > gpos:
    #         p = Pos(1, gpos-self.cds_beg)
    #         c = self._init_codon_(1)
    #         c.seq = self.seq[:3]
    #         c.locs = np[:3]
    #         reg = '5-UTR' if self.strand == '+' else '3-UTR'
    #         return c, p, reg

    #     if self.cds_end < gpos:
    #         p = Pos(len(self.seq), gpos-self.cds_end)
    #         c = self._init_codon_((len(self.seq)+2)/3)
    #         c.seq = self.seq[c.index*3-3:c.index*3]
    #         c.locs = np[c.index*3-3:c.index*3]
    #         reg = '3-UTR' if self.strand == '+' else '5-UTR'
    #         return c, p, reg
    #     return None

    def _in_exon(self, gpos):

        for i, exon in enumerate(self.exons):
            exind = i+1 if self.strand == '+' else len(self.exons) - i
            if i == 0 and exon[0] > gpos:
                return 'IntergenicUpstream'
            if i == len(self.exons)-1 and exon[1] < gpos:
                return 'IntergenicDownstream'
            if exon[0] <= gpos and exon[1] >= gpos:
                if gpos >= self.cds_beg and gpos <= self.cds_end:
                    return 'Coding_%d' % exind
                else:
                    return 'Exonic_%d' % exind
            if i > 0:
                pexon = self.exons[i-1]
                if gpos > pexon[1] and gpos < exon[0]:
                    if self.strand == '+':
                        return 'Intronic_%d_%d' % (exind-1, exind)
                    else:
                        return 'Intronic_%d_%d' % (exind, exind+1)

        raise Exception()       # you shouldn't reach here
    
    def _gpos2codon_p(self, gpos, np):

        if gpos < self.cds_beg:
            p = Pos(1, gpos-self.cds_beg)
            c = self._init_codon_(1)
            c.seq = self.seq[:3]
            c.locs = np[:3]
            reg = '5-UTR;'+self._in_exon(gpos)
            return c, p, reg

        if gpos > self.cds_end:
            p = Pos(len(self.seq), gpos-self.cds_end)
            c = self._init_codon_((len(self.seq)+2)/3)
            c.seq = self.seq[c.index*3-3:c.index*3]
            c.locs = np[c.index*3-3:c.index*3]
            reg = '3-UTR;'+self._in_exon(gpos)
            return c, p, reg
        
        for i, pos in enumerate(np):
            if gpos == pos:
                c = self._init_codon_(i/3+1)
                c.seq    = self.seq[i-i%3:i-i%3+3]
                c.locs   = np[i-i%3:i-i%3+3]
                p = Pos(i+1, 0)
                return c, p, self._in_exon(gpos)
            if gpos < pos:
                if gpos - np[i-1] < pos - gpos:
                    p = Pos(i, gpos-np[i-1])
                    ci = i/3+1
                else:
                    p = Pos(i+1, gpos-pos)
                    ci = (i+1)/3+1
                c = self._init_codon_(ci)
                c.seq = self.seq[ci*3-3:ci*3]
                c.locs = np[ci*3-3:ci*3]
                return c, p, self._in_exon(gpos)

    def _gpos2codon_n(self, gpos, np):

        if gpos < self.cds_beg:
            p = Pos(1, gpos-self.cds_beg)
            c = self._init_codon_(1)
            c.seq = self.seq[:3]
            c.locs = np[:3]
            reg = '3-UTR;'+self._in_exon(gpos)
            return c, p, reg

        if gpos > self.cds_end:
            p = Pos(len(self.seq), gpos-self.cds_end)
            c = self._init_codon_((len(self.seq)+2)/3)
            c.seq = self.seq[c.index*3-3:c.index*3]
            c.locs = np[c.index*3-3:c.index*3]
            reg = '5-UTR;'+self._in_exon(gpos)
            return c, p, reg

        for i, pos in enumerate(np):
            if gpos == pos:
                c = self._init_codon_(i/3+1)
                c.seq = self.seq[i-i%3:i-i%3+3]
                c.locs = tuple(reversed(np[i-i%3:i-i%3+3]))
                p = Pos(i+1, 0)
                return c, p, self._in_exon(gpos)
            
            if gpos > pos:
                if np[i-1] - gpos < gpos - pos:
                    p = Pos(i, np[i-1] - gpos)
                    ci = i/3+1
                else:
                    p = Pos(i+1, pos - gpos)
                    ci = (i+1)/3+1
                c = self._init_codon_(ci)
                c.seq = self.seq[ci*3-3:ci*3]
                c.locs = np[ci*3-3:ci*3]
                return c, p, self._in_exon(gpos)

    
    def gpos2codon(self, chrm, gpos):
        """ return Codon as well as (tnuc) Pos """
        self.ensure_seq()
        gpos = int(gpos)

        # no check chrm == self.chrm, due to differential
        # naming convention: chr12 vs 12.
        np = self.position_array()
        assert len(np) == len(self.seq)

        # ret = self._gpos2codon_UTR(gpos, np)
        # if ret: return ret
        
        if self.strand == "+":
            return self._gpos2codon_p(gpos, np)
        else:
            return self._gpos2codon_n(gpos, np)

    def overlap_region(self, beg, end):

        if self.beg >= beg and self.end <= end:
            return 'whole'

        coding = False
        intronic = False
        UTR5 = False
        UTR3 = False
        p_ex_end = None
        for ex_beg, ex_end in self.exons:
            if ex_end >= beg and ex_beg <= end:
                if beg < self.cds_beg:
                    if self.strand == '+': UTR5 = True
                    else: UTR3 = True
                if end > self.cds_beg: coding = True
                if end > self.cds_end:
                    if self.strand == '+': UTR3 = True
                    else: UTR5 = True
                if beg < self.cds_end: coding = True
            if p_ex_end and p_ex_end < end and beg < ex_beg:
                # p_ex_end---ex_beg vs beg---end
                intronic = True
            p_ex_end = ex_end

        regc = []
        if self.strand == '+':
            if UTR5: regc.append('UTR5')
            if coding: regc.append('coding')
            if intronic: regc.append('intronic')
            if UTR3: regc.append('UTR3')
        else:
            if UTR3: regc.append('UTR3')
            if coding: regc.append('coding')
            if intronic: regc.append('intronic')
            if UTR5: regc.append('UTR5')

        return ','.join(regc)


class Gene():

    def __init__(self, name=''):

        self.name    = name
        self.alias   = ''
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


def parse_ucsc_refgene(map_file, name2gene):
    """ start 1-based, end 1-based """

    cnt = 0
    for line in opengz(map_file):
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        if fields[13] != 'cmpl' or fields[14] != 'cmpl':
            continue
        gene_name = fields[12].upper()
        if gene_name in name2gene:
            gene = name2gene[gene_name]
        else:
            gene = Gene(gene_name)
            name2gene[gene_name] = gene
        t = Transcript()
        t.name = fields[1]
        t.chrm = fields[2]
        t.strand = fields[3]
        t.beg    = int(fields[4])+1
        t.end    = int(fields[5])
        t.cds_beg = int(fields[6])+1
        t.cds_end = int(fields[7])
        t.source = 'UCSC_refGene'
        ex_begs, ex_ends = fields[9], fields[10]

        for ex_beg, ex_end in zip(map(lambda x: int(x)+1, ex_begs.strip(',').split(',')),
                                  map(int, ex_ends.strip(',').split(','))):
            t.exons.append((ex_beg, ex_end))
            
        t.exons.sort() # keep exons sorted
        gene.tpts.append(t)
        t.gene = gene
        cnt += 1

    sys.stderr.write('[%s] Loaded %d transcripts from UCSC refgene.\n' % (__name__, cnt))

    return


def parse_ucsc_refgene_customized(map_file, name2gene):

    """ start 1-based, end 1-based """
    cnt = 0
    for line in open(map_file):
        fields = line.strip().split()
        gene_name = fields[0].upper()
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
        t.source = 'custom'
        t.name = '.'
        ex_begs, ex_ends = fields[8], fields[9]

        for ex_beg, ex_end in zip(map(int, ex_begs.split(',')),
                                  map(int, ex_ends.split(','))):
            t.exons.append((ex_beg, ex_end))
            
        t.exons.sort() # keep exons sorted
        gene.tpts.append(t)
        t.gene = gene
        cnt += 1

    sys.stderr.write('[%s] Loaded %d transcripts from customized table.\n' % (__name__, cnt))

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
            if 'chromosome' in info:
                reg = Region(info['chromosome'], int(fields[3]), int(fields[4]))
            # else:
            # reg = None
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
    """ gtf file is gffv2
    parser does not assume order in the GTF file
    """

    gtf_fh = opengz(gtf_fn)
    id2ent = {}
    cnt = 0
    for line in gtf_fh:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        info = dict(re.findall(r'\s*([^"]*) "([^"]*)";', fields[8]))
        # info = dict([_.split('=') for _ in fields[8].split(';')])
        if fields[2] == 'gene' and info['gene_biotype'] == 'protein_coding':
            gene_id = info['gene_id']
            if gene_id not in id2ent: id2ent[gene_id] = Gene()
            g = id2ent[gene_id]
            g.name = info['gene_name'].upper()
            if g.name not in name2gene: name2gene[g.name] = g
            g.beg = int(fields[3])
            g.end = int(fields[4])
            
        elif fields[2] == 'transcript' and info['gene_biotype'] == 'protein_coding':
            tid = info['transcript_id']
            if tid not in id2ent: id2ent[tid] = Transcript()
            t = id2ent[tid]
            t.chrm = fields[0]
            t.strand = fields[6]
            t.beg = int(fields[3])
            t.end = int(fields[4])
            t.name = info['transcript_id']
            gene_id = info['gene_id']
            if gene_id not in id2ent: id2ent[gene_id] = Gene()
            t.gene = id2ent[gene_id]
            t.gene.tpts.append(t)
            t.source = 'Ensembl'
            cnt += 1
        elif fields[2] == 'exon' and info['gene_biotype'] == 'protein_coding':
            tid = info['transcript_id']
            if tid not in id2ent: id2ent[tid] = Transcript()
            t = id2ent[tid]
            t.exons.append((int(fields[3]), int(fields[4])))
        elif fields[2] == 'CDS' and info['gene_biotype'] == 'protein_coding':
            tid = info['transcript_id']
            if tid not in id2ent: id2ent[tid] = Transcript()
            t = id2ent[tid]
            t.cds.append((int(fields[3]), int(fields[4])))

    sys.stderr.write("[%s] Loaded %d transcripts from Ensembl GTF file.\n" % (__name__, cnt))

def parse_ccds_table(ccds_fn, name2gene):

    """ start 0-based end 0-based """

    ccds_fh = open(ccds_fn)
    ccds_fh.readline()
    cnt = 0
    for line in ccds_fh:
        fields = line.strip().split('\t')
        if fields[5] != 'Public':
            continue
        gene_name = fields[2].upper()
        if gene_name not in name2gene:
            name2gene[gene_name] = Gene(gene_name)

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
        if cnt > 1000:
            break
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        info = dict(re.findall(r'\s*([^"]*) "([^"]*)";', fields[8]))
        if fields[2] == 'gene' and info['gene_type'] != 'pseudogene':
            gene_name = info['gene_name'].upper()
            gid = info['gene_id']
            if gene_name in name2gene:
                g = name2gene[gene_name]
                id2ent[gid] = g
            else:
                if gid not in id2ent: id2ent[gid] = Gene(gene_name)
                g = id2ent[gid]
                name2gene[gene_name] = g
            g.beg = int(fields[3])
            g.end = int(fields[4])
            # if info['gene_type'] == 'pseudogene':
            #     g.pseudo = True
        elif fields[2] == 'transcript' and info['gene_type'] == 'protein_coding':
            tid = info['transcript_id']
            if tid not in id2ent: id2ent[tid] = Transcript()
            t = id2ent[tid]
            t.chrm = fields[0]
            t.strand = fields[6]
            t.beg = int(fields[3])
            t.end = int(fields[4])
            t.name = tid
            gid = info['gene_id']
            if gid not in id2ent: id2ent[gid] = Gene()
            t.gene = id2ent[gid]
            t.gene.tpts.append(t)
            t.source = 'GENCODE'
            id2ent[t.name] = t
            cnt += 1
        elif fields[2] == 'exon' and info['gene_type'] == 'protein_coding':
            tid = info['transcript_id']
            if tid not in id2ent: id2ent[tid] = Transcript()
            t = id2ent[tid]
            t.exons.append((int(fields[3]), int(fields[4])))
        elif fields[2] == 'CDS' and info['gene_type'] == 'protein_coding':
            tid = info['transcript_id']
            if tid not in id2ent: id2ent[tid] = Transcript()
            t = id2ent[tid]
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
            gene_name = info['gene_id'].upper()
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
            gene_name = info['gene_id'].upper()
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

def parse_uniprot_mapping(fn):

    tid2uniprot = {}
    for line in opengz(fn):
        fields = line.strip().split('\t')
        tid2uniprot[fields[2]] = fields[0]

    err_print('[%s] Loaded %d transcript with UniProt mapping.' % (__name__, len(tid2uniprot)))

    return tid2uniprot

def extend_taa_seq(taa_pos_base, old_seq, new_seq, tpt):

    taa_pos = None
    termlen = None
    seq_end = tpt.cds_end
    i = 0
    while True:
        ci = i*3
        old_codon_seq = old_seq[ci:ci+3]
        new_codon_seq = new_seq[ci:ci+3]
        # if sequence comes to ends, extend sequence from reference file
        if (old_codon_seq not in standard_codon_table or 
            new_codon_seq not in standard_codon_table):
            seq_inc = faidx.refgenome.fetch_sequence(tpt.chrm, seq_end+1, seq_end+100)
            old_seq += seq_inc
            new_seq += seq_inc
            old_codon_seq = old_seq[ci:ci+3]
            new_codon_seq = new_seq[ci:ci+3]
            seq_end += 100

        taa_ref_run = codon2aa(old_codon_seq)
        taa_alt_run = codon2aa(new_codon_seq)
        #print i, old_codon_seq, new_codon_seq, taa_ref_run, taa_alt_run
        if taa_pos == None and taa_ref_run != taa_alt_run:
            taa_pos = i
            taa_ref = taa_ref_run
            taa_alt = taa_alt_run
        if taa_alt_run == '*':
            if taa_pos == None:
                # Terminating codon encountered before difference
                return None     # nothing occur to protein level
            termlen = i + 1 - taa_pos
            break
        i += 1

    if taa_pos == None:
        print 'oldseq', old_seq
        print 'newseq', new_seq
    taa_pos += taa_pos_base
    return taa_pos, taa_ref, taa_alt, str(termlen)

def translate_seq(seq):

    if len(seq) % 3 != 0:
        raise IncompatibleTranscriptError('translated coding sequence not multiplicative of 3, most likely a truncated sequence.')

    aa_seq = []
    for i in xrange(len(seq)/3):
        aa = codon2aa(seq[i*3:i*3+3])
        aa_seq.append(aa)
        if aa == '*':
            break

    return ''.join(aa_seq)

# def tgpos2codon(thash, chrm, pos):

#     for t in thash.get_transcripts(chrm, pos):
#         c, p, r = t.gpos2codon(chrm, pos)
#         yield t, c, p, r

# def tgpos2codon_longest(thash, chrm, pos):

#     l = None
#     for t, c, p, r in gpos2codon(thash,chrm,pos):
#         if (not l) or len(l) < len(t):
#             l = t
#             res = (t, c, p, r)

#     if l:
#         yield res

