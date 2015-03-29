import sys, re
from utils import *
from io_utils import *
from err import *
import faidx
from record import *
from collections import deque

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

def tnuc_region_in_exon(np, beg, end):
    """ region in tnuc positions """

    if beg.tpos != 0: return False
    if end.tpos != 0: return False
    for i in xrange(beg.pos, end.pos-1):
        if abs(np[i] - np[i+1]) != 1:
            return False
    return True

def tnuc_region_in_intron(np, beg, end):
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

        if (not seq) or (len(seq) != self.end - self.beg + 1):
            raise SequenceRetrievalError()
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

    def taa_range2tnuc_seq(self, taa_beg, taa_end):

        if taa_beg*3 > len(self) or taa_end*3 > len(self):
            raise IncompatibleTranscriptError('codon nonexistent')

        self.ensure_seq()
        return self.seq[taa_beg*3-3:taa_end*3]

    def taa_range2aa_seq(self, taa_beg, taa_end):

        return translate_seq(self.taa_range2tnuc_seq(taa_beg, taa_end))

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

    def gnuc_seq2tnuc(self, gnuc_seq):
        if self.strand == '+':
            return gnuc_seq
        else:
            return reverse_complement(gnuc_seq)

    def cpos2aa(self, cpos):
        self.ensure_seq()
        return translate_seq(self.seq[cpos*3-3:cpos*3])

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

    def describe(self, gpos):
        """ determine the position of a single site """

        rg = RegAnno()
        if gpos < self.cds_beg:
            rg.UTR = '5' if self.strand == '+' else '3'
        if gpos > self.cds_end:
            rg.UTR = '3' if self.strand == '+' else '5'

        for i, exon in enumerate(self.exons):
            exind = i+1 if self.strand == '+' else len(self.exons) - i
            if i == 0 and exon[0] > gpos:
                rg.intergenic = 'Upstream'
                return rg
            if i == len(self.exons)-1 and exon[1] < gpos:
                rg.intergenic = 'Downstream'
                return rg
            if exon[0] <= gpos and exon[1] >= gpos: # exonic
                rg.exonic = True
                if gpos >= self.cds_beg and gpos <= self.cds_end:
                    rg.cds = True
                if gpos == self.cds_beg:
                    rg.start = True
                if gpos == self.cds_end:
                    rg.stop = True
                if gpos == exon[1]:
                    rg.splice = 'NextToDonor' if self.strand == '+' else 'NextToAcceptor'
                if gpos == exon[0]:
                    rg.splice = 'NextToAcceptor' if self.strand == '+' else 'NextToDonor'
                rg.exon = exind
                return rg
            if i > 0:
                pexon = self.exons[i-1]
                if gpos > pexon[1] and gpos < exon[0]: # intronic
                    rg.intronic = True
                    if self.strand == '+':
                        rg.intron_exon1 = exind-1
                        rg.intron_exon2 = exind
                    else:
                        rg.intron_exon1 = exind
                        rg.intron_exon2 = exind+1
                    if gpos in [pexon[1]+1, pexon[1]+2]:
                        rg.splice = 'Donor' if self.strand == '+' else 'Acceptor'
                    if gpos in [exon[0]-2, exon[0]-1]:
                        rg.splice = 'Acceptor' if self.strand == '-' else 'Donor'
                    return rg

        raise Exception()       # you shouldn't reach here

    def describe_span(self, gnuc_beg, gnuc_end):

        rg = RegSpanAnno()
        rg.b1 = self.describe(gnuc_beg)
        rg.b2 = self.describe(gnuc_end)
        rg.transcript_regs = self.overlap_region(gnuc_beg, gnuc_end)

        return rg

    def _gpos2codon_p(self, gpos, np):

        if gpos < self.cds_beg:
            p = Pos(1, gpos-self.cds_beg)
            c = self._init_codon_(1)
            c.seq = self.seq[:3]
            c.locs = np[:3]
            reg = self.describe(gpos)
            return c, p, reg

        if gpos > self.cds_end:
            p = Pos(len(self.seq), gpos-self.cds_end)
            c = self._init_codon_((len(self.seq)+2)/3)
            c.seq = self.seq[c.index*3-3:c.index*3]
            c.locs = np[c.index*3-3:c.index*3]
            reg = self.describe(gpos)
            return c, p, reg
        
        for i, pos in enumerate(np):
            if gpos == pos:
                c = self._init_codon_(i/3+1)
                c.seq    = self.seq[i-i%3:i-i%3+3]
                c.locs   = np[i-i%3:i-i%3+3]
                p = Pos(i+1, 0)
                reg = self.describe(gpos)
                return c, p, reg
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
                reg = self.describe(gpos)
                return c, p, reg

    def _gpos2codon_n(self, gpos, np):

        if gpos < self.cds_beg:
            p = Pos(len(self.seq), self.cds_beg-gpos)
            c = self._init_codon_(1)
            c.seq = self.seq[:3]
            c.locs = np[:3]
            reg = self.describe(gpos)
            return c, p, reg

        if gpos > self.cds_end:
            p = Pos(1, self.cds_end-gpos)
            c = self._init_codon_((len(self.seq)+2)/3)
            c.seq = self.seq[c.index*3-3:c.index*3]
            c.locs = np[c.index*3-3:c.index*3]
            reg = self.describe(gpos)
            return c, p, reg

        for i, pos in enumerate(np):
            if gpos == pos:
                c = self._init_codon_(i/3+1)
                c.seq = self.seq[i-i%3:i-i%3+3]
                c.locs = tuple(reversed(np[i-i%3:i-i%3+3]))
                p = Pos(i+1, 0)
                reg = self.describe(gpos)
                return c, p, reg
            
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
                reg = self.describe(gpos)
                return c, p, reg

    def gpos2codon(self, gpos):
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

        return regc

    def taa_roll_left_ins(self, index, taa_insseq):

        """ index is the position where the insertion comes after
        """

        self.ensure_seq()
        _taa_insseq_ = deque(taa_insseq)
        while True:
            if index <= 1:
                break
            rightmost = _taa_insseq_[-1]
            left_aa = translate_seq(
                self.seq[(index-1)*3:index*3])
            if rightmost != left_aa:
                break
            _taa_insseq_.pop()
            _taa_insseq_.appendleft(left_aa)
            index -= 1

        return index, ''.join(_taa_insseq_)

    def taa_roll_right_ins(self, index, taa_insseq):

        """ index is the position where the insertion comes after
        """

        self.ensure_seq()
        _taa_insseq_ = deque(taa_insseq)
        taa_len = len(self.seq) / 3
        while True:
            if index + 1 >= taa_len:
                break
            leftmost = _taa_insseq_[0]
            right_aa = translate_seq(
                self.seq[index*3:(index+1)*3])
            # print leftmost, right_aa, index
            if leftmost != right_aa:
                break
            _taa_insseq_.popleft()
            _taa_insseq_.append(right_aa)
            index += 1

        return index, ''.join(_taa_insseq_)

    def taa_roll_3p_ins(self, index, insseq):

        """ roll to 3' """
        if self.strand == '+':
            return self.taa_roll_left_ins(index, insseq)
        else:
            return self.taa_roll_right_ins(index, insseq)

    def taa_roll_left_del(self, taa_beg, taa_end):

        self.ensure_seq()
        while True:
            if taa_beg <= 1:
                break
            left_aa = self.cpos2aa(taa_beg-1)
            rightmost = self.cpos2aa(taa_end)
            if left_aa != rightmost:
                break
            taa_beg -= 1
            taa_end -= 1

        return taa_beg, taa_end

    def taa_roll_right_del(self, taa_beg, taa_end):

        self.ensure_seq()
        taa_len = len(self.seq) / 3
        while True:
            if taa_end + 1 >= taa_len:
                break
            right_aa = self.cpos2aa(taa_end+1)
            leftmost = self.cpos2aa(taa_beg)
            if leftmost != right_aa:
                break
            taa_beg += 1
            taa_end += 1

        return taa_beg, taa_end

    def tnuc_roll_left_ins(self, p, tnuc_insseq):

        """ p is the position where insertion comes after """

        self.ensure_seq()
        _tnuc_insseq_ = deque(tnuc_insseq)
        while True:
            if p <= 1:
                break
            left_base = self.seq[p-1]
            right_most = _tnuc_insseq_[-1]
            # print p, left_base, right_most
            if left_base != right_most:
                break
            _tnuc_insseq_.pop()
            _tnuc_insseq_.appendleft(left_base)
            p -= 1

        return p, ''.join(_tnuc_insseq_)

    def tnuc_roll_right_ins(self, p, tnuc_insseq):

        self.ensure_seq()
        _tnuc_insseq_ = deque(tnuc_insseq)
        tnuc_len = len(self.seq)
        while True:
            if p + 1 >= tnuc_len:
                break
            right_base = self.seq[p]
            left_most = _tnuc_insseq_[0]
            # print p, left_most, right_base
            if right_base != left_most:
                break
            _tnuc_insseq_.popleft()
            _tnuc_insseq_.append(right_base)
            p += 1

        return p, ''.join(_tnuc_insseq_)

    def tnuc_roll_left_del(self, beg, end):

        self.ensure_seq()
        while True:
            if beg <= 1:
                break
            left_base = self.seq[beg-2]
            right_most = self.seq[end-1]
            if left_base != right_most:
                break
            beg -= 1
            end -= 1

        return beg, end

    def tnuc_roll_right_del(self, beg, end):

        self.ensure_seq()
        tnuc_len = len(self.seq)
        while True:
            if end >= tnuc_len - 1:
                break
            right_base = self.seq[end]
            left_most = self.seq[beg-1]
            if right_base != left_most:
                break
            beg += 1
            end += 1

        return beg, end

    def getseq(self, beg, end):

        self.ensure_seq()
        return self.seq[beg-1:end]

    def tnuc_del_id(self, pbeg, pend):

        if pbeg == pend:
            tnuc_posstr = str(pbeg)
        else:
            tnuc_posstr = '%s_%s' % (pbeg, pend)

        tnuc_delseq = self.getseq(pbeg, pend)
        if len(tnuc_delseq) > delrep_len:
            tnuc_delrep = str(len(tnuc_delseq))
        else:
            tnuc_delrep = tnuc_delseq

        return '%sdel%s' % (tnuc_posstr, tnuc_delrep)

    def taa_del_id(self, taa_beg, taa_end):

        if taa_beg == taa_end:
            taa_posstr = '%s%d' % (
                self.cpos2aa(taa_beg), taa_beg)
        else:
            taa_posstr = '%s%d_%s%d' % (
                self.cpos2aa(taa_beg), taa_beg,
                self.cpos2aa(taa_end), taa_end)
        taa_del_len = taa_end - taa_beg + 1
        if taa_del_len > delrep_len:
            taa_delrep = str(taa_del_len)
        else:
            taa_delrep = self.taa_range2aa_seq(taa_beg, taa_end)

        return '%sdel%s' % (taa_posstr, taa_delrep)

    def taa_ins_id(self, index, taa_insseq):

        aa = self.cpos2aa(index)
        aa2 = self.cpos2aa(index+1)
        return '%s%d_%s%dins%s' % (aa, index, aa2, index+1, taa_insseq)

    def extend_taa_seq(self, taa_pos_base, old_seq, new_seq):

        taa_pos = None
        termlen = None
        seq_end = self.cds_end
        i = 0
        while True:
            ci = i*3
            old_codon_seq = old_seq[ci:ci+3]
            new_codon_seq = new_seq[ci:ci+3]
            # if sequence comes to ends, extend sequence from reference file
            if (old_codon_seq not in standard_codon_table or 
                new_codon_seq not in standard_codon_table):
                seq_inc = faidx.refgenome.fetch_sequence(self.chrm, seq_end+1, seq_end+100)
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


    def tnuc_mnv_coding_inframe(self, beg, end, altseq, r):

        """ beg and end are integer tnuc positions
        altseq follows the tnuc (cDNA) order
        set taa range
        """

        beg_codon_index = (beg + 2) / 3
        end_codon_index = (end + 2) / 3

        beg_codon_beg = beg_codon_index*3 - 2
        end_codon_end = end_codon_index*3 # 1 past the last codon

        old_seq = self.seq[beg_codon_beg-1:end_codon_end]
        new_seq = self.seq[beg_codon_beg-1:beg-1]+altseq+self.seq[end:end_codon_end]
        old_taa_seq = translate_seq(old_seq)
        new_taa_seq = translate_seq(new_seq)
        if old_taa_seq == new_taa_seq:
            return '(=)'

        # block substitution in nucleotide level may end up
        # an insertion or deletion on the protein level
        old_taa_seq1, new_taa_seq1, head_trim, tail_trim = double_trim(old_taa_seq, new_taa_seq)
        if not old_taa_seq1:
            _beg_index = beg_codon_index + head_trim - 1
            _end_index = beg_codon_index + head_trim
            _beg_aa = codon2aa(self.seq[_beg_index*3-3:_beg_index*3])
            _end_aa = codon2aa(self.seq[_end_index*3-3:_end_index*3])
            return '%s%d_%s%dins%s' % (
                _beg_aa, _beg_index,
                _end_aa, _end_index, new_taa_seq1)

        if not new_taa_seq1:
            if len(old_taa_seq1) == 1:
                return '%s%ddel' % (old_taa_seq1[0], beg_codon_index + head_trim)
            else:
                return '%s%d_%s%ddel' % (
                    old_taa_seq1[0], beg_codon_index + head_trim, 
                    old_taa_seq1[1], end_codon_index - tail_trim)
        if len(old_taa_seq1) == 1:
            if len(new_taa_seq1) == 1:
                return '%s%d%s' % (
                    old_taa_seq1[0], beg_codon_index + head_trim, new_taa_seq1)
            else:
                return '%s%ddelins%s' % (
                    old_taa_seq1[0], beg_codon_index + head_trim, new_taa_seq1)
        return '%s%d_%s%ddelins%s' % (
            old_taa_seq1[0], beg_codon_index + head_trim,
            old_taa_seq1[-1], end_codon_index - tail_trim, new_taa_seq1)

    def tnuc_mnv_coding_frameshift(self, beg, end, altseq, r):

        beg_codon_index = (beg + 2) / 3
        beg_codon_beg = beg_codon_index * 3 - 2
        old_seq = self.seq[beg_codon_beg-1:]
        new_seq = self.seq[beg_codon_beg-1:beg-1]+altseq+self.seq[end:]

        ret = self.extend_taa_seq(beg_codon_index, old_seq, new_seq)
        if ret:
            taa_pos, taa_ref, taa_alt, termlen = ret
            return '%s%d%sfs*%s' % (taa_ref, taa_pos, taa_alt, termlen)
        else:
            return '(=)'

    def tnuc_mnv_coding(self, beg, end, altseq, r):

        if (len(altseq) - (end-beg+1)) % 3 == 0:
            return self.tnuc_mnv_coding_inframe(beg, end, altseq, r)
        else:
            return self.tnuc_mnv_coding_frameshift(beg, end, altseq, r)

def gnuc_del_id(chrm, beg, end):

    if beg == end:
        gnuc_posstr = str(beg)
    else:
        gnuc_posstr = '%d_%d' % (beg, end)
    del_len = end - beg + 1
    if del_len > delrep_len:
        gnuc_delrep = str(del_len)
    else:
        gnuc_delrep = faidx.getseq(chrm, beg, end)

    return '%sdel%s' % (gnuc_posstr, gnuc_delrep)


def gnuc_roll_left_del(chrm, beg, end):

    """ beg and end are 1st and last base in the deleted sequence """

    sb = faidx.SeqBuf(chrm, beg)
    while True:
        if beg <= 1:
            break
        left_base = sb.get_base(chrm, beg-1)
        rightmost = sb.get_base(chrm, end)
        if left_base != rightmost:
            break
        beg -= 1
        end -= 1

    return beg, end

def gnuc_roll_right_del(chrm, beg, end):

    """ beg and end are 1st and last base in the deleted sequence """

    sb = faidx.SeqBuf(chrm, end)
    chrmlen = faidx.refgenome.chrm2len(chrm)
    while True:
        # check end of chromosome
        if end + 1 >= chrmlen:
            break

        right_base = sb.get_base(chrm, end+1)
        leftmost = sb.get_base(chrm, beg)
        if right_base != leftmost:
            break
        beg += 1
        end += 1

    return beg, end

def gnuc_roll_left_ins(chrm, pos, gnuc_insseq):

    """ pos is where insertion occur after """

    sb = faidx.SeqBuf(chrm, pos)
    _gnuc_insseq_ = deque(gnuc_insseq)
    while True:
        if pos <= 1:
            break
        left_base = sb.get_base(chrm, pos)
        rightmost = gnuc_insseq[-1]
        if left_base != rightmost:
            break
        _gnuc_insseq_.pop()
        _gnuc_insseq_.appendleft(left_base)
        pos -= 1

    return pos, ''.join(_gnuc_insseq_)

def gnuc_roll_right_ins(chrm, pos, gnuc_insseq):

    """ pos is where insertion occur after """

    sb = faidx.SeqBuf(chrm, pos)
    chrmlen = faidx.refgenome.chrm2len(chrm)
    _gnuc_insseq_ = deque(gnuc_insseq)
    while True:
        if pos + 1 >= chrmlen:
            break
        right_base = sb.get_base(chrm, pos+1)
        leftmost = gnuc_insseq[0]
        if right_base != leftmost:
            break
        _gnuc_insseq_.popleft()
        _gnuc_insseq_.append(right_base)
        pos += 1

    return pos, ''.join(_gnuc_insseq_)

class Gene():

    def __init__(self, name=''):

        self.name    = name
        self.dbxref  = ''       # for storing GENEID etc.
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

    def beg(self):
        if hasattr(self, 'beg'):
            return self.beg
        else:
            return gene.longest_tpt().beg

    def end(self):
        if hasattr(self, 'end'):
            return self.end
        else:
            return gene.longest_tpt().end

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
        t.chrm = normalize_chrm(fields[2])
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
        t.chrm = normalize_chrm(fields[1])
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
                if hasattr(g, '_gene_id') and g._gene_id != info['ID']:
                    continue   # if a gene_name appears twice, then all the subsequent occurrences are all ignored.
            else:
                g = Gene(gene_name)
                name2gene[gene_name] = g
            g._gene_id = info['ID']
            g.beg = int(fields[3])
            g.end = int(fields[4])
            id2ent[info['ID']] = g
            if 'Dbxref' in info:
                g.dbxref = info['Dbxref']
        elif fields[2] == 'mRNA' and 'Parent' in info and info['Parent'] in id2ent:
            t = Transcript()
            t.chrm = normalize_chrm(reg.name)
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
            t.chrm = normalize_chrm(fields[0])
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
        t.chrm = normalize_chrm(fields[0])
        t.strand = fields[6]
        t.cds_beg = int(fields[7])+1
        t.cds_end = int(fields[8])+1

        # without UTR information, take CDS boundary as the exon boundary
        t.beg = t.cds_beg
        t.end = t.cds_end

        t.name = fields[4]
        # note that CCDS do not annotate UTR, so all the exons are equivalently cds
        t.exons = [(int(b)+1, int(e)+1) for b,e in re.findall(r"[\s\[](\d+)-(\d+)[,\]]", fields[9])]
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
        t.chrm = normalize_chrm(fields[1])
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
        # if cnt > 1000:
        #     break
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
            t.chrm = normalize_chrm(fields[0])
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
        if len(fields) < 9: continue # the old transcript definition (hg18) from AceView is a bit corrupted.
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
            if gene_name in name2gene:
                g = name2gene[gene_name]
            else:
                g = Gene(gene_name)
                name2gene[gene_name] = g

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

