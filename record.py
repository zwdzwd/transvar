import re
from utils import *
from err import *

class Pos():

    def __init__(self, pos='', tpos=0):

        self.pos = pos
        self.tpos = tpos         # respect to exon boundary, non-zero value indicates the position is relative to exon boundary

    def __repr__(self):
        if self.tpos < 0:
            return '%s%d' % (str(self.pos), self.tpos)
        elif self.tpos > 0:
            return '%s+%d' % (str(self.pos), self.tpos)
        else: return str(self.pos)

def parse_pos(posstr):

    if posstr.isdigit():
        p = Pos()
        p.pos = int(posstr)
        p.tpos = 0
    else:
        m = re.match(r'(\d+)([+-]\d+)', posstr)
        if not m: err_raise(InvalidInputError, 'Invalid position string %s.' % posstr, __name__)
        p = Pos()
        p.pos = int(m.group(1))
        p.tpos = int(m.group(2))

    return p

class Query(object):

    def __init__(self):

        self.op = None
        self.is_codon = True
        self.gn_name = None
        self.tpt = ''

    def set_pos(self, pos_str):

        if (pos_str.isdigit() and int(pos_str) > 0):
            self.pos = int(pos_str)
            return True
        else:
            err_warn('Abnormal position %s. skip.\n' % pos_str, __name__)
            return False

class QuerySNV(Query):

    def __init__(self):

        super(QuerySNV, self).__init__()
        self.pos = ''
        self.ref = ''
        self.alt = ''

    def cpos(self):
        return self.pos.pos

class QueryDEL(Query):

    def __init__(self):

        super(QueryDEL, self).__init__()
        self.beg = None
        self.end = None
        self.delseq = ''

class QueryINS(Query):

    def __init__(self):

        super(QueryINS, self).__init__()
        self.pos = ''           # position of base before
        self.insseq = ''

class QueryMNV(Query):

    def __init__(self):

        super(QueryMNV, self).__init__()
        self.beg = None
        self.end = None
        self.refseq = ''
        self.altseq = ''

template = "{r.chrm}\t{r.pos}\t{r.tname}\t{r.reg}\t{gnuc}/{tnuc}/{taa}\t{r.info}"
class Record():

    def __init__(self, muttype=''):

        self.tname = '.'        # transcript name
        self.chrm = '.'         # genomic chromosome
        self.pos = '.'          # genomic position string
        self.reg = '.'          # region
        self.info = '.'         # ;-separated key=value pair
        self.muttype = muttype

    def tnuc(self):
        """ format in HGVS nomenclature e.g., c.12345A>T """
        s = 'c.'
        if self.muttype:
            if hasattr(self, 'tnuc_range') and self.tnuc_range: s += self.tnuc_range
            if s == 'c.': return '.'
        else:
            if hasattr(self, 'tnuc_pos') and self.tnuc_pos: s += str(self.tnuc_pos)
            if hasattr(self, 'tnuc_ref') and self.tnuc_ref: s += self.tnuc_ref
            s += '>'
            if hasattr(self, 'tnuc_alt') and self.tnuc_alt: s += self.tnuc_alt
            if s == 'c.>': return '.'
        return s

    def gnuc(self):
        """ format in chr1:A12345T """
        s = self.chrm+':'
        if self.muttype:
            if hasattr(self, 'gnuc_range') and self.gnuc_range: s += self.gnuc_range
            if s == '.:': return '.'
        else:
            if hasattr(self, 'gnuc_ref') and self.gnuc_ref: s += self.gnuc_ref
            if hasattr(self, 'gnuc_pos') and self.gnuc_pos: s += str(self.gnuc_pos)
            if hasattr(self, 'gnuc_alt') and self.gnuc_alt: s += self.gnuc_alt
            if s == '.:': return '.'
        return s

    def taa(self):
        """ format in HGVS nomenclature e.g., p.E545K """
        s = 'p.'
        if self.muttype:
            if hasattr(self, 'taa_range') and self.taa_range: s += self.taa_range
            if s == 'p.': return '.'
        else:
            if hasattr(self, 'taa_ref') and self.taa_ref: s += self.taa_ref
            if hasattr(self, 'taa_pos') and self.taa_pos: s += str(self.taa_pos)
            if hasattr(self, 'taa_alt') and self.taa_alt: s += self.taa_alt
            if s == 'p.': return '.'
        return s

    def format(self, op):
        s = op+'\t' if op else ''
        s += template.format(r=self,
                             gnuc=self.gnuc(), tnuc = self.tnuc(), taa = self.taa())
        print s
