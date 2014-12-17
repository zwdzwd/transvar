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

    def __eq__(self, other):
        if (self.pos == other.pos and
            self.tpos == other.tpos):
            return True
        else:
            return False

    def included_plus(self):
        if self.tpos > 0:
            return self.pos + 1
        else:
            return self.pos

    def included_minus(self):
        if self.tpos < 0:
            return self.pos - 1
        else:
            return self.pos


class RegAnno():

    """ annotation of a single site """

    def __init__(self):
        self.intronic = False
        self.exonic = False 
        self.cds = False
        self.UTR = None         # '3' or '5'
        self.exon = None
        self.intron_exon1 = None
        self.intron_exon2 = None
        self.intergenic = None  # 'Upstream' or 'Downstream'
        self.splice = None      # 'NextToDonor' | 'Donor' | 'Acceptor' | 'NextToAcceptor'

    def append(self, f, a):
        if f:
            return f+';'+a
        else:
            return a

    def format(self):
        f = ''
        if self.UTR: f = self.append(f, '%s-UTR' % self.UTR)
        if self.intronic:
            f = self.append(f, 'Intronic_%d_%d' %
                            (self.intron_exon1, self.intron_exon2))
        elif self.exonic:
            if self.cds:
                f = self.append(f, 'Exonic_%d' % self.exon)
            else:
                f = self.append(f, 'CDS_%d' % self.exon)
        elif self.intergenic:
            f = self.append(f, 'Intergenic%s' % self.intergenic)

        return f

class RegSpanAnno():

    """ annotation of a span """

    def __init__(self):

        self.whole_gene = False
        self.bp1 = None         # an object of RegAnno
        self.bp2 = None         # an object of RegAnno
        self.span_features = []

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

        """ for a region by default, no mutation information included """
        self.beg = ''
        self.end = ''
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

class QueryREG(Query):

    def __init__(self):

        super(QueryREG, self).__init__()
        self.beg = ''
        self.end = ''
        self.refseq = ''

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
        self.beg = ''
        self.end = ''
        self.delseq = ''
        # for amino acid
        self.beg_aa = ''
        self.end_aa = ''

class QueryFrameShift(Query):

    def __init__(self):

        super(QueryFrameShift, self).__init__()
        self.pos = None
        self.ref = ''
        self.alt = ''
        self.stop_index = ''

class QueryINS(Query):

    def __init__(self):

        super(QueryINS, self).__init__()
        self.pos = ''           # position of base before
        self.insseq = ''
        # for amino acid level query
        self.beg = ''
        self.beg_aa = ''
        self.end = ''
        self.end_aa = ''

class QueryMNV(Query):

    def __init__(self):

        super(QueryMNV, self).__init__()
        self.beg = ''
        self.beg_aa = ''
        self.end = ''
        self.end_aa = ''
        self.refseq = ''
        self.altseq = ''

class QueryDUP(Query):

    def __init__(self):

        super(QueryDUP, self).__init__()
        self.beg = ''
        self.beg_aa = ''
        self.end = ''
        self.end_aa = ''
        self.dupseq = ''


def normalize_reg(q):

    """ create a sensable region 
    respect to the length of the chromosome """

    if q.beg > reflen(q.tok):
        err_print('Region beg %d greater than chromosome length %d, truncated.'
                  % (q.beg, reflen(q.tok)))
        q.beg = reflen(q.tok)

    if q.end > reflen(q.tok):
        err_print('Region end %d greater than chromosome length %d, truncated.'
                  % (q.end, reflen(q.tok)))
        q.end = reflen(q.tok)
    if q.beg < 0:
        err_print('Region beg %d negative, truncated to 0.')
        q.beg = 0    
    if q.end < 0:
        err_print('Region end %d negative, truncated to 0.')
        q.end = 0    

template = "{r.chrm}\t{r.pos}\t{r.tname}\t{r.reg}\t{gnuc}/{tnuc}/{taa}\t{r.info}"
class Record():

    def __init__(self):

        self.tname = '.'        # transcript name
        self.chrm = '.'         # genomic chromosome
        self.pos = '.'          # genomic position string
        self.reg = '.'          # region
        self.info = '.'         # ;-separated key=value pair

    def tnuc(self):
        """ format in HGVS nomenclature e.g., c.12345A>T """
        s = 'c.'
        if hasattr(self, 'tnuc_range') and self.tnuc_range:
            s += self.tnuc_range
            if s == 'c.': return '.'
        else:
            if hasattr(self, 'tnuc_pos') and self.tnuc_pos: s += str(self.tnuc_pos)
            if hasattr(self, 'tnuc_ref') and self.tnuc_ref: s += self.tnuc_ref
            s += '>'
            if hasattr(self, 'tnuc_alt') and self.tnuc_alt: s += self.tnuc_alt
            if s == 'c.>': return '.'
        return s

    def append_info(self, app):
        if self.info and self.info != '.':
            self.info += ';'+app
        else:
            self.info = app

    def gnuc(self):
        
        """ format in chr1:A12345T """
        s = self.chrm+':g.'
        if hasattr(self, 'gnuc_range') and self.gnuc_range:
            s += self.gnuc_range
        else:
            if hasattr(self, 'gnuc_pos') and self.gnuc_pos: s += str(self.gnuc_pos)
            if hasattr(self, 'gnuc_ref') and self.gnuc_ref: s += self.gnuc_ref
            s += '>'
            if hasattr(self, 'gnuc_alt') and self.gnuc_alt: s += self.gnuc_alt
        if s == '.:g.>': return '.'
        return s

    def taa(self):
        """ format in HGVS nomenclature e.g., p.E545K """
        s = 'p.'
        if hasattr(self, 'taa_range') and self.taa_range:
            s += self.taa_range
        else:
            if hasattr(self, 'taa_ref') and self.taa_ref: s += self.taa_ref
            if hasattr(self, 'taa_pos') and self.taa_pos: s += str(self.taa_pos)
            if hasattr(self, 'taa_alt') and self.taa_alt: s += self.taa_alt
        if s == 'p.': return '.'
        return s

    def format_id(self):
        return '%s/%s/%s' % (self.gnuc(), self.tnuc(), self.taa())

    def format(self, op):
        s = op+'\t' if op else ''
        s += template.format(r=self,
                             gnuc=self.gnuc(), tnuc = self.tnuc(), taa = self.taa())
        print s

