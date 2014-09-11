class Query():

    def __init__(self):

        self.op = None
        self.is_codon = True
        self.gn_name = None
        self.pos = ''
        self.ref = ''
        self.alt = ''
        self.tpt = ''

    def set_pos(self, pos_str):

        if (pos_str.isdigit() and int(pos_str) > 0):
            self.pos = int(pos_str)
            return True
        else:
            sys.stderr.write('[Warning] abnormal position %s. skip.\n' % pos_str)
            return False

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
        if hasattr(self, 'tnuc_pos') and self.tnuc_pos: s += str(self.tnuc_pos)
        if hasattr(self, 'tnuc_ref') and self.tnuc_ref: s += self.tnuc_ref
        s += '>'
        if hasattr(self, 'tnuc_alt') and self.tnuc_alt: s += self.tnuc_alt
        if s == 'c.>': return '.'
        return s

    def gnuc(self):
        """ format in chr1:A12345T """
        s = self.chrm+':'
        if hasattr(self, 'gnuc_ref') and self.gnuc_ref: s += self.gnuc_ref
        if hasattr(self, 'gnuc_pos') and self.gnuc_pos: s += str(self.gnuc_pos)
        if hasattr(self, 'gnuc_alt') and self.gnuc_alt: s += self.gnuc_alt
        if s == '.:': return '.'
        return s

    def taa(self):
        """ format in HGVS nomenclature e.g., p.E545K """
        s = 'p.'
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
