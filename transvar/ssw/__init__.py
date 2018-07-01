import ctypes, os
_DIRNAME=os.path.abspath(os.path.dirname(__file__))
so_file = [f for f in os.listdir(_DIRNAME) if f.startswith('_sswlib') and f.endswith('.so')][0]
ssw=ctypes.CDLL(os.path.join(_DIRNAME, so_file))

class SSWAlign(ctypes.Structure):
    _fields_ = [('score', ctypes.c_uint32),
                ('qbeg', ctypes.c_int),
                ('qend', ctypes.c_int),
                ('tbeg', ctypes.c_int),
                ('tend', ctypes.c_int),
                ('ctype', ctypes.POINTER(ctypes.c_int)),
                ('clen', ctypes.POINTER(ctypes.c_int)),
                ('ncigar', ctypes.c_int)]

ssw.aln_gap.restype = ctypes.POINTER(SSWAlign)
ssw.aln.restype = ctypes.POINTER(SSWAlign)

# call ssw.aln:
# aln = ssw.aln(qseq, tseq)
# print aln.contents.score
# print aln.contents.ncigar
# print aln.contents.tbeg
# print aln.contents.tend
# ssw.free_aln(aln)

class SSWAln:
    def __init__(self):
        
        self.score = -1
        self.qbeg = -1
        self.qend = -1
        self.rbeg = -1
        self.rend = -1
        self.cigar = []
        
    def __repr__(self):
        return "<Align sc:%d q:[%d-%d] r:[%d-%d] cigar:%s>" % (self.score, self.qbeg, self.qend, self.rbeg, self.rend, str(self.cigar))

from builtins import bytes    
def ssw_aln(qseq, rseq, gap=False):
    if gap:
        _aln = ssw.aln_gap(qseq.encode('ascii'), rseq.encode('ascii'))
    else:
        _aln = ssw.aln(qseq.encode('ascii'), rseq.encode('ascii'))

    aln = SSWAln()
    aln.score = _aln.contents.score
    aln.qbeg  = _aln.contents.qbeg
    aln.qend  = _aln.contents.qend
    aln.rbeg  = _aln.contents.tbeg
    aln.rend  = _aln.contents.tend
    if aln.qbeg > 0:
        aln.cigar.append((4, aln.qbeg))
    qlen = aln.qbeg
    for i in range(_aln.contents.ncigar):
        ct = _aln.contents.ctype[i]
        cl = _aln.contents.clen[i]
        aln.cigar.append((ct, cl))
        if ct == 0 or ct == 1:
            qlen += cl
    if qlen < len(qseq):
        aln.cigar.append((4, len(qseq)-qlen))
        
    ssw.free_aln(_aln)
    return aln
