"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Wanding Zhou, Tenghui Chen, Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""
from __future__ import division
import sys
from .err import *
from . import tabix

MAXCHRMLEN=300000000
def normalize_chrm(chrm):

    if chrm == '23' or chrm == 'chr23': chrm = 'X'
    if chrm == '24' or chrm == 'chr24': chrm = 'Y'
    if chrm == '25' or chrm == 'chr25': chrm = 'M'
    if chrm == 'MT' or chrm == 'chrMT': chrm = 'M'
    if chrm.isdigit() or chrm in ['X', 'Y', 'M']:
        # not chrm.startswith('chr')
        chrm = 'chr'+chrm

    return chrm

def normalize_chrm_dbsnp(chrm):

    if chrm == '23' or chrm == 'chr23': chrm = 'X'
    if chrm == '24' or chrm == 'chr24': chrm = 'Y'
    if chrm == '25' or chrm == 'chr25': chrm = 'MT'
    if chrm == 'M' or chrm == 'chrM': chrm = 'MT'
    if chrm.startswith('chr'):
        chrm = chrm[3:]

    return chrm

aa_3to1_table = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Sec": "U",
    "Pyl": "O",
    "X": "*",
    # ambiguous aa
    "Glx": "Z",                 # Q or E
    "Asx": "B",                 # N or D
}

def aa_is_stop(aa):

    if aa == '*' or aa == 'X':
        return True
    else:
        return False

def aa_has_stop(aa):

    if '*' in aa or 'X' in aa:
        return True
    else:
        return False

aa_1to3_table = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "E": "Glu",
    "Q": "Gln",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "U": "Sec",
    "O": "Pyl",
    "*": "X",
    # ambiguous aa
    "Z": "Glx",                 # Q or E
    "B": "Asx",                 # N or D
}

def aa_3to1(aaseq):

    aaseq1 = ''
    for i in range(len(aaseq)//3):
        aa = aaseq[3*i:3*i+3]
        aaseq1 += aa_3to1_table[aa]

    return aaseq1

def aa_1to3(aaseq1, use_list=False):

    if use_list:
        aaseq3 = []
        for aa in aaseq1:
            aaseq3.append(aa_1to3_table[aa])
    else:
        aaseq3 = ''
        for aa in aaseq1:
            aaseq3 += aa_1to3_table[aa]

    return aaseq3

def aaf(aaseq, args, use_list=False):
    if args.aa3:
        return aa_1to3(aaseq, use_list=use_list)
    else:
        return aaseq

def printseq(seq, args):

    if len(seq) > args.seqmax:
        return '%s..%s' % (seq[:3], seq[-3:])
    else:
        return seq

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
        chrm = normalize_chrm(t.chrm)
        for ki in range(t.cds_beg//self.binsize, t.cds_end//self.binsize+1):
            k = (chrm, ki)
            self.add_transcript_by_key(k, t)
        # k1 = (chrm, t.cds_beg // self.binsize)
        # k2 = (chrm, t.cds_end // self.binsize)
        # if k1 == k2:
        #     self.add_transcript_by_key(k1, t)
        # else:
        #     self.add_transcript_by_key(k1, t)
        #     self.add_transcript_by_key(k2, t)

    def get_transcripts_cds(self, chrm, beg, end=None, flanking=0):

        """ get transcript if between CDS beginning and end """
        if not end: end = beg
        chrm = normalize_chrm(chrm)
        kbeg = int(beg) // self.binsize
        kend = int(end) // self.binsize
        ts = []
        for ki in range(kbeg, kend+1):
            k = (chrm, ki)
            if k in self.key2transcripts:
                for t in self.key2transcripts[k]:
                    if t.cds_beg-flanking <= end and t.cds_end+flanking >= beg:
                        ts.append(t)
        return ts

    def get_transcripts(self, chrm, beg, end=None, flanking=0):

        """ get transcript if between beginning and end """

        if not end: end = beg
        chrm = normalize_chrm(chrm)
        kbeg = int(beg) // self.binsize
        kend = int(end) // self.binsize
        ts = []
        for ki in range(kbeg, kend+1):
            k = (chrm, ki)
            # print ki, kbeg, kend, k
            if k in self.key2transcripts:
                for t in self.key2transcripts[k]:
                    if t.beg-flanking <= end and t.end+flanking >= beg:
                        ts.append(t)

        return ts

    def get_closest_transcripts_upstream(self, chrm, pos):
        pos = int(pos)
        chrm = normalize_chrm(chrm)
        for ki in range(pos//self.binsize, -1, -1):
            k = (chrm, ki)
            if k in self.key2transcripts:
                tpts = [t for t in self.key2transcripts[k] if t.end < pos]
                if tpts:
                    tpts.sort(key=lambda t: t.end, reverse=True)
                    return tpts[0]
        return None

    def get_closest_transcripts_downstream(self, chrm, pos):
        pos = int(pos)
        chrm = normalize_chrm(chrm)
        for ki in range(pos//self.binsize, MAXCHRMLEN//self.binsize, 1):
            k = (chrm, ki)
            if k in self.key2transcripts:
                tpts = [t for t in self.key2transcripts[k] if t.beg > pos]
                if tpts:
                    tpts.sort(key=lambda t: t.beg)
                    return tpts[0]
        return None

    # def get_closest_transcripts(self, chrm, beg, end):

    #     """ closest transcripts upstream and downstream """
    #     return (self.get_closest_transcripts_upstream(chrm, beg),
    #             self.get_closest_transcripts_downstream(chrm, end))

    # def __init__(self):
    #     # chrm => bin => list of transcripts
    #     self.chrm2bin = {}
    #     self.bins = 3000
    #     self.binsize = 100000

    # def insert(self, tpt):

    #     if tpt.chrm not in self.chrm2bin:
    #         self.chrm2bin[tpt.chrm] = [[]]*self.bins
    #     else:
    #         bin1 = self.chrm2bin[tpt.chrm][tpt.cds_beg//self.binsize]
    #         bin2 = self.chrm2bin[tpt.chrm][tpt.cds_end//self.binsize]
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
    #     for tpt in self.chrm2bin[chrm][pos // self.binsize]:
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

def get_config(config, option, rv=None):

    if not config.has_section(rv):
        err_warn('%s (%s) has no default value, please specify' % (option, rv))
        return None

    if config.has_option(rv, option):
        return config.get(rv, option)
    else:
        err_warn('%s (%s) has no default value, please specify' % (option, rv))
        return None

def replace_defaults(args, config):
    if args.refversion:
        rv = args.refversion
    elif 'refversion' in config.defaults():
        rv = config.get('DEFAULT', 'refversion')
    else:
        rv = 'hg19'

    args.refversion = rv

    def _set_arg_(argname, rv):
        # replace __DEF__ by the actual default from config file
        if getattr(args, argname) == '_DEF_':
            setattr(args, argname, get_config(config, argname, rv))

    argnames = ['ensembl', 'reference', 'refseq', 'ccds',
                'gencode', 'ucsc', 'kg', 'aceview']
    for argname in argnames:
        _set_arg_(argname, rv)

    # _set_arg_('uniprot', 'idmap')

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
            result.extend([lst[_] for _ in range(start, end)])

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
        fh = gzip.open(fn, 'rt')
    else:
        fh = open(fn, 'rt')

    return fh

def double_trim(seq1, seq2):

    # trim head
    head_trim = 0
    while seq1 and seq2:
        if seq1[0] == seq2[0]:
            head_trim += 1
            seq1 = seq1[1:]
            seq2 = seq2[1:]
        else:
            break

    # trim tail
    tail_trim = 0
    while seq1 and seq2:
        if seq1[-1] == seq2[-1]:
            tail_trim += 1
            seq1 = seq1[:-1]
            seq2 = seq2[:-1]
        else:
            break

    return (seq1, seq2, head_trim, tail_trim)

def tabix_query(index, chrm, beg, end):

    beg = max(0, beg)
    # the following is necessary when the end is magnitudes greater
    # than the chromosome length
    # end = min(end, faidx.refgenome.chrm2len(normalize_chrm(chrm)))
    try:
        return index.query(chrm, beg, end)
    except tabix.TabixError: # when the chromosome is unavailable
        return []
