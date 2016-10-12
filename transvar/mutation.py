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

import re, sys, argparse
from .utils import *
from .record import *
from .err import *

def _parse_gdna_mutation(s):

    # m = re.match(r'(g\.)?(\d+)(_(\d+))?(\.)?(del([atgcATGC\d]*))?(ins([atgcATGC]*))?(([atgcATGC?]*)>([atgcATGC?]*))?(dup([atgcATGC\d]*))?$', s)
    m = re.match(r'(g\.)?(\d+)(_(\d+))?(\.)?(del([atgcnATGCN\d]*))?(ins([atgcnATGCN]*))?(([atgcnATGCN?]*)>([atgcnATGCN?]*))?(dup([atgcnATGCN\d]*))?$', s)

    if not m:
        raise InvalidInputError('invalid_mutation_string_%s' % s)

    (_, _beg, _end_s, _end, _, _is_del, _d,
     _is_ins, _i, _is_sub, _ref, _alt, _is_dup, _dupseq) = m.groups()

    if _is_sub and len(_ref) <= 1 and len(_alt) <= 1:
        q = QuerySNV()
        q.pos = int(_beg)
        q.ref = _ref if _ref and _ref != '?' else ''
        q.alt = _alt if _alt and _alt != '?' else ''
    elif (_is_del and _is_ins and
          (_d == '1' or len(_d) == 1) and len(_i) == 1):
        q = QuerySNV()
        q.pos = int(_beg)
        q.ref = '' if _d.isdigit() else _d.upper()
        q.alt = _i.upper() if _i else ''
    elif _is_del and not _is_ins:
        q = QueryDEL()
        q.beg = int(_beg)
        q.end = int(_end) if _end else q.beg
        q.delseq = '' if _d.isdigit() else _d.upper()
    elif _is_ins and not _is_del:
        q = QueryINS()
        q.pos = int(_beg)
        if _i:
            q.insseq = _i.upper()
        else:
            raise InvalidInputError('insertion_without_inserted_sequence_%s' % s)
    elif _is_ins and _is_del:
        q = QueryMNV()
        q.beg = int(_beg)
        q.end = int(_end) if _end else q.beg
        if _d and not _d.isdigit(): q.refseq = _d.upper()
        q.altseq = _i.upper() if _i else ''
    elif _is_sub:
        q = QueryMNV()
        q.beg = int(_beg)
        q.end = int(_end) if _end else q.beg
        if _ref and len(_ref) != q.end-q.beg+1:
            q.end = q.beg + len(_ref) - 1
        q.refseq = _ref.upper() if _ref else ''
        q.altseq = _alt.upper() if _alt else ''
    elif _is_dup:
        q = QueryDUP()
        q.beg = int(_beg)
        q.end = int(_end) if _end else q.beg
        q.dupseq = _dupseq.upper() if _dupseq else ''
    else:
        # use only the beg and end
        # treat input as a region
        q = QueryREG()
        q.beg = int(_beg)
        q.end = int(_end) if _end else q.beg
        q.refseq = _ref.upper() if _ref else ''

    return q

def _parse_cdna_mutation(s):

    # m = re.match(r'(c\.)?([\d+-]+)(_([\d+-]+))?(\.)?(del([atgcnATGCN\d]*))?(ins([atgcnATGCN]*))?(([atgcnATGCN?]*)>([atgcnATGCN?]*))?(dup([atgcnATGCN\d]*))?$', s)
    m = re.match(r'(c\.)?([\d*+-]+)(_([\d*+-]+))?(\.)?(del([atgcnATGCN\d]*))?(ins([atgcnATGCN]*))?(([atgcnATGCN?]*)>([atgcnATGCN?]*))?(dup([atgcnATGCN\d]*))?$', s)
    if not m:
        raise InvalidInputError('invalid_mutation_string_%s' % s)

    (_, _beg, _end_s, _end, _, _is_del, _d,
     _is_ins, _i, _is_sub, _ref, _alt, _is_dup, _dupseq) = m.groups()

    if _is_sub and len(_ref) <= 1 and len(_alt) <= 1:
        q = QuerySNV()
        q.pos = parse_pos(_beg)
        q.ref = _ref if _ref and _ref != '?' else ''
        q.alt = _alt if _alt and _alt != '?' else ''
    elif (_is_del and _is_ins and
          (_d == '1' or len(_d) == 1) and len(_i) == 1):
        q = QuerySNV()
        q.pos = parse_pos(_beg)
        q.ref = '' if _d.isdigit() else _d.upper()
        q.alt = _i.upper() if _i else ''
    elif _is_del and not _is_ins:
        q = QueryDEL()
        q.beg = parse_pos(_beg)
        q.end = parse_pos(_end) if _end else q.beg
        q.delseq = '' if _d.isdigit() else _d.upper()
    elif _is_ins and not _is_del:
        q = QueryINS()
        q.pos = parse_pos(_beg)
        if _i: q.insseq = _i.upper()
        else: err_die('insertion without inserted sequence: %s.' % s)
    elif _is_ins and _is_del:
        q = QueryMNV()
        q.beg = parse_pos(_beg)
        q.end = parse_pos(_end) if _end else q.beg
        if _d and not _d.isdigit(): q.refseq = _d.upper()
        q.altseq = _i.upper() if _i else ''
    elif _is_sub:
        q = QueryMNV()
        q.beg = parse_pos(_beg)
        q.end = parse_pos(_end) if _end else q.beg
        q.refseq = _ref.upper() if _ref else ''
        q.altseq = _alt.upper() if _alt else ''
    elif _is_dup:
        q = QueryDUP()
        q.beg = parse_pos(_beg)
        q.end = parse_pos(_end) if _end else q.beg
        q.dupseq = _dupseq.upper() if _dupseq else ''
    else:
        # use only the beg and end
        # treat input as a region
        q = QueryREG()
        q.beg = parse_pos(_beg)
        q.end = parse_pos(_end) if _end else q.beg
        q.refseq = _ref.upper() if _ref else ''

    return q

def read_aa(aaseq):

    if (not aaseq) or aaseq.isdigit():
        return ''
    if len(aaseq) < 3:
        return aaseq.upper()
    if aaseq[0].isupper() and aaseq[1].islower() and aaseq[2].islower() and len(aaseq) % 3 == 0:
        try:
            aaseq1 = aa_3to1(aaseq)
        except KeyError:
            return aaseq.upper()
        return aaseq1
    else:
        return aaseq.upper()

def _parse_protein_mutation(s):

    success = False
    # test fs without alternative allele
    m = re.match(r'(p\.)?([A-Za-z*?]*)(\d+)([A-Za-z*?]+)?fs((Ter|[\*Xx])(\d+))?$', s)
    if m:
        _, _beg_aa, _beg_i, _alt, _hasterlen, _tersymbol, _stop_i, = m.groups()
        _is_fs = True
        _end_s, _end_aa, _end_i = (None, None, None)
        _is_del, _d, _is_ins, _i = (None, None, None, None)
        _has_ref, _ref = (None, None)
        success = True

    if not success:
        m = re.match(r'(p\.)?([A-Za-z*?]*)(\d+)(_([A-Za-z*?]*)(\d+))?(del([^i][A-Za-z*?\d]*)?)?(ins([A-Za-z*?]+))?>?([A-Za-z*?]+)?(fs((Ter|[\*Xx])(\d+))?)?(ref([A-Za-zx*]*))?$', s)

        if m:
            (_, _beg_aa, _beg_i, _end_s, _end_aa, _end_i, 
             _is_del, _d, _is_ins, _i, _alt, _is_fs, _hasterlen,
             _tersymbol, _stop_i, _has_ref, _ref) = m.groups()
            success = True

    if not success:
        raise InvalidInputError('invalid_mutation_string_%s' % s)

    _beg_aa = read_aa(_beg_aa)
    _end_aa = read_aa(_end_aa)
    _alt = read_aa(_alt)
    _ref = read_aa(_ref)
    _i = read_aa(_i)
    _d = read_aa(_d)

    if _is_fs:
        #print 'fs'
        q = QueryFrameShift()
        q.pos = int(_beg_i)
        q.ref = _beg_aa
        q.alt = _alt
        q.stop_index = int(_stop_i) if _hasterlen else -1
    elif _is_del and not _is_ins:
        #print 'del'
        q = QueryDEL()
        q.beg = int(_beg_i)
        q.end = int(_end_i) if _end_i else q.beg
        q.beg_aa = _beg_aa
        q.end_aa = _end_aa
        q.delseq = _d
    elif (_is_del and _is_ins and (_d == '1' or (_d and len(_d) == 1) or not _end_s) and len(_i) == 1):
        #print 'snv'
        q = QuerySNV()
        q.pos = int(_beg_i)
        q.ref = _beg_aa
        q.alt = _alt

    elif _is_del and _is_ins:
        #print 'mnv'
        q = QueryMNV()
        q.beg = int(_beg_i)
        q.end = int(_end_i) if _end_i else q.beg
        q.beg_aa = _beg_aa
        q.end_aa = _end_aa
        q.refseq = _d
        q.altseq = _i
    elif _is_ins and not _is_del:
        # print 'ins'
        q = QueryINS()
        q.beg = int(_beg_i)
        q.beg_aa = _beg_aa
        q.end = int(_end_i)
        q.end_aa = _end_aa
        q.insseq = _i
    elif not _end_s and (not _beg_aa or len(_beg_aa) == 1) and _alt:
        # print 'snv'
        q = QuerySNV()
        q.pos = int(_beg_i)
        q.ref = _beg_aa
        q.alt = _alt
    elif _i:
        # print 'mnv'
        q = QueryMNV()
        q.beg = int(_beg_i)
        q.end = int(_end_i) if _end_i else q.beg
        q.beg_aa = _beg_aa
        q.end_aa = _end_aa
        q.refseq = _d
        q.altseq = _i
    else:
        # print 'reg'
        q = QueryREG()
        q.beg = int(_beg_i)
        q.end = int(_end_i) if _end_i else q.beg
        q.beg_aa = _beg_aa
        q.end_aa = _end_aa
        q.refseq = _ref
        if (not q.refseq) and q.beg == q.end:
            q.refseq = q.beg_aa

    return q


def parse_mutation_str(mut_str, mut):

    mut_str = mut_str.strip()
    if mut == 'g':
        return _parse_gdna_mutation(mut_str)
    elif mut == 'c':
        return _parse_cdna_mutation(mut_str)
    elif mut == 'p':
        return _parse_protein_mutation(mut_str)
    else:
        raise InvalidInputError('invalid_mutation_string_%s' % mut_str)

def parse_tok_mutation_str(s, muttype):

    s = s.strip()
    if s.find(':') >= 0:
        tok, mut_str = s.split(':', 1)
        q = parse_mutation_str(mut_str, muttype)
        q.tok = tok
    else:
        q = QueryGENE()
        q.tok = s
        
    # if s.find(":") >= 0:
    #     tok, mut_str = s.split(':', 1)
    # elif s.find(".") >= 0:
    #     tok, mut_str = s.split('.', 1)
    # else:
    #     pairs = s.split(None, 1)
    #     if len(pairs) == 2:
    #         tok, mut_str = pairs
    #     elif len(pairs) == 1:
    #         q = QueryGENE()
    #         q.tok = s
    #         return q
    #     else:
    #         err_warn('invalid mutation string: %s' % s)
    #         raise InvalidInputError
    # q = parse_mutation_str(mut_str, muttype)
    # q.tok = tok

    return q

def _list_parse_mutation(args, fields, indices, muttype):

    # protein
    if args.g > 0 and args.p > 0: # gene, position, ref, alt in separate columns for SNV

        q = QuerySNV()
        q.op = '|'.join(indices.extract(fields))
        q.tok = fields[args.g-1].strip()
        q.pos = parse_pos(fields[args.p-1].strip())
        if args.r > 0: q.ref = fields[args.r-1].strip()
        if args.a > 0: q.alt = fields[args.a-1].strip()
        if args.t > 0: q.tpt = fields[args.t-1].strip()
        q.is_codon = True

    elif args.g > 0 and args.n > 0:
        if muttype == 'g':      # gDNA
            q = QuerySNV()
            q.op = '|'.join(indices.extract(fields))
            q.tok = fields[args.g-1].strip()
            q.pos = int(fields[args.n-1].strip())
            if args.r > 0: q.ref = fields[args.r-1].strip()
            if args.a > 0: q.alt = fields[args.a-1].strip()
            if args.t > 0: q.tpt = fields[args.t-1].strip()
            q.is_codon = False
        else:                   # cDNA
            q = QuerySNV()
            q.op = '|'.join(indices.extract(fields))
            q.tok = fields[args.g-1].strip()
            q.pos = parse_pos(fields[args.n-1].strip())
            if args.r > 0: q.ref = fields[args.r-1].strip()
            if args.a > 0: q.alt = fields[args.a-1].strip()
            if args.t > 0: q.tpt = fields[args.t-1].strip()
            q.is_codon = False

    elif args.g > 0 and args.m > 0: # gene and mutation string

        q = parse_mutation_str(fields[args.m-1].strip(), muttype)
        q.op = '|'.join(indices.extract(fields))
        q.tok = fields[args.g-1].strip()
        if args.t > 0: q.tpt = fields[args.t-1].strip()

    elif args.m > 0:

        q = parse_tok_mutation_str(fields[args.m-1].strip(), muttype)
        q.op = '|'.join(indices.extract(fields))
        if args.t > 0: q.tpt = fields[args.t-1].strip()

    else:
        q.Query()
        q.msg = "InvalidInputLine"
        err_warn("Invalid line: %s" % line.strip('\n'))

    return q

def vcf_parse_mutation(args, at='g'):

    nrec = 0
    for line in opengz(args.vcf):

        if nrec % 1000 == 0:
            sys.stderr.write("\rProcessed %d records\033[K" % nrec)
        nrec += 1

        if line.startswith('##'):
            sys.stdout.write(line)
            continue
        if line.startswith('#CHROM'):
            sys.stdout.write(line.strip()+'\t'+print_header_s()+'\n')
            continue

        fields = line.strip().split('\t')
        chrm = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alts = fields[4]
        for alt in alts.split(','):
            if alt == '<DEL>':
                q = QueryDEL()
                q.tok = chrm
                m = re.match(r'.*END=(\d+)', fields[7])
                q.beg = pos
                q.end = int(m.group(1))
            elif len(ref) == 1 and len(alt) == 1:
                q = QuerySNV()
                q.tok = chrm
                q.pos = pos
                q.ref = ref
                q.alt = alt
            elif len(ref) > 1 and len(alt) == 1 and ref[0] == alt:
                q = QueryDEL()
                q.tok = chrm
                q.beg = pos + 1
                q.end = pos + len(ref) - 1
                q.delseq = ref[1:]
            elif len(ref) == 1 and len(alt) > 1 and ref == alt[0]:
                q = QueryINS()
                q.tok = chrm
                q.pos = pos
                q.insseq = alt[1:]
            elif len(ref) > 1 and len(alt) > 1:
                q = QueryMNV()
                q.tok = chrm
                q.beg = pos
                q.end = pos + len(ref) - 1
                q.refseq = ref.upper()
                q.altseq = alt.upper()
            else:
                q = Query()     # parsing failure
                err_warn("Invalid VCF line: %s" % line.strip('\n'))

            q.op = line.strip()

            yield q, line

    sys.stderr.write("\nProcessed %d records\033[K\n" % nrec)

def list_parse_mutation(args, muttype):

    indices = parse_indices(args.o)
    if args.skipheader:
        args.l.readline()

    for line in args.l:
        # print line.strip()
        if args.d == 's':
            fields = line.strip().split()
        else:
            fields = line.strip().split(args.d)
        try:
            q = _list_parse_mutation(args, fields, indices, muttype)
        # except IndexError as e:
        #     err_print(str(e))
        #     err_print("may be different delimiter in the input file? (try -d s)")
        #     sys.exit(1)
        # except InvalidInputError as e:
        #     err_print(str(e))
        #     continue
        except Exception as e:
            wrap_exception(e, '|'.join(indices.extract(fields)), args)
            continue

        yield q, line

def parser_add_mutation(parser):

    parser.add_argument('--noheader', action='store_true', help='repress header print')
    parser.add_argument('-i', default=None,
                        help='<gene/chrm>:<mutation>, E.g., MET:1010, PIK3CA:E545K, PIK3CA:c.1633G>A, chr12:25398285')
    parser.add_argument('-l', default=None, type = argparse.FileType('r'), help = 'mutation list file')
    parser.add_argument('--vcf', default=None, help = 'vcf input file')
    parser.add_argument('-d', default="\t",
                        help="table delimiter [\\t], use 's' for space.")
    parser.add_argument('-g', type=int,
                        default=-1, help='column for gene/chromosome (1-based)')
    parser.add_argument('-p', type=int, default=-1,
                        help='column for amino acid position (1-based)')
    parser.add_argument('-n', type=int, default=-1,
                        help='column for nucleotide position (1-based)')
    parser.add_argument('-r', type=int, default=-1,
                        help='column for reference base/amino acid (1-based)')
    parser.add_argument('-a', type=int, default=-1,
                        help='column for variant base/amino acid (1-based)')
    parser.add_argument('-t', type=int, default=-1,
                        help='columns for preferred transcript (1-based)')
    parser.add_argument('-m', type=int, default=1,
                        help='column for <gene/chrm>:<mutation> (1-based)')
    parser.add_argument('-o', default='-', 
                        help='columns to be printed in output (1-based), e.g., 3,4,5-10')
    parser.add_argument('--skipheader', action='store_true',
                        help='skip header')
    parser.add_argument('--seqmax', type=int, default=10, help='maximum reference sequence to output (10), use -1 for infinity')
    parser.add_argument('--max-candidates', type=int, dest='nc', default=10, help='maximum candidate output for fuzzy search') 
    parser.add_argument('--oneline', action='store_true', help='output one line for each query')
    parser.add_argument('--aa3', action='store_true', help='use 3 letter code for protein output')
    parser.add_argument('--aacontext', type=int, default=0, help='output amino acid context')
    parser.add_argument('--haplotype', action='store_true', help='use haplotype mode for mnv')
    parser.add_argument('--print-protein', dest='pp', action='store_true', help='print protein sequence')
    parser.add_argument('--print-protein-pretty', dest='ppp', action='store_true', help='print protein sequence in a human-readable format')
