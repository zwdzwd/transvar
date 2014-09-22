"""
annotate nucleotide position(s) or mutations
"""
import sys, argparse, re
from transcripts import *
from utils import *
from mutation import parse_tok_mutation_str
from record import *

def pos2codon(thash, chrm, pos):

    for t in thash.get_transcripts(chrm, pos):
        c = t.npos2codon(chrm, pos)
        yield t, c

def pos2codon_longest(thash, chrm, pos):

    longest = None
    longest_c = None
    for t, c in pos2codon(thash,chrm,pos):
        if (not longest) or len(longest) < len(t):
            longest = t
            longest_c = c

    if longest:
        yield longest, longest_c

def _main_core_(args, thash, q):

    if args.longest:
        tc_iter = pos2codon_longest(thash, q.chrm, q.pos)
    else:
        tc_iter = pos2codon(thash, q.chrm, q.pos)

    found = False
    for t, c in tc_iter:
        if isinstance(c, Codon):
            found = True

            r = Record()
            r.chrm = t.chrm
            r.tname = t.name
            r.reg = '%s (%s, coding)' % (t.gene.name, t.strand)
            r.pos = '-'.join(map(str, c.locs))

            # at the ends of retained intron transcripts from ENSEMBL,
            # codon sequence is not always of length 3
            if c.seq in standard_codon_table:
                r.taa_ref = standard_codon_table[c.seq]
            r.taa_pos = c.index
            if q.alt:
                if c.strand == "+":
                    alt_seq = set_seq(c.seq, c.locs.index(q.pos), q.alt)
                else:
                    alt_seq = set_seq(c.seq, 2-c.locs.index(q.pos), complement(q.alt))
                r.taa_alt = standard_codon_table[alt_seq]

            r.gnuc_pos = q.pos
            r.gnuc_ref = c.refseq()[c.locs.index(q.pos)]
            r.gnuc_alt = q.alt

            if c.strand == '+':
                r.tnuc_ref = r.gnuc_ref
                r.tnuc_alt = r.gnuc_alt
                r.tnuc_pos = (c.index-1)*3 + c.locs.index(q.pos) + 1
            else:
                r.tnuc_ref = complement(r.gnuc_ref)
                r.tnuc_alt = complement(r.gnuc_alt) if r.gnuc_alt else ''
                r.tnuc_pos = c.index*3 - c.locs.index(q.pos)

            r.format(q.op)

        elif isinstance(c, NonCoding):
            found = True

            r = Record()
            r.chrm = t.chrm
            r.gnuc_pos = q.pos
            r.tname = t.name
            r.reg = '%s (%s noncoding)' % (t.gene.name, t.strand)
            r.info = c.format()
            r.format(q.op)

    if not found:
        r = Record()
        r.gnuc_ref = q.ref
        r.gnuc_alt = q.alt
        r.gnuc_pos = q.pos
        r.info = 'status=NoValidTranscriptFound'
        r.format(q.op)

def list_parse_genomic_mutation(args):

    indices = parse_indices(args.o)
    if args.skipheader:
        args.l.readline()

    for line in args.l:
        fields = line.strip().split(args.d)
        q = Query()
        q.op = '\t'.join(indices.extract(fields))
        if args.c > 0 and args.p > 0:
            q.chrm = fields[args.c-1].strip()
            q.pos = int(fields[args.p-1].strip())
            if args.r > 0: q.ref = fields[args.r-1].strip()
            if args.v > 0: q.alt = fields[args.v-1].strip()
            if args.t > 0: q.tpt = fields[args.t-1].strip()
            yield q
        elif args.c > 0 and args.m > 0:
            q.chrm = fields[args.c-1].strip()
            if args.t > 0: q.tpt = fields[args.t-1].strip()
            ret = parse_mutation_str(fields[args.m-1].strip())
            if ret:
                q.chrm, q.is_codon, q.pos, q.ref, q.alt = ret
                yield q

        elif args.m > 0:
            ret = parse_tok_mutation_str(fields[args.m-1].strip())
            if ret:
                q.chrm, q.is_codon, q.pos, q.ref, q.alt = ret
                if args.t > 0: q.tpt = fields[args.t-1].strip()
                yield q

def main_list(args, thash):

    for q in list_parse_genomic_mutation(args):
        _main_core_(args, thash, q)

def main_one(args, thash):
    q = Query()
    ret = parse_tok_mutation_str(args.i)
    if not ret: return
    q.op = args.i
    q.chrm, q.is_codon, q.pos, q.ref, q.alt = ret
    _main_core_(args, thash, q)

def oldmainone():
    m = re.match(r'([^:]*):([ATGC]?)(\d+)([ATGC]?)', args.npos)
    chrm = m.group(1)
    ref = m.group(2)
    pos = int(m.group(3))
    alt = m.group(4)

    tpts = thash.get_transcripts(chrm, pos, args.standard)
    if not tpts:
        return

    if not args.alltrans:
        tpt = tpts[0]
        if len(tpts) > 1:
            for _tpt in tpts:
                if _tpt.is_standard():
                    tpt = _tpt
                    break
        tpts = [tpt]

    for tpt in tpts:
        prnstr = args.npos
        c = tpt.npos2codon(chrm, pos)
        if isinstance(c, Codon):
            if alt:
                prnstr += '\t'
                prnstr += nuc_mutation(c, pos, ref, alt)
            else:
                prnstr += '\t'
                prnstr += c.format()
        elif isinstance(c, NonCoding):
            prnstr += '\t'+c.format()
            if alt:
                prnstr += "\t%s\t%s" % (ref, alt)

        print prnstr


def main(args):

    name2gene, thash = parse_annotation(args)

    if args.l:
        main_list(args, thash)

    if args.i:
        main_one(args, thash)
    
def add_parser_anno(subparsers, d):

    parser = subparsers.add_parser('anno', help=__doc__)
    parser_add_annotation(parser, d)
    parser.add_argument('-i', default=None,
                        help="<chrm>:[<ref>]<pos>[<alt>], E.g., chr12:25398285")
    parser.add_argument('-l', default=None,
                        type = argparse.FileType('r'),
                        help='mutation list file')
    parser.add_argument('-d', default="\t", 
                        help="table delimiter of mutation list [\\t]")
    parser.add_argument('-c', type=int, default=-1,
                        help='column for chromosome (1-based)')
    parser.add_argument("-p", type=int, default=-1,
                        help='column for position (1-based)')
    parser.add_argument('-r', type=int, default=-1,
                        help='column for reference base (1-based)')
    parser.add_argument('-v', type=int, default=-1,
                        help='column for variant base (1-based)')
    parser.add_argument('-t', type=int, default=-1,
                        help='columns for preferred transcript (1-based)')
    parser.add_argument("-m", type=int, default=-1,
                        help="column for mutation in format <chrm>:[<ref>]<pos>[<alt>] (1-based)")
    parser.add_argument('-o', default='-',
                        help='columns in the original table to be output (1-based)')
    parser.add_argument('--skipheader', action='store_true',
                        help='skip header')
    parser.add_argument('--longest', action="store_true", help='use longest transcript instead of reporting all transcripts')

    parser.set_defaults(func=main)
