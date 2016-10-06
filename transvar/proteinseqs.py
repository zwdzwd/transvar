"""
The MIT License

Copyright (c) 2016
Wanding Zhou

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

from .utils import *

def variant_protein_seq_sub(r, t, args, taa_beg, taa_end, taa_alt):

    """ r.taa_alt, r.taa_ref and r.taa_pos must be accessible """
    if args.pp or args.ppp:
        pp = list(aaf(t.get_proteinseq(), args, use_list=True))

        if aa_has_stop(taa_alt):
            if args.ppp:
                pp[taa_beg-1] = '__[%s>%s]' % (''.join(pp[taa_beg:]), taa_alt)
                pp = pp[:taa_beg]
            else:
                pp[taa_beg-1] = taa_alt
                pp = pp[:taa_beg]
        else:
            if args.ppp:
                pp[taa_beg-1] = '__[%s>%s]__' % (
                    ''.join(pp[taa_beg-1:taa_end]), taa_alt)
                del pp[taa_beg:taa_end]
            else:
                pp[taa_beg-1] = taa_alt
                del pp[taa_beg:taa_end]
                
        r.append_info('variant_protein_seq=%s' % ''.join(pp))

def variant_protein_seq_del(r, t, args, taa_beg, taa_end):

    if args.pp or args.ppp:
        pp = list(aaf(t.get_proteinseq(), args, use_list=True))
        if args.pp:
            del pp[taa_beg-1:taa_end]
        elif args.ppp:
            delseq = ''.join(pp[taa_beg-1:taa_end])
            del pp[taa_beg-1:taa_end]
            pp.insert(taa_beg-1, '__[%s_deletion]__' % delseq)
        r.append_info('variant_protein_seq=%s' % ''.join(pp))
        
def variant_protein_seq_fs(r, t, aae, args):

    if args.pp or args.ppp:
        pp = list(aaf(t.get_proteinseq(), args, use_list=True))
        if args.ppp:
            pp = pp[:aae.taa_pos-1] + \
                 ['__[frameshift_%s>%s]' % (
                     ''.join(pp[aae.taa_pos-1:]), aae.new_aa_seq)]
        else:
            pp = pp[:aae.taa_pos]
            pp.append(aaf(aae.new_aa_seq, args))

        r.append_info('variant_protein_seq=%s' % ''.join(pp))

def variant_protein_seq_ins(r, t, args, taa_pos, insseq):

    if args.pp or args.ppp:
        pp = list(aaf(t.get_proteinseq(), args, use_list=True))
        if args.ppp:
            pp = pp[:taa_pos]+['__[insert_%s]__' % aaf(insseq, args)]+pp[taa_pos:]
        else:
            pp = pp[:taa_pos]+list(aaf(insseq, args, use_list=True))+pp[taa_pos:]
        r.append_info('variant_protein_seq=%s' % ''.join(pp))
