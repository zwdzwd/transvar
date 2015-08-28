#!/usr/bin/env python
"""
The MIT License

Copyright (c) 2015 by The University of Texas MD Anderson Cancer Center
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
import os, sys, re

from subprocess import Popen, PIPE

def wrap(line):

    lw = 80
    fields = line.strip().split('\t')
    yield '\t'.join(fields[:4])
    yield ' '*3+'\t'.join(fields[4:6])
    line = '\t'.join(fields[6:])

    while len(line) > 0:
        k = lw-3
        yield ' '*3+line[:k]
        line = line[k:]

fh = open(sys.argv[2], 'w') # open('README.md.temp', 'w')

import transvar

result = ''
for line in open(sys.argv[1]):

    line = line.replace('@VERSION', transvar.__version__)

    if line.startswith('$'):    # exexcute sentinel

        fh.write(line)
        line = line[2:].strip()
        ar = re.split(r'[\'"]', line.strip())

        B = []
        for i, e in enumerate(ar):
            if i % 2 == 0:
                if len(e.strip())>0:
                    B.extend(e.strip().split())
            else:
                B.append(e.strip())

        p = Popen(B, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        result, err = p.communicate()
        result = result.strip()

        print
        print '======'+line+'======'
        print '\n'.join([rr for rr in result.split('\n') if not (len(rr.strip()) == 0 or rr.startswith('[') or rr.startswith('input'))])

    elif line.startswith('@'):

        for rr in result.split('\n'):

            if len(rr.strip()) == 0 or rr.startswith('[') or rr.startswith('input'):
                continue

            for r in wrap(rr):
                fh.write(r+'\n')

    else:

        fh.write(line)
