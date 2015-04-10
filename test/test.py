#!/usr/bin/env python
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

result = ''
for line in open(sys.argv[1]):

    if line.startswith('$'):

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
