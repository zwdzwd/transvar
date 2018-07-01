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

import inspect
import sys

class InvalidInputError(Exception):
    pass

class IncompatibleTranscriptError(Exception):
    pass

class UnImplementedError(Exception):
    pass

class UndeterminedStopCodon(Exception):
    pass

class ReferenceUnavailableError(Exception):
    pass

class SequenceRetrievalError(Exception):
    pass

class WrongReferenceError(Exception):
    pass


def err_die(msg):
    fn = inspect.stack()[1][3]
    for line in msg.splitlines():
        sys.stderr.write('[%s] %s\n' % (fn, line))
    sys.stderr.write('[%s] abort\n' % fn)
    sys.exit(1)

def err_warn(msg):
    fn = inspect.stack()[1][3]
    sys.stderr.write('\r[%s] warning: %s\n' % (fn, msg))

def err_raise(cls, msg):
    fn = inspect.stack()[1][3]
    raise cls('[%s] exception: %s' % (fn, msg))

def err_print(msg):
    fn = inspect.stack()[1][3]
    sys.stderr.write('[%s] %s\n' % (fn, str(msg)))

