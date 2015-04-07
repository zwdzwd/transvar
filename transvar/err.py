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

class UnknownChromosomeError(Exception):
    pass


class SequenceRetrievalError(Exception):
    pass

class WrongReferenceError(Exception):
    pass


def err_die(msg):
    fn = inspect.stack()[1][3]
    sys.stderr.write('[%s] %s\n' % (fn, msg))
    sys.stderr.write('[%s] abort\n' % fn)
    sys.exit(1)

def err_warn(msg):
    fn = inspect.stack()[1][3]
    sys.stderr.write('[%s] warning: %s\n' % (fn, msg))

def err_raise(cls, msg):
    fn = inspect.stack()[1][3]
    raise cls('[%s] exception: %s' % (fn, msg))

def err_print(msg):
    fn = inspect.stack()[1][3]
    sys.stderr.write('[%s] %s\n' % (fn, str(msg)))
