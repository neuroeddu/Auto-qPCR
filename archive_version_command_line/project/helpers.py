import string
from pprint import PrettyPrinter

LOG_QUIET   = 0
LOG_INFO    = 1
LOG_VERBOSE = 5
LOG_DEBUG   = 9

verbosity = LOG_QUIET
pp = PrettyPrinter()

# Useful regular expression patterns
RXP = {
    #@<number> denotes a placeholder that must be replaced by an actual value before use
    'FileExt': r'\.@1$',
    'IntRange': r'^[0-9]+\-[0-9]+$'
}

def log(vlevel, msg, stay_on_line = False):
    if (verbosity >= 0) and (vlevel <= verbosity):
        if stay_on_line:
            print(msg, end='', flush=True)
        else:
            print(msg)

def bprint(msg, blen = None):
    l = len(msg) if blen is None else blen
    if l > 0:
        print("\n{}\n{}\n{}".format('=' * l, msg, '-' * l))
