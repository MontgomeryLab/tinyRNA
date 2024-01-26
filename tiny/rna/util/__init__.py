import functools
import time
import gzip
import re

from .argparse import *
from .formatting import *
from .csv import *


class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


def report_execution_time(step_name: str):
    def timer(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            start = time.time()
            return_val = func(*args, **kwargs)
            end = time.time()
            hrs, rem = divmod(end - start, 3600)
            min, sec = divmod(rem, 60)
            units = [(hrs, '%dh'), (min, '%dm'), (sec, '%.2fs')]
            print(f"{step_name} took {' '.join([u % q for q, u in units if q > 0])}")
            return return_val
        return wrapper
    return timer


def sorted_natural(lines, key=None, reverse=False):
    """Sorts alphanumeric strings with entire numbers considered in the sorting order,
    rather than the default behavior which is to sort by the individual ASCII values
    of the given number. Returns a sorted copy of the list, just like sorted().

    Not sure who to credit... it seems this snippet has been floating around for quite
    some time. Strange that there isn't something in the standard library for this."""

    convert = lambda text: int(text) if text.isdigit() else text.lower()
    extract = (lambda data: key(data)) if key is not None else lambda x: x
    alphanum_key = lambda elem: [convert(c) for c in re.split(r'(\d+)', extract(elem))]
    return sorted(lines, key=alphanum_key, reverse=reverse)


def append_to_exception(e, msg):
    """Appends to an exception's error message while preserving its provenance and traceback"""

    if type(e) is KeyError:
        e.args += (msg,)
    else:
        primary_msg = "%s\n%s" % (str(e.args[0]), msg)
        e.args = (primary_msg,) + e.args[1:]


# File IO interface for reading and writing Gzip files
gzip_open = functools.partial(gzip.GzipFile, compresslevel=6, fileobj=None, mtime=0)
