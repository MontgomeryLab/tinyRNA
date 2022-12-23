import argparse
import functools
import textwrap
import gzip
import time
import os
import re

from datetime import datetime


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


class SmartFormatter(argparse.HelpFormatter):
    """Properly formats argparse help string fields that contain newlines.

    Newlines are only preserved in Argparse description fields if they are
    followed by whitespace (e.g. \n\n, \n\t, etc.). This is because docstrings
    can be lengthy and require single newlines to maintain an appropriate width,
    but in the terminal we want this text to instead fill the terminal width.
    Leading whitespace (indentation) is also preserved.

    Newlines are only preserved in argument help fields that begin with the
    token "R|". Otherwise, all newlines are replaced just as they are in
    vanilla argparse."""

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            lines = []
            for line in text[2:].splitlines():
                fill = textwrap.fill(line, width, subsequent_indent="  ")
                lines.extend(fill.splitlines())
            return lines
        return argparse.HelpFormatter._split_lines(self, text, width)

    def _fill_text(self, text, width, indent):
        input_sections = text.split('\n\n')
        output_lines = []

        for section in input_sections:
            section_lines = []
            for subsection in re.split(r'\n(?=\s)', section):
                leading_ws = re.match(r'\s+', subsection)
                sub_indent = indent if not leading_ws else indent + leading_ws.group()
                subsection = textwrap.TextWrapper.wordsep_simple_re.sub(' ', subsection)
                section_lines.append(
                    textwrap.fill(
                        subsection, width,
                        initial_indent=sub_indent,
                        subsequent_indent=sub_indent,
                        replace_whitespace=False,
                        break_on_hyphens=False,
                ))

            output_lines.append('\n'.join(section_lines))

        return '\n\n'.join(output_lines)


def from_here(config_file, input_file):
    """Calculates paths relative to the config file which contains them"""
    config_file, input_file = (os.path.expanduser(p) for p in [config_file, input_file])

    if not os.path.isabs(input_file):
        from_here = os.path.dirname(config_file)
        input_file = os.path.normpath(os.path.join(from_here, input_file))

    return input_file


def make_filename(args, ext='.csv'):
    return '_'.join([str(chnk) for chnk in args if chnk is not None]) + ext


r_reserved_keywords = [
    "if", "else", "repeat", "while", "function",
    "for", "in", "next", "break", "TRUE", "FALSE",
    "NULL", "Inf", "NaN", "NA", "NA_integer_",
    "NA_real_", "NA_complex_", "NA_character_"]


def get_r_safename(name: str) -> str:
    """Converts a string to a syntactically valid R name

    This can be used as the Python equivalent of R's make.names() function.
    https://stat.ethz.ch/R-manual/R-devel/library/base/html/make.names.html
    """

    # If the name starts with a non-letter character or a dot
    # followed by a number, the character "X" is prepended
    leading_char = lambda x: re.sub(r"^(?=[^a-zA-Z.]+|\.\d)", "X", x)

    # If the name contains characters that aren't (locale based) letters,
    # numbers, dot, or underscore, the characters are replaced with a dot
    special_char = lambda x: re.sub(r"[^\w.]", ".", x)

    # If the name contains R keywords, a dot is appended to the keyword
    reserved = "|".join(r_reserved_keywords)
    reserved_wrd = lambda x: re.sub(fr"^({reserved})$", r'\1.', x)

    return reserved_wrd(special_char(leading_char(name)))


class ReadOnlyDict(dict):
    """A very simple "read-only" wrapper for dictionaries. This will primarily be used
    for passing around argparse's command line args as a dictionary, but ensuring that
    preferences cannot be accidentally modified as they are passed around to a menagerie
    of class constructors. All methods except __setitem__() are deferred to the base class."""

    def __init__(self, rw_dict):
        super().__init__(rw_dict)

    def __setitem__(self, *_):
        raise RuntimeError("Attempted to modify read-only dictionary after construction.")

def sorted_natural(lines, reverse=False):
    """Sorts alphanumeric strings with entire numbers considered in the sorting order,
    rather than the default behavior which is to sort by the individual ASCII values
    of the given number. Returns a sorted copy of the list, just like sorted().

    Not sure who to credit... it seems this snippet has been floating around for quite
    some time. Strange that there isn't something in the standard library for this."""

    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split(r'(\d+)', key)]
    return sorted(lines, key=alphanum_key, reverse=reverse)


# File IO interface for reading and writing Gzip files
gzip_open = functools.partial(gzip.GzipFile, compresslevel=6, fileobj=None, mtime=0)


# For timestamp matching and creation
timestamp_format = re.compile(r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}")
def get_timestamp():
    return datetime.now().strftime('%Y-%m-%d_%H-%M-%S')


# Allows -h/--help without listing it within the optional arguments section
def add_transparent_help(parser):
    parser.add_argument('-h', '--help', action="help", help=argparse.SUPPRESS)

