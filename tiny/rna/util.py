import argparse
import functools
import textwrap
import time
import os
import re


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
    # Properly formats argparse helpstring fields that contain newlines
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            lines = []
            for line in text[2:].splitlines():
                fill = textwrap.fill(line, width, subsequent_indent="  ")
                lines.extend(fill.splitlines())
            return lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def from_here(config_file, input_file):
    """Calculates paths relative to the config file which contains them"""
    config_file, input_file = (os.path.expanduser(p) for p in [config_file, input_file])

    if not os.path.isabs(input_file):
        from_here = os.path.dirname(config_file)
        input_file = os.path.normpath(os.path.join(from_here, input_file))

    return input_file


def make_filename(args, ext='.csv'):
    return '_'.join([str(chnk) for chnk in args if chnk is not None]) + ext


def get_r_safename(name: str) -> str:
    """Converts a string to a syntactically valid R name

    This can be used to match names along axes of DataFrames produced by R,
    assuming that the R script takes no measures to preserve names itself.
    https://stat.ethz.ch/R-manual/R-devel/library/base/html/make.names.html
    """

    leading_char = lambda x: re.sub(r"^(?=[^a-zA-Z.]+|\.\d)", "X", x)
    special_char = lambda x: re.sub(r"[^a-zA-Z0-9_.]", ".", x)
    return special_char(leading_char(name))