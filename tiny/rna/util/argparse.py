import argparse
import textwrap
import re


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


class ReadOnlyDict(dict):
    """A very simple "read-only" wrapper for dictionaries. This will primarily be used
    for passing around argparse's command line args as a dictionary, but ensuring that
    preferences cannot be accidentally modified as they are passed around to a menagerie
    of class constructors. All methods except __setitem__() are deferred to the base class."""

    def __init__(self, rw_dict):
        super().__init__(rw_dict)

    def __setitem__(self, *_):
        raise RuntimeError("Attempted to modify read-only dictionary after construction.")


def add_transparent_help(parser):
    # Allows -h/--help without listing it within the optional arguments section
    parser.add_argument('-h', '--help', action="help", help=argparse.SUPPRESS)