import re

from datetime import datetime


class ReportFormatter:
    error_header   = "========[ ERROR ]=============================================="
    warning_header = "========[ WARNING ]============================================"

    def __init__(self, key_mapper):
        self.key_mapper = key_mapper
        self.warnings = []
        self.errors = []

    def print_report(self):
        for error_section in self.errors:
            print(self.error_header)
            print(error_section)
            print()

        for warning_section in self.warnings:
            print(self.warning_header)
            print(warning_section)
            print()

    def add_error_section(self, header, body=None):
        self.add_section(header, body, self.errors)

    def add_warning_section(self, header, body=None):
        self.add_section(header, body, self.warnings)

    def add_section(self, header, body, dest):
        if isinstance(body, dict):
            body = self.nested_dict(body)
        elif isinstance(body, list):
            body = '\n'.join(body)
        if body:
            dest.append('\n'.join([header, body]))
        else:
            dest.append(header)

    def nested_dict(self, summary, indent='\t'):
        report_lines = self.recursive_indent(summary, indent)
        return '\n'.join(report_lines)

    def recursive_indent(self, mapping: dict, indent: str):
        """Converts a nested dictionary into a properly indented list of lines

        Args:
            mapping: the nested dictionary to be formatted
            indent: the indent string to prepend to lines on the current level
        """

        lines = []
        for key, val in mapping.items():
            key_header = f"{indent}{self.key_mapper.get(key, key)}: "
            if isinstance(val, dict):
                lines.append(key_header)
                if len(val):
                    lines.extend(self.recursive_indent(val, indent + '\t'))
                else:
                    lines.append(indent + "\t" + "{ NONE }")
            elif isinstance(val, (list, set)):
                lines.append(key_header)
                if len(val):
                    lines.extend([indent + '\t' + line for line in map(str, val)])
                else:
                    lines.append(indent + "\t" + "[ NONE ]")
            else:
                lines.append(key_header + str(val))
        return lines


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


# For timestamp matching and creation
timestamp_format = r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}"
def get_timestamp():
    return datetime.now().strftime('%Y-%m-%d_%H-%M-%S')


def make_filename(args, ext='.csv'):
    return '_'.join([str(chnk) for chnk in args if chnk]) + ext
