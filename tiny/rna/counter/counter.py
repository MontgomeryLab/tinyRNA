"""tiny-count is a precision counting tool for hierarchical classification and quantification of small RNA-seq reads"""

import multiprocessing as mp
import traceback
import argparse
import shutil
import sys
import os

from collections import defaultdict
from typing import Tuple, List, Dict
from pkg_resources import resource_filename

from tiny.rna.counter.validation import GFFValidator
from tiny.rna.counter.features import Features, FeatureCounter
from tiny.rna.counter.statistics import MergedStatsManager
from tiny.rna.util import report_execution_time, from_here, ReadOnlyDict, get_timestamp, add_transparent_help
from tiny.rna.configuration import CSVReader, PathsFile, get_templates

# Global variables for multiprocessing
counter: FeatureCounter


def get_args():
    """Get input arguments from the user/command line."""

    arg_parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__,
        add_help=False,
    )

    required_args = arg_parser.add_argument_group(
        title="Required arguments",
        description="You must either provide a Paths File or request templates for detailing your configuration.")
    optional_args = arg_parser.add_argument_group(
        title="Optional arguments",
        description="These options can be used in conjunction with the Paths File (-pf) argument mentioned above.")

    # Required arguments
    mutex_top_grp = required_args.add_mutually_exclusive_group(required=True)
    mutex_top_grp.add_argument('-pf', '--paths-file', metavar='FILE', help='your Paths File')
    mutex_top_grp.add_argument('--get-templates', action='store_true',
                               help='Copies the template configuration files required by '
                                    'tiny-count into the current directory.')

    # Optional arguments
    optional_args.add_argument('-o', '--out-prefix', metavar='PREFIX', default='tiny-count_{timestamp}',
                               help='The output prefix to use for file names. All occurrences of the '
                                    'substring {timestamp} will be replaced with the current date and time.')
    optional_args.add_argument('-nh', '--normalize-by-hits', metavar='T/F', default='T',
                               help='If T/true, normalize counts by (selected) overlapping feature counts.')
    optional_args.add_argument('-dc', '--decollapse', action='store_true',
                               help='Create a decollapsed copy of all SAM files listed in your Samples Sheet. '
                                    'This option is ignored for non-collapsed inputs.')
    optional_args.add_argument('-sv', '--stepvector', choices=['Cython', 'HTSeq'], default='Cython',
                               help='Select which StepVector implementation is used to find '
                                    'features overlapping an interval.')
    optional_args.add_argument('-a', '--all-features', action='store_true', help=argparse.SUPPRESS)
                               #help='Represent all features in output counts table, '
                               #     'even if they did not match in Stage 1 selection.')
    optional_args.add_argument('-p', '--is-pipeline', action='store_true',
                               help='Indicates that tiny-count was invoked as part of a pipeline run '
                                    'and that input files should be sourced as such.')
    optional_args.add_argument('-d', '--report-diags', action='store_true',
                               help='Produce diagnostic information about uncounted/eliminated '
                                    'selection elements.')

    add_transparent_help(arg_parser)
    args = arg_parser.parse_args()

    if args.get_templates:
        get_templates("tiny-count")
        sys.exit(0)
    else:
        args_dict = vars(args)
        args_dict['out_prefix'] = args.out_prefix.replace('{timestamp}', get_timestamp())
        args_dict['normalize_by_hits'] = args.normalize_by_hits.lower() in ['t', 'true']
        return ReadOnlyDict(args_dict)


def load_samples(samples_csv: str, is_pipeline: bool) -> List[Dict[str, str]]:
    """Parses the Samples Sheet to determine library names and alignment files for counting

    Sample files may have a .fastq(.gz) extension (i.e. when tiny-count is called as part of a
    pipeline run) or a .sam extension (i.e. when tiny-count is called as a standalone tool).

    Args:
        samples_csv: a csv file which defines sample group, replicate, and file location
        is_pipeline: helps locate sample SAM files. If true, files are assumed to reside
            in the working directory. If false, files are assumed to reside in the same
            directory as their source FASTQ files with '_aligned_seqs.sam' appended
            (i.e. /dir/sample1.fastq -> /dir/sample1_aligned_seqs.sam).

    Returns:
        inputs: a list of dictionaries for each library, where each dictionary defines the
        library name and the location of its SAM file for counting.
    """

    def get_library_filename(csv_row_file: str, samples_csv: str) -> str:
        """The input samples.csv may contain either fastq or sam files"""

        file_ext = os.path.splitext(csv_row_file)[1].lower()

        # If the sample file has a fastq(.gz) extension, infer the name of its pipeline-produced .sam file
        if file_ext in [".fastq", ".gz"]:
            # Fix relative paths to be relative to sample_csv's path, rather than relative to cwd
            csv_row_file = os.path.basename(csv_row_file) if is_pipeline else from_here(samples_csv, csv_row_file)
            csv_row_file = os.path.splitext(csv_row_file)[0] + "_aligned_seqs.sam"
        elif file_ext == ".sam":
            if not os.path.isabs(csv_row_file):
                raise ValueError("The following file must be expressed as an absolute path:\n%s" % (csv_row_file,))
        else:
            raise ValueError("The filenames defined in your Samples Sheet must have a .fastq(.gz) or .sam extension.\n"
                             "The following filename contained neither:\n%s" % (csv_row_file,))
        return csv_row_file

    inputs = list()

    for row in CSVReader(samples_csv, "Samples Sheet").rows():
        library_name = f"{row['Group']}_rep_{row['Replicate']}"
        library_file_name = get_library_filename(row['File'], samples_csv)
        library_normalization = row['Normalization']

        record = {
            "Name": library_name,
            "File": library_file_name,
            "Norm": library_normalization
        }

        if record not in inputs: inputs.append(record)

    return inputs


def load_config(features_csv: str, is_pipeline: bool) -> List[dict]:
    """Parses the Features Sheet to provide inputs to FeatureSelector and ReferenceTables

    Args:
        features_csv: a csv file which defines feature sources and selection rules
        is_pipeline: not currently used; helps properly locate input files

    Returns:
        rules: a list of dictionaries, each representing a parsed row from input.
            Note that these are just rule definitions which FeatureSelector will
            further digest to produce its rules table.
    """

    rules = list()

    for rule in CSVReader(features_csv, "Features Sheet").rows():
        rule['nt5end'] = rule['nt5end'].upper().translate({ord('U'): 'T'})  # Convert RNA base to cDNA base
        rule['Identity'] = (rule.pop('Key'), rule.pop('Value'))             # Create identity tuple
        rule['Hierarchy'] = int(rule['Hierarchy'])                          # Convert hierarchy to number
        rule['Overlap'] = rule['Overlap'].lower()                           # Built later in ReferenceTables

        # Duplicate rule entries are not allowed
        if rule not in rules: rules.append(rule)

    return rules


def load_gff_files(paths: PathsFile, libraries: List[dict], rules: List[dict]) -> Dict[str, list]:
    """Retrieves the appropriate file path and alias attributes for each GFF,
    then validates

    Args:
        paths: a loaded PathsFile with a populated gff_files parameter
        libraries: libraries parsed from Samples Sheet, each as a dict with a 'File' key
        rules: selection rules parsed from Features Sheet

    Returns:
        gff: a dict of GFF files with lists of alias attribute keys
    """

    gff_files = defaultdict(list)
    ignore_alias = ["id"]

    # Build dictionary of unique files and allowed aliases
    for gff in paths['gff_files']:
        gff_files[gff['path']].extend(
            alias for alias in gff.get('alias', ())
            if alias.lower() not in ignore_alias
        )

    # Remove duplicate aliases (per file), keep original order
    for file, alias in gff_files.items():
        gff_files[file] = sorted(set(alias), key=alias.index)

    # Prepare supporting file inputs for GFF validation
    sam_files = [lib['File'] for lib in libraries]

    # Todo: update GFFValidator so that genomes and ebwt can be safely passed here
    # ref_genomes = [map_path(fasta) for fasta in paths['reference_genome_files']]
    # ebwt_prefix = map_path(paths['ebwt'])

    GFFValidator(gff_files, rules, alignments=sam_files).validate()
    return gff_files


@report_execution_time("Counting and merging")
def map_and_reduce(libraries, paths, prefs):
    """Assigns one worker process per library and merges the statistics they report"""

    # MergedStatsManager handles final output files, regardless of multiprocessing status
    summary = MergedStatsManager(Features, paths['features_csv'], prefs)

    # Use a multiprocessing pool if multiple sam files were provided
    if len(libraries) > 1:
        mp.set_start_method("fork")
        with mp.Pool(len(libraries)) as pool:
            async_results = pool.imap_unordered(counter.count_reads, libraries)

            for stats_result in async_results:
                summary.add_library(stats_result)
    else:
        # Only one library, multiprocessing not beneficial for task
        summary.add_library(counter.count_reads(libraries[0]))

    return summary


@report_execution_time("tiny-count's overall runtime")
def main():
    # Get command line arguments.
    args = get_args()
    is_pipeline = args['is_pipeline']

    try:
        # Load and validate config files and input files
        paths = PathsFile(args['paths_file'], is_pipeline)
        libraries = load_samples(paths['samples_csv'], is_pipeline)
        selection = load_config(paths['features_csv'], is_pipeline)
        gff_files = load_gff_files(paths, libraries, selection)

        # global for multiprocessing
        global counter
        counter = FeatureCounter(gff_files, selection, **args)

        # Assign and count features using multiprocessing and merge results
        merged_counts = map_and_reduce(libraries, paths, args)

        # Write final outputs
        merged_counts.write_report_files()
    except Exception as e:
        if type(e) is SystemExit: return
        traceback.print_exception(*sys.exc_info())
        if args['is_pipeline']:
            print("\n\ntiny-count encountered an error. Don't worry! You don't have to start over.\n"
                  "You can resume the pipeline at tiny-count. To do so:\n\t"
                  "1. cd into your Run Directory\n\t"
                  '2. Run "tiny recount --config your_run_config.yml"\n\t'
                  '   (that\'s the processed run config) ^^^\n\n', file=sys.stderr)


if __name__ == '__main__':
    main()
