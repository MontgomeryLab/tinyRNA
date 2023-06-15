"""tiny-count is a precision counting tool for hierarchical classification and quantification of small RNA-seq reads"""

import multiprocessing as mp
import traceback
import argparse
import sys
import os

from typing import List, Dict, Tuple

from tiny.rna.counter.validation import GFFValidator, AlignmentSqValidator
from tiny.rna.counter.features import Features, FeatureCounter
from tiny.rna.counter.statistics import MergedStatsManager
from tiny.rna.counter.hts_parsing import ReferenceFeatures, ReferenceSeqs, ReferenceBase
from tiny.rna.configuration import PathsFile, SamplesSheet, FeaturesSheet, get_templates
from tiny.rna.util import (
    report_execution_time,
    add_transparent_help,
    get_timestamp,
    make_filename,
    ReadOnlyDict
)

# Global variables for multiprocessing
counter: FeatureCounter


def get_args() -> 'ReadOnlyDict':
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
    optional_args.add_argument('-o', '--out-prefix', metavar='PREFIX',
                               help='The output prefix to use for file names.')
    optional_args.add_argument('-ng', '--normalize-by-genomic-hits', metavar='T/F', default='T',
                               help='Normalize counts by genomic hits.')
    optional_args.add_argument('-nf', '--normalize-by-feature-hits', metavar='T/F', default='T',
                               help='Normalize counts by feature hits.')
    optional_args.add_argument('-vs', '--verify-stats', metavar='T/F', default='T',
                               help='Verify that all reported stats are internally consistent.')
    optional_args.add_argument('-dc', '--decollapse', action='store_true',
                               help='Create a decollapsed SAM copy of all files listed in your Samples Sheet. '
                                    'This option is ignored for non-collapsed inputs.')
    optional_args.add_argument('-sv', '--stepvector', choices=['Cython', 'HTSeq'], default='Cython',
                               help='Select which StepVector implementation is used to find '
                                    'features overlapping an interval.')
    optional_args.add_argument('-a', '--all-features', action='store_true', help=argparse.SUPPRESS)  # deprecated
    optional_args.add_argument('-p', '--in-pipeline', action='store_true',
                               help='All file inputs and outputs will be read from and written to the working '
                                    'directory regardless of the exact paths listed in configuration files. '
                                    'This is convenient when working with workflow runners.')
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
        prefix, autodoc_dir = get_output_dest(args_dict)
        args_dict['out_prefix'] = prefix
        args_dict['autodoc_dir'] = autodoc_dir
        for tf in ('normalize_by_feature_hits', 'normalize_by_genomic_hits', 'verify_stats'):
            args_dict[tf] = args_dict[tf].lower() in ['t', 'true']
        return ReadOnlyDict(args_dict)


def get_output_dest(args: dict) -> Tuple[str, str]:
    """Creates the run directory and autodoc directory if running in standalone mode"""

    prefix = args['out_prefix']
    if "{timestamp}" in prefix:
        print("Support for {timestamp} substitution in tiny-count "
              "prefixes was removed in v1.6", file=sys.stderr)

    if not args['in_pipeline']:
        run_ts = get_timestamp()
        run_directory = f"tiny-count_{run_ts}"
        prefix_tail = make_filename([prefix, run_ts], ext='')
        prefix = os.path.join(run_directory, prefix_tail)
        config = os.path.join(run_directory, "config")
        os.makedirs(config)
    else:
        config = None  # Handled by workflow runner

    return prefix, config


def load_paths(args: ReadOnlyDict) -> 'PathsFile':
    """Parses the Paths File which provides the locations of config files and reference files.

    Args: A ReadOnlyDict of argparse commandline arguments
    Returns: A loaded PathsFile object
    """

    context = "Pipeline Step" if args['in_pipeline'] else "Standalone Run"
    paths = PathsFile(args['paths_file'], args['in_pipeline'])

    if context == "Standalone Run":
        paths.save_run_profile(args['autodoc_dir'])

    return paths


def load_samples(samples_csv: str, args: ReadOnlyDict) -> List[Dict[str, str]]:
    """Parses the Samples Sheet to determine library names and alignment files for counting

    Sample files can have a .fastq(.gz) extension (i.e. when tiny-count is called as part of a
    pipeline run) or a .sam/.bam extension (i.e. when tiny-count is called as a standalone tool).

    Args:
        samples_csv: a csv file which defines sample group, replicate, and file location
        in_pipeline: helps locate sample alignment files. If true, files are assumed to
            reside in the working directory.

    Returns:
        inputs: a list of dictionaries for each library. Each dictionary holds
        the library's name, file location, and normalization preferences.
    """

    context = "Pipeline Step" if args['in_pipeline'] else "Standalone Run"
    samples = SamplesSheet(samples_csv, context=context)
    inputs = list()

    for file, group_rep, norm in zip(samples.hts_samples, samples.groups_reps, samples.normalizations):
        record = {
            "Name": "_rep_".join(group_rep),
            "File": file,
            "Norm": norm
        }

        if record not in inputs: inputs.append(record)

    return inputs


def load_config(features_csv: str, args: ReadOnlyDict) -> List[dict]:
    """Parses the Features Sheet to provide inputs to FeatureSelector and reference parsers

    Args:
        features_csv: a csv file which defines feature sources and selection rules
        in_pipeline: not currently used; helps properly locate input files

    Returns:
        rules: a list of dictionaries, each representing a parsed row from input.
            Note that these are just rule definitions which FeatureSelector will
            further digest to produce its rules table.
    """

    context = "Pipeline Step" if args['in_pipeline'] else "Standalone Run"
    features = FeaturesSheet(features_csv, context=context)

    return features.rules


def load_references(paths: PathsFile, libraries: List[dict], rules: List[dict], prefs) -> ReferenceBase:
    """Determines the reference source (GFF or alignment @SQ headers) and constructs the appropriate object

    Args:
        paths: a PathsFile object that optionally contains a `gff_files` configuration
        libraries: libraries parsed from Samples Sheet, each as a dict with a 'File' key
        rules: selection rules parsed from Features Sheet
        prefs: command line arguments to pass to the ReferenceBase subclass

    Returns:
        references: a ReferenceBase object, subclassed based on
            the presence of GFF files in the Paths File
    """

    gff_files = paths.get_gff_config()
    aln_files = [lib['File'] for lib in libraries]

    if gff_files:
        GFFValidator(gff_files, rules, alignments=aln_files).validate()
        references = ReferenceFeatures(gff_files, **prefs)
    else:
        sq_validator = AlignmentSqValidator(aln_files)

        # Reuse sq_validator's parsing results to save time
        sq_validator.validate()
        sequences = sq_validator.reference_seqs
        references = ReferenceSeqs(sequences, **prefs)

    return references


@report_execution_time("Counting and merging")
def map_and_reduce(libraries, paths, prefs):
    """Assigns one worker process per library and merges the statistics they report"""

    # MergedStatsManager handles final output files, regardless of multiprocessing status
    summary = MergedStatsManager(Features, paths['features_csv'], prefs)

    # Use a multiprocessing pool if multiple alignment files were provided
    if len(libraries) > 1:
        mp.set_start_method("fork")
        with mp.Pool(len(libraries)) as pool:
            async_results = pool.imap_unordered(counter.count_reads, libraries)

            for result in async_results:
                summary.add_library_stats(result)
    else:
        # Only one library, multiprocessing not beneficial for task
        summary.add_library_stats(counter.count_reads(libraries[0]))

    return summary


@report_execution_time("tiny-count's overall runtime")
def main():
    # Get command line arguments.
    args = get_args()

    try:
        # Load and validate config and input files
        paths = load_paths(args)
        libraries = load_samples(paths['samples_csv'], args)
        selection = load_config(paths['features_csv'], args)
        reference = load_references(paths, libraries, selection, args)

        # global for multiprocessing
        global counter
        counter = FeatureCounter(reference, selection, **args)

        # Assign and count features using multiprocessing and merge results
        merged_counts = map_and_reduce(libraries, paths, args)

        # Write final outputs
        merged_counts.write_report_files()
    except Exception as e:
        if type(e) is SystemExit: return
        traceback.print_exception(*sys.exc_info())
        if args['in_pipeline']:
            print("\n\ntiny-count encountered an error. Don't worry! You don't have to start over.\n"
                  "You can resume the pipeline at tiny-count. To do so:\n\t"
                  "1. cd into your Run Directory\n\t"
                  '2. Run "tiny recount --config your_run_config.yml"\n\t'
                  '   (that\'s the processed run config) ^^^\n\n', file=sys.stderr)


if __name__ == '__main__':
    main()
