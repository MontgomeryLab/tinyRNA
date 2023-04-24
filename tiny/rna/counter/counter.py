"""tiny-count is a precision counting tool for hierarchical classification and quantification of small RNA-seq reads"""

import multiprocessing as mp
import traceback
import argparse
import sys
import os

from typing import List, Dict

from tiny.rna.counter.validation import GFFValidator, SamSqValidator
from tiny.rna.counter.features import Features, FeatureCounter
from tiny.rna.counter.statistics import MergedStatsManager
from tiny.rna.counter.hts_parsing import ReferenceFeatures, ReferenceSeqs, ReferenceBase
from tiny.rna.configuration import PathsFile, SamplesSheet, CSVReader, get_templates
from tiny.rna.util import (
    report_execution_time,
    add_transparent_help,
    append_to_exception,
    get_timestamp,
    ReadOnlyDict
)

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
    optional_args.add_argument('-a', '--all-features', action='store_true', help=argparse.SUPPRESS)
                               #help='Represent all features in output counts table, '
                               #     'even if they did not match in Stage 1 selection.')
    optional_args.add_argument('-p', '--in-pipeline', action='store_true',
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
        for tf in ('normalize_by_feature_hits', 'normalize_by_genomic_hits', 'verify_stats'):
            args_dict[tf] = args_dict[tf].lower() in ['t', 'true']
        return ReadOnlyDict(args_dict)


def load_samples(samples_csv: str, in_pipeline: bool) -> List[Dict[str, str]]:
    """Parses the Samples Sheet to determine library names and alignment files for counting

    Sample files may have a .fastq(.gz) extension (i.e. when tiny-count is called as part of a
    pipeline run) or a .sam/.bam extension (i.e. when tiny-count is called as a standalone tool).

    Args:
        samples_csv: a csv file which defines sample group, replicate, and file location
        in_pipeline: helps locate sample SAM files. If true, files are assumed to reside
            in the working directory. If false, files are assumed to reside in the same
            directory as their source FASTQ files with '_aligned_seqs.sam' appended
            (i.e. /dir/sample1.fastq -> /dir/sample1_aligned_seqs.sam).

    Returns:
        inputs: a list of dictionaries for each library, where each dictionary defines the
        library name and the location of its SAM file for counting.
    """

    context = "Pipeline Step" if in_pipeline else "Standalone Run"
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


def load_config(features_csv: str, in_pipeline: bool) -> List[dict]:
    """Parses the Features Sheet to provide inputs to FeatureSelector and reference parsers

    Args:
        features_csv: a csv file which defines feature sources and selection rules
        in_pipeline: not currently used; helps properly locate input files

    Returns:
        rules: a list of dictionaries, each representing a parsed row from input.
            Note that these are just rule definitions which FeatureSelector will
            further digest to produce its rules table.
    """

    sheet = CSVReader(features_csv, "Features Sheet")
    rules = list()

    try:
        for rule in sheet.rows():
            rule['nt5end'] = rule['nt5end'].upper().translate({ord('U'): 'T'})  # Convert RNA base to cDNA base
            rule['Identity'] = (rule.pop('Key'), rule.pop('Value'))             # Create identity tuple
            rule['Hierarchy'] = int(rule['Hierarchy'])                          # Convert hierarchy to number
            rule['Overlap'] = rule['Overlap'].lower()                           # Built later in reference parsers

            # Duplicate rule entries are not allowed
            if rule not in rules: rules.append(rule)
    except Exception as e:
        msg = f"Error occurred on line {sheet.line_num} of {os.path.basename(features_csv)}"
        append_to_exception(e, msg)
        raise

    return rules


def load_references(paths: PathsFile, libraries: List[dict], rules: List[dict], prefs) -> ReferenceBase:
    """Determines the reference source (GFF or SAM @SQ headers) and constructs the appropriate object

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
    sam_files = [lib['File'] for lib in libraries]

    if gff_files:
        GFFValidator(gff_files, rules, alignments=sam_files).validate()
        references = ReferenceFeatures(gff_files, **prefs)
    else:
        sq_validator = SamSqValidator(sam_files)

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

    # Use a multiprocessing pool if multiple sam files were provided
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
    in_pipeline = args['in_pipeline']

    try:
        # Load and validate config files and input files
        paths = PathsFile(args['paths_file'], in_pipeline)
        libraries = load_samples(paths['samples_csv'], in_pipeline)
        selection = load_config(paths['features_csv'], in_pipeline)
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
        if in_pipeline:
            print("\n\ntiny-count encountered an error. Don't worry! You don't have to start over.\n"
                  "You can resume the pipeline at tiny-count. To do so:\n\t"
                  "1. cd into your Run Directory\n\t"
                  '2. Run "tiny recount --config your_run_config.yml"\n\t'
                  '   (that\'s the processed run config) ^^^\n\n', file=sys.stderr)


if __name__ == '__main__':
    main()
