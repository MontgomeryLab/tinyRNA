import ruamel.yaml
import argparse
import shutil
import errno
import time
import copy
import sys
import csv
import re
import os

from pkg_resources import resource_filename
from collections import Counter, OrderedDict, defaultdict
from typing import Union, Any, Dict, Optional
from glob import glob

from tiny.rna.compatibility import RunConfigCompatibility
from tiny.rna.counter.validation import GFFValidator
from tiny.rna.util import get_timestamp, get_csv_dialect, get_r_safename, append_to_exception


class ConfigBase:
    """Base class for YAML-based tinyRNA configurations

    Args:
        config_file: the path of the YAML config file
        BackwardCompatibility: a class for providing compatibility with old config files

    Attributes:
        yaml: the YAML interface for reading config and writing processed config
        config: the configuration object produced by loading the config file
        inf: the input configuration file name
        basename: the basename of config_file
        dir: parent directory of the input file
        dt: a date-time string to use in output prefixes
    """

    def __init__(self, config_file: str, BackwardCompatibility = None):
        if '~' in config_file: config_file = os.path.expanduser(config_file)
        self.dir = os.path.dirname(os.path.abspath(config_file))
        self.inf = config_file
        self.basename = os.path.basename(config_file)
        self.dt = ''

        self.yaml = ruamel.yaml.YAML()
        with open(config_file, 'r') as f:
            self.config = self.yaml.load(f)
            if BackwardCompatibility is not None:
                self.config = BackwardCompatibility(self.config).upgrade()

    def __getitem__(self, key: str) -> Any:
        return self.get(key)

    def __setitem__(self, key: str, val: Union[str, list, dict, bool, None]) -> Union[str, list, dict, bool, None]:
        return self.set(key, val)

    def __contains__(self, key: str) -> bool:
        return key in self.config

    def get(self, key: str, default=None) -> Any:
        return self.config.get(key, default)

    def set(self, key: str, val: Union[str, list, dict, bool]) -> Union[str, list, dict, bool]:
        self.config[key] = val
        return val

    def set_if_not(self, key: str, val: Union[str, list, dict, bool]) -> Any:
        """Apply the setting if it has not been previously set. This differs from
        dict.setdefault() in that existing keys are overwritten if the associated
        value is false in boolean context (e.g. None, False, empty container, etc.)"""

        if not self[key]:
            self[key] = val
            return val
        else:
            return self[key]

    def set_default_dict(self, setting_dict: dict) -> None:
        """Apply each setting in setting_dict if it has not been previously set"""
        for key, val in setting_dict.items():
            # Can't use config.setdefault(), it considers None and [] "already set"
            self.set_if_not(key, val)

    def append_to(self, key: str, val: Any) -> list:
        """Append a list-type setting (per-library settings)"""
        target = self[key]
        if type(target) is list:
            target.append(val)
            return target
        else:
            raise ValueError(f"Tried appending to a non-existent key: {key}")

    def append_if_absent(self, key: str, val: Any) -> list:
        """Append to list-type setting if the value is not already present"""
        target = self[key]
        if val not in target:
            return self.append_to(key, val)

    """========== HELPERS =========="""

    @staticmethod
    def prefix(path: str) -> str:
        """Returns everything from path except the file extension"""
        return os.path.splitext(path)[0]

    @staticmethod
    def cwl_file(file: Union[str,dict], verify=True) -> dict:
        """Returns a minimal File object as defined by CWL"""
        if isinstance(file, dict): file = file['path']
        if '~' in file: file = os.path.expanduser(file)
        if verify and not os.path.exists(file):
            raise(FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file))
        return {'class': 'File', 'path': file}

    @staticmethod
    def cwl_dir(dir: str) -> dict:
        """Returns a minimal Directory object as defined by CWL"""
        return {'class': 'Directory', 'path': dir}

    @staticmethod
    def joinpath(path1: str, path2: str) -> str:
        """Combines two relative paths intelligently"""
        path1 = os.path.expanduser(path1) if path1 is not None else ""
        path2 = os.path.expanduser(path2) if path2 is not None else ""
        if os.path.isabs(path2): return path2
        return os.path.normpath(os.path.join(path1, path2))

    def from_here(self, destination: Union[str,dict,None], origin: Union[str,dict,None] = None) -> Union[str,dict,None]:
        """Calculates paths relative to the input config file. An empty destination returns an empty path.
        If destination is an absolute path, then destination is returned without any further joining.
        Args:
            destination: the path to evaluate relative to the config file. Can be a string or CWL file dictionary.
            origin: the starting point (default: the config file's directory)
        """

        if isinstance(origin, dict): origin = origin.get('path')
        if origin is None: origin = self.dir

        if self.is_path_dict(destination):
            return dict(destination, path=self.joinpath(origin, destination["path"]))
        elif self.is_path_str(destination):
            return self.joinpath(origin, destination)
        else:
            return destination

    @staticmethod
    def is_path_dict(val, empty_ok=False):
        return (isinstance(val, dict) and
                isinstance(val.get("path"), (str, bytes)) and
                (len(val['path']) or empty_ok))

    @staticmethod
    def is_path_str(val, empty_ok=False):
        return (isinstance(val, (str, bytes)) and
                (len(val) or empty_ok))

    @classmethod
    def is_path(cls, val, empty_ok=False):
        return (cls.is_path_dict(val, empty_ok) or
                cls.is_path_str(val, empty_ok))

    def setup_step_inputs(self):
        """For now, only tiny-plot requires additional setup for step inputs
        This function is called at both startup and resume"""

        def setup_tiny_plot_inputs():
            cs_filter = 'plot_class_scatter_filter'
            style_req = ['include', 'exclude']
            classes = (self.get(cs_filter) or {}).get('classes')  # backward compatibility
            if not classes: return

            # Validate filter style
            style = self[cs_filter]['style'].lower()
            assert style in style_req, \
                f'{cs_filter} -> style: must be {" or ".join(style_req)}.'

            # Assign the workflow key and reset the other filter(s)
            self[f"{cs_filter}_{style}"] = classes.copy()
            style_req.remove(style)
            for style in style_req:
                self[f"{cs_filter}_{style}"] = []

        setup_tiny_plot_inputs()

    def create_run_directory(self) -> str:
        """Create the destination directory for pipeline outputs"""
        run_dir = self["run_directory"]
        if not os.path.isdir(run_dir):
            os.makedirs(run_dir)

        return run_dir

    def get_outfile_path(self, infile: str = None) -> str:
        """Returns the path and file for the processed run config"""
        if infile is None: infile = self.inf
        return self.joinpath(self['run_directory'], os.path.basename(infile))

    def write_processed_config(self, filename: str = None) -> str:
        """Writes the current configuration to disk"""
        if filename is None:
            filename = self.get_outfile_path(self.inf)

        if "processed_run_config" in self:
            # The CWL specification doesn't provide a way to get this info,
            # but it's needed for run_directory/config, so we store it here
            self['processed_run_config'] = self.cwl_file(filename, verify=False)

        with open(filename, 'w') as outconf:
            self.yaml.dump(self.config, outconf)

        return filename


class Configuration(ConfigBase):
    """A class for processing and updating a YAML config file for CWL

    Ultimately, this class populates workflow settings and per-library settings. This
    is a convenience to the user as it is tedious to define inputs and outputs pertaining
    to each workflow step in a manner compatible with CWL. Settings are determined by the
    Run Config, Paths File, Samples Sheet, and Features Sheet.
    Users can provide both relative and absolute paths.

    IMPORTANT: Paths provided in any config file are evaluated relative to the containing config file.

    Args:
        config_file: The Run Config file path
        validate_gffs: If true, validate GFFs (not all contexts need to)
        skip_setup: If true, only load the Run Config and Paths File with
            no further processing (useful for testing)

    Attributes:
        paths: the configuration object from processing the paths_config file.
            This holds path info for other config files and prefixes, and is updated
            appropriately if 'run_bowtie_index' is set to 'true'
    """

    def __init__(self, config_file: str, validate_gffs=False, skip_setup=False):
        # Parse YAML configuration file
        super().__init__(config_file, RunConfigCompatibility)

        self.paths = self.load_paths_config()
        self.samples_sheet = self.load_samples_config()
        self.features_sheet = self.load_features_config()

        if skip_setup: return
        self.setup_pipeline()
        self.setup_ebwt_idx()
        self.setup_step_inputs()
        if validate_gffs: self.validate_inputs()

    def load_paths_config(self) -> 'PathsFile':
        """Returns a PathsFile object and updates keys related to the Paths File path"""

        # Resolve relative path to the Paths File and construct
        resolved = self.from_here(self['paths_config'])
        paths = PathsFile(resolved)

        # Absorb PathsFile object keys into the configuration
        for key in [*PathsFile.single, *PathsFile.groups]:
            self[key] = paths.as_cwl_file_obj(key)
        for key in PathsFile.prefix:
            self[key] = paths[key]
        return paths

    def load_samples_config(self) -> 'SamplesSheet':
        samples_sheet_path = self.paths['samples_csv']
        samples_sheet = SamplesSheet(samples_sheet_path, context="Pipeline Start")

        self['sample_basenames'] = samples_sheet.sample_basenames
        self['control_condition'] = samples_sheet.control_condition
        self['run_deseq'] = samples_sheet.is_compatible_df

        self['in_fq'] = [self.cwl_file(fq, verify=False) for fq in samples_sheet.hts_samples]
        self['fastp_report_titles'] = [f"{g}_rep_{r}" for g, r in samples_sheet.groups_reps]

        return samples_sheet

    def load_features_config(self) -> 'FeaturesSheet':
        """Retrieves GFF Source and Type Filter definitions for use in GFFValidator"""

        features_sheet_path = self.paths['features_csv']
        return FeaturesSheet(features_sheet_path, context="Pipeline Start")
            
    def setup_pipeline(self):
        """Overall settings for the whole pipeline"""

        self.dt = get_timestamp()
        self['run_date'], self['run_time'] = self.dt.split('_')

        # Ensure compatible string joins while preserving 0
        for key in ('user', 'run_name', 'run_directory'):
            self[key] = str(self[key]) if self[key] is not None else ''

        default_run_name = '_'.join(x for x in [self['user'], "tinyrna"] if x)
        self['run_name'] = f"{self['run_name'] or default_run_name}_{self.dt}"

        # Prefix Run Directory basename while preserving subdirectory structure
        rd_head, rd_tail = os.path.split(self['run_directory'].rstrip(os.sep))
        basename = '_'.join(x for x in [self['run_name'], rd_tail] if x)
        self['run_directory'] = self.joinpath(rd_head or os.getcwd(), basename)

        self.templates = resource_filename('tiny', 'templates/')

    def setup_ebwt_idx(self):
        """Determines Bowtie index prefix and whether bowtie-build should run.
        self['ebwt'] is used for the bowtie commandline argument (see note below)
        self.paths['ebwt'] is the actual prefix path
        """

        # Empty values for ebwt (''/~/None) trigger bowtie-build
        self['run_bowtie_build'] = not bool(self.paths['ebwt'])

        if self['run_bowtie_build']:
            # Set the prefix to the run directory outputs. This is necessary
            # because workflow requires bt_index_files to be a populated list.
            self.paths['ebwt'] = self.get_ebwt_prefix()

        # verify_bowtie_build_outputs() will check if these end up being long indexes
        self['bt_index_files'] = self.get_bt_index_files()

        # When CWL copies bt_index_files for the bowtie.cwl InitialWorkDirRequirement, it does not
        # preserve the prefix path. What the workflow "sees" is the ebwt files at working dir root
        self['ebwt'] = os.path.basename(self.paths['ebwt'])

    def get_ebwt_prefix(self):
        """Determines the output prefix path for bowtie indexes that haven't been built yet. The basename
        of the prefix path is simply the basename of the reference genome sans file extension"""

        genome_files = self['reference_genome_files']
        if not genome_files:
            raise ValueError("If your Paths File doesn't have a value for \"ebtw:\", then bowtie indexes "
                             "will be built, but you'll need to provide your reference genome files under "
                             '"reference_genome_files:" (also in your Paths File)')

        genome_basename = os.path.basename(genome_files[0]['path'])
        return self.prefix(os.path.join( # prefix path:
            self['run_directory'], self['dir_name_bt_build'], genome_basename
        ))

    def get_bt_index_files(self):
        """Builds the list of expected bowtie index files from the ebwt prefix. If an index file
        doesn't exist then they will be automatically rebuilt from the user's reference genomes.
        File existence isn't checked if bowtie-build is already scheduled for this run."""

        try:
            verify_file_paths = not bool(self['run_bowtie_build'])
            prefix = self.paths['ebwt']
            ext = "ebwt"

            return [
                self.cwl_file(f"{prefix}.{subext}.{ext}", verify=verify_file_paths)
                for subext in ['1', '2', '3', '4', 'rev.1', 'rev.2']
            ]
        except FileNotFoundError as e:
            problem = "The following Bowtie index file couldn't be found:\n\t%s\n\n" % (e.filename,)
            rebuild = "Indexes will be built from your reference genome files during this run."
            userfix = "Please either correct your ebwt prefix or add reference genomes in the Paths File."

            if self['reference_genome_files']:
                print(problem + rebuild, file=sys.stderr)
                self.paths['ebwt'] = self.get_ebwt_prefix()
                self['run_bowtie_build'] = True
                return self.get_bt_index_files()
            else:
                sys.exit(problem + userfix)

    def validate_inputs(self):
        """For now, only GFF files are validated here"""

        gff_files = self.paths.get_gff_config()

        if gff_files:
            GFFValidator(
                gff_files,
                self.features_sheet.get_source_type_filters(),
                self.paths['ebwt'] if not self['run_bowtie_build'] else None,
                self.paths['reference_genome_files']
            ).validate()

    def execute_post_run_tasks(self, return_code):
        if self['run_bowtie_build']:
            self.verify_bowtie_build_outputs()

    def verify_bowtie_build_outputs(self):
        """Ensures that bowtie indexes were produced before saving the new ebwt prefix to the Paths File.
        If large indexes were produced, paths under bt_index_files need to be updated in the processed Run Config"""

        indexes = glob(os.path.join(self['run_directory'], self['dir_name_bt_build'], "*.ebwt*"))
        large_indexes = [f for f in indexes if f.endswith(".ebwtl")]

        # Update Paths File
        if indexes:
            self.paths.write_processed_config(self.paths.inf)

        # Update Run Config
        if large_indexes:
            for expected in self['bt_index_files']:
                expected['path'] += "l"
                assert expected['path'] in large_indexes
            self.write_processed_config()

    def save_run_profile(self, config_file_name=None) -> str:
        """Saves Samples Sheet and processed run config to the Run Directory for record keeping"""

        run_dir = self['run_directory']
        self.paths.save_run_profile(run_dir)
        self.samples_sheet.save_run_profile(run_dir)
        self.features_sheet.save_run_profile(run_dir)

        # The paths_* keys should now point to the copy produced above
        self['paths_file'] = self.cwl_file(os.path.join(run_dir, self.paths.basename))  # CWL file object
        self['paths_config'] = self.paths.basename                                      # User-facing value

        return self.write_processed_config(config_file_name)

    """========== COMMAND LINE =========="""

    @staticmethod
    def main():
        """Main routine to process the run information."""

        # Get input config file
        parser = argparse.ArgumentParser()
        required_group = parser.add_argument_group("required arguments")
        required_group.add_argument('-i', '--input-file', metavar='CONFIG', required=True,
                            help="The Run Config file to be processed")

        args = parser.parse_args()

        file_basename = os.path.basename(args.input_file)
        config_object = Configuration(args.input_file)
        config_object.write_processed_config(f"processed_{file_basename}")


def get_templates(context: str):
    """Copies a context-based subset of configuration file templates to the CWD

        Args:
            context: The command for which template files should be provided
                (currently: "tiny" or "tiny-count" or "tiny-plot")

        Returns: None
    """

    tiny_count = ['paths.yml', 'samples.csv', 'features.csv']
    tiny_plot =  ['tinyrna-light.mplstyle']
    tiny =       ['run_config_template.yml', *tiny_count, *tiny_plot]

    files_to_copy = {
        'tiny': tiny,
        'tiny-count': tiny_count,
        'tiny-plot': tiny_plot
    }.get(context)

    if files_to_copy is None:
        raise ValueError(f"Invalid template file context: {context}")
    else:
        conflicts = [f for f in files_to_copy if os.path.exists(f)]
        if conflicts:
            sys.exit(
                "The following files already exist in the current directory:\n\t"
                + '\n\t'.join(conflicts) + "\nPlease remove or rename them and try again."
            )

    # Copy template files to the current working directory
    templates_path = resource_filename('tiny', 'templates')
    for template in files_to_copy:
        shutil.copyfile(f"{templates_path}/{template}", f"{os.getcwd()}/{template}")


class PathsFile(ConfigBase):
    """A configuration class for managing and validating Paths Files.
    Relative paths are automatically resolved on lookup and list types are enforced.
    While this is convenient, developers should be aware of the following caveats:
        - Lookups that return list values do not return the original object; don't
          expect modifications to stick. If appending, use append_to().
        - Chained assignments can produce unexpected results.

    Args:
        file: The path to the Paths File, as a string
        in_pipeline: True only when utilized by a step in the workflow,
            in which case input files are sourced from the working directory
            regardless of the path indicated in the Paths File
    """

    # Parameter types
    required = ('samples_csv', 'features_csv')
    single =   ('samples_csv', 'features_csv', 'plot_style_sheet', 'adapter_fasta')
    groups =   ('reference_genome_files', 'gff_files')
    prefix =   ('ebwt', 'run_directory', 'tmp_directory')

    # Parameters that should be held constant between resume runs
    resume_forbidden = ('run_directory', 'ebwt', 'reference_genome_files')

    def __init__(self, file: str, in_pipeline=False):
        super().__init__(file)
        self.in_pipeline = in_pipeline
        self.map_path = self.from_pipeline if in_pipeline else self.from_here
        self.check_backward_compatibility()
        self.validate_paths()

    @staticmethod
    def from_pipeline(value):
        """When tiny-count runs as a pipeline step, all file inputs are
        sourced from the working directory regardless of original path."""

        if ConfigBase.is_path_dict(value, empty_ok=True):
            return dict(value, path=os.path.basename(value['path']))
        elif ConfigBase.is_path_str(value, empty_ok=True):
            return os.path.basename(value)
        else:
            return value

    def __getitem__(self, key: str):
        """Automatically performs path resolution for both single and group parameters.
        Note that only keys listed in self.groups are guaranteed to be returned as lists."""

        value = self.config.get(key)
        if key in self.groups:
            if value is None: return []
            return [self.map_path(p) for p in value]
        else:
            return self.map_path(value)

    def as_cwl_file_obj(self, key: str):
        """Returns the specified parameter with file paths converted to CWL file objects."""

        val = self[key]

        if not val:
            return val
        elif key in self.single:
            return self.cwl_file(val)
        elif key in self.groups:
            return [self.cwl_file(sub) for sub in val if self.is_path(sub)]
        elif key in self.prefix:
            raise ValueError(f"The parameter {key} isn't meant to be a CWL file object.")
        else:
            raise ValueError(f'Unrecognized parameter: "{key}"')

    def validate_paths(self):
        assert all(self[req] for req in self.required), \
            "The following parameters are required in {selfname}: {params}" \
            .format(selfname=self.basename, params=', '.join(self.required))

        # The availability of these file entries in the working directory will vary by step.
        # This is determined by the step's CWL CommandLineTool specification.
        # Instead of checking file existence within each step that uses this class,
        # check only at pipeline startup and let the workflow runner worry about files from there.
        if self.in_pipeline: return

        for key in self.single:
            resolved_path = self[key]
            if not resolved_path: continue
            assert os.path.isfile(resolved_path), \
                "The file provided for {key} in {selfname} could not be found:\n\t{file}" \
                .format(key=key, selfname=self.basename, file=resolved_path)

        for key in self.groups:
            for entry in self[key]:
                if isinstance(entry, dict): entry = entry.get('path')
                if not entry: continue
                assert os.path.isfile(entry), \
                    "The following file provided under {key} in {selfname} could not be found:\n\t{file}" \
                    .format(key=key, selfname=self.basename, file=entry)

    def check_backward_compatibility(self):
        assert 'gff_files' in self, \
            "The gff_files parameter was not found in your Paths File. This likely means " \
            "that you are using a Paths File from an earlier version of tinyRNA. Please " \
            "check the release notes and update your configuration files."

        missing_keys = [key for key in (*self.single, *self.groups, *self.prefix)
                        if key not in self.config]

        assert not missing_keys, \
            "The following expected keys were missing in {selfname}:\n\t{missing}" \
            .format(selfname=self.basename, missing="\n\t".join(missing_keys))

    def get_gff_config(self) -> Dict[str, list]:
        """Restructures GFF input info so that it can be more easily handled.
        To be clear, the Paths File YAML could be structured to match the desired output,
        but the current format was chosen because it's more readable with more forgiving syntax.

        The YAML format is [{"path": "gff_path_1", "alias": [alias1, alias2, ...]}, { ... }, ...]
        The output format is {"gff_path_1": [alias1, alias2, ...], ...}
        """

        gff_files = defaultdict(list)
        id_filter = lambda alias: alias.lower() != 'id'

        # Build dictionary of files and allowed aliases
        for gff in self['gff_files']:
            if not self.is_path_dict(gff): continue
            alias = gff.get('alias')
            path = gff['path']

            # Allow for some user error in YAML syntax
            if isinstance(alias, str): alias = [alias]
            if not isinstance(alias, list): alias = []
            gff_files[path].extend(filter(id_filter, alias))

        # Remove duplicate aliases per file, keep order
        for file, alias in gff_files.items():
            gff_files[file] = sorted(set(alias), key=alias.index)

        if not len(gff_files) and not self.in_pipeline:
            self.print_sequence_counting_notice()

        return dict(gff_files)

    def print_sequence_counting_notice(self):
        print("No GFF files were specified in {selfname}.\n"
              "Reads will be quantified by sequence "
              "rather than by feature."
              .format(selfname=self.basename))

        for s in reversed(range(6)):
            print(f"Proceeding in {s}s", end='\r')
            time.sleep(1)

    # Override
    def append_to(self, key: str, val: Any):
        """Overrides method in the base class. This is necessary due to automatic
        path resolution in __getitem__ which returns a *temporary* list value, and
        items appended to the temporary list would otherwise be lost."""

        assert key in self.groups, "Tried appending to a non-list type parameter"

        target = self.config.get(key)
        if not isinstance(target, list):
            self.config[key] = target = []

        target.append(val)
        return target

    def save_run_profile(self, run_directory):
        """Saves a copy of the Paths File to the Run Directory with amended paths.
        Note the distinction between out_obj[key] and self[key]. The latter performs
        automatic path resolution, whereas out_obj is essentially just a dict."""

        out_obj = copy.deepcopy(self.config)
        out_file = os.path.join(run_directory, self.basename)

        adjacent_paths = self.required
        absolute_paths = [path for path in (*self.single, *self.prefix)
                          if path not in ("run_directory", *self.required)]

        for adjacent in adjacent_paths:
            out_obj[adjacent] = os.path.basename(self[adjacent])

        for key in absolute_paths:
            if not self.is_path_str(self[key]): continue
            out_obj[key] = os.path.abspath(self[key])

        for key in self.groups:
            for i, entry in enumerate(self[key]):
                if self.is_path_dict(entry):
                    out_obj[key][i]['path'] = os.path.abspath(entry['path'])
                elif self.is_path_str(entry):
                    out_obj[key][i] = os.path.abspath(entry)

        with open(out_file, 'w') as f:
            self.yaml.dump(out_obj, f)


class SamplesSheet:
    def __init__(self, file, context):
        self.csv = CSVReader(file, "Samples Sheet")
        self.basename = os.path.basename(file)
        self.dir = os.path.dirname(file)
        self.file = file

        self.hts_samples = []
        self.groups_reps = []
        self.normalizations = []
        self.sample_basenames = []
        self.control_condition = None
        self.is_compatible_df = False

        map_path, validate_path = self.get_context_methods(context).values()
        self.validate_sample_path = validate_path
        self.map_path = map_path
        self.context = context

        self.read_csv()

    def get_context_methods(self, context):
        return {
            "Standalone Run": {'map_path': self.from_here, 'validate_path': self.validate_alignments_filepath},
            "Pipeline Start": {'map_path': self.from_here, 'validate_path': self.validate_fastq_filepath},
            "Pipeline Step":  {'map_path': self.from_pipeline, 'validate_path': lambda _: True},  # skip validation
        }[context]

    def read_csv(self):
        reps_per_group = Counter()
        try:
            for row in self.csv.rows():
                hts_sample = self.map_path(row['File'])
                group_name = row['Group']
                rep_number = row['Replicate']
                is_control = row['Control'].lower() == 'true'
                norm_prefs = row['Normalization']
                basename = self.get_sample_basename(hts_sample)

                self.validate_sample_path(hts_sample)
                self.validate_group_rep(group_name, rep_number)
                self.validate_control_group(is_control, group_name)
                self.validate_normalization(norm_prefs)

                self.hts_samples.append(hts_sample)
                self.sample_basenames.append(basename)
                self.groups_reps.append((group_name, rep_number))
                self.normalizations.append(norm_prefs)
                reps_per_group[group_name] += 1

                if is_control: self.control_condition = group_name
        except Exception as e:
            if not isinstance(e, AssertionError):
                # Validation steps already indicate row, exception likely from something else
                msg = f"Error occurred on line {self.csv.row_num} of {self.basename}"
                append_to_exception(e, msg)
            raise

        self.is_compatible_df = self.validate_deseq_compatibility(reps_per_group)

    def from_here(self, input_file):
        """Resolves .sam/.bam paths in pipeline startup and standalone context"""

        return ConfigBase.joinpath(self.dir, input_file)

    def from_pipeline(self, input_file):
        """Resolves .fastq(.gz) file paths in pipeline step context"""

        basename_root, ext = os.path.splitext(os.path.basename(input_file))
        return f"{basename_root}_aligned_seqs.sam"

    def validate_alignments_filepath(self, file: str):
        """Checks file existence, extension, and duplicate entries in standalone context"""

        root, ext = os.path.splitext(file)

        assert os.path.isfile(file), \
            "The file on row {row_num} of {selfname} was not found:\n\t{file}" \
            .format(row_num=self.csv.row_num, selfname=self.basename, file=file)

        assert ext in (".sam", ".bam"), \
            "Files in {selfname} must have a .sam or .bam extension (row {row_num})" \
            .format(selfname=self.basename, row_num=self.csv.row_num)

        assert file not in self.hts_samples, \
            "Alignment files cannot be listed more than once in {selfname} (row {row_num})" \
            .format(selfname=self.basename, row_num=self.csv.row_num)

    def validate_fastq_filepath(self, file: str):
        """Checks file existence, extension, and duplicate entries in pipeline startup context"""

        root, ext = os.path.splitext(file)

        assert os.path.isfile(file), \
            "The fastq file on row {row_num} of {selfname} was not found:\n\t{file}" \
            .format(row_num=self.csv.row_num, selfname=self.basename, file=file)

        assert ext in (".fastq", ".gz"), \
            "Files in {selfname} must have a .fastq(.gz) extension (row {row_num})" \
            .format(selfname=self.basename, row_num=self.csv.row_num)

        assert file not in self.hts_samples, \
            "Fastq files cannot be listed more than once in {selfname} (row {row_num})" \
            .format(selfname=self.basename, row_num=self.csv.row_num)

    def validate_group_rep(self, group:str, rep:str):
        assert (group, rep) not in self.groups_reps, \
            "The same group and replicate number cannot appear on " \
            "more than one row in {selfname} (row {row_num})" \
            .format(selfname=self.basename, row_num=self.csv.row_num)

    def validate_control_group(self, is_control: bool, group: str):

        if not is_control or self.context != "Pipeline Start":
            return

        assert self.control_condition in (group, None), \
            "Experiments with multiple control conditions aren't supported " \
            "(row {row_num} in {selfname}).\nHowever, if the control condition " \
            "is unspecified, all possible comparisons will be made and this " \
            "should accomplish your goal." \
            .format(row_num=self.csv.row_num, selfname=self.basename)

    def validate_normalization(self, norm):
        assert re.fullmatch(r"\s*((?:\d+(?:\.\d*)?|\.\d+)|(rpm))\s*", norm, re.IGNORECASE) or not norm, \
            "Invalid normalization value in {selfname} (row {row_num})" \
            .format(selfname=self.basename, row_num=self.csv.row_num)

    def validate_deseq_compatibility(self, sample_groups: Counter) -> Optional[bool]:

        if self.context != "Pipeline Start":
            return None

        SamplesSheet.validate_r_safe_sample_groups(sample_groups)

        total_samples = sum(sample_groups.values())
        total_coefficients = len(sample_groups)
        degrees_of_freedom = total_samples - total_coefficients

        if degrees_of_freedom < 1:
            print("Your experiment design has less than one degree of freedom, which is incompatible "
                  "with DESeq2. The tiny-deseq step will be skipped and DGE plots will not be produced.",
                  file=sys.stderr)
            return False
        else:
            return True

    @staticmethod
    def validate_r_safe_sample_groups(sample_groups: Counter):
        """Determine the "syntactically valid" translation of each group name to ensure
        that two groups won't share the same name once translated in tiny-deseq.r"""

        safe_names = defaultdict(list)
        for group in sample_groups:
            safe_names[get_r_safename(group)].append(group)

        collisions = [' ≈ '.join(cluster) for cluster in safe_names.values() if len(cluster) > 1]

        assert len(collisions) == 0, \
            "The following group names are too similar and will cause a namespace collision in R:\n" \
            + '\n'.join(collisions)

    def save_run_profile(self, run_directory):
        """Writes a copy of the CSV with absolute paths"""

        outfile = os.path.join(run_directory, self.basename)
        header = CSVReader.tinyrna_sheet_fields['Samples Sheet'].keys()
        coldata = zip(self.hts_samples, self.groups_reps, self.normalizations)

        with open(outfile, 'w', newline='') as out_csv:
            csv_writer = csv.writer(out_csv)
            csv_writer.writerow(header)
            for sample, (group, rep), norm in coldata:
                control = (group == self.control_condition) or ""
                sample = os.path.abspath(sample)
                csv_writer.writerow([sample, group, rep, control, norm])

    @staticmethod
    def get_sample_basename(filename):
        root, _ = os.path.splitext(filename)
        return os.path.basename(root)


class FeaturesSheet:
    def __init__(self, file, context):
        self.csv = CSVReader(file, "Features Sheet")
        self.basename = os.path.basename(file)
        self.dir = os.path.dirname(file)
        self.context = context
        self.file = file

        self.rules = []
        self.read_csv()

    def read_csv(self):
        try:
            rules, hierarchies = [], []
            for rule in self.csv.rows():
                rule['nt5end'] = rule['nt5end'].upper().translate({ord('U'): 'T'})  # Convert RNA base to cDNA base
                rule['Identity'] = (rule.pop('Key'), rule.pop('Value'))             # Create identity tuple
                rule['Overlap'] = rule['Overlap'].lower()                           # Built later in reference parsers
                hierarchy = int(rule.pop('Hierarchy'))                              # Convert hierarchy to number

                # Duplicate rules are screened out here
                # Equality check omits hierarchy value
                if rule not in rules:
                    rules.append(rule)
                    hierarchies.append(hierarchy)
        except Exception as e:
            msg = f"Error occurred on line {self.csv.row_num} of {self.basename}"
            append_to_exception(e, msg)
            raise

        # Reunite hierarchy values with their rules
        self.rules = [
            dict(rule, Hierarchy=hierarchy)
            for rule, hierarchy in zip(rules, hierarchies)
        ]

    def get_source_type_filters(self):
        """Returns only the Source Filter and Type Filter columns"""

        interests = ("Filter_s", "Filter_t")
        return [{selector: rule[selector] for selector in interests}
                for rule in self.rules]

    def save_run_profile(self, run_directory):
        """Copies the Features Sheet to the run directory"""

        outfile = os.path.join(run_directory, self.basename)
        shutil.copyfile(self.file, outfile)


class CSVReader(csv.DictReader):
    """A simple wrapper class for csv.DictReader

    This makes field labels consistent across the project, simplifies the code, and
    allows for validation and reordering of column names. We also keep track of the
    row number for diagnostic outputs; the base class offers the line_num attribute,
    but line_num != row_num if a record spans multiple lines in the csv.
    """

    # user-facing name -> internal short name
    tinyrna_sheet_fields = {
        "Features Sheet": OrderedDict({
           "Select for...":     "Key",
           "with value...":     "Value",
           "Classify as...":    "Class",
           "Source Filter":     "Filter_s",
           "Type Filter":       "Filter_t",
           "Hierarchy":         "Hierarchy",
           "Strand":            "Strand",
           "5' End Nucleotide": "nt5end",
           "Length":            "Length",
           "Overlap":           "Overlap",
           "Mismatches":        "Mismatch"
        }),
        "Samples Sheet": OrderedDict({
            "Input Files":       "File",
            "Sample/Group Name": "Group",
            "Replicate Number":  "Replicate",
            "Control":           "Control",
            "Normalization":     "Normalization"
        })
    }

    def __init__(self, filename: str, doctype: str = None):
        self.doctype = doctype
        self.tinyrna_file = filename
        self.row_num = 0
        try:
            self.tinyrna_fields = tuple(CSVReader.tinyrna_sheet_fields[doctype].values())
        except KeyError as ke:
            raise ValueError("w-HH-at 'n t'heck are you doin")

    def rows(self):
        self.replace_excel_ellipses()
        with open(os.path.expanduser(self.tinyrna_file), 'r', encoding='utf-8-sig', newline='') as f:
            super().__init__(f, fieldnames=self.tinyrna_fields, dialect=get_csv_dialect(f))
            header = next(self)

            # Compatibility check. Column headers are still often changed at this stage
            # and it doesn't make sense to infer column identity
            self.validate_csv_header(header)

            # For skipping rows that are empty or all whitespace cells
            non_ws_cell = lambda cell: re.sub(r'^\s+$', '', cell)
            empty_row = lambda row: not any(non_ws_cell(x) for x in row.values())

            for row in self:
                self.row_num += 1
                if empty_row(row): continue
                yield row

    def replace_excel_ellipses(self):
        """Excel has an autocorrect setting that converts three periods to an ellipsis character.
        The resulting ellipsis is not UTF-8 compatible and causes the decoder to fail. Switching
        to the correct encoding is still problematic because the character, while visually correct,
        doesn't compare equal to "..." so header validation fails. Just replace it."""

        try:
            with open(self.tinyrna_file, 'r', encoding='utf-8') as f:
                next(f)
        except UnicodeDecodeError as e:
            if e.reason == "invalid continuation byte" and e.object[e.start] == 0xc9:
                with open(self.tinyrna_file, 'r', encoding='latin-1') as f:
                    lines = f.readlines()
                    lines[0] = lines[0].replace("\xc9", '...')
                    lines[0].encode('latin-1').decode('utf-8-sig')  # sanity check
                # Encode with utf-8 rather than utf-8-sig for Excel compatibility
                with open(self.tinyrna_file, 'w', encoding='utf-8') as f:
                    f.writelines(lines)
                print(f"Replaced invalid character in {self.doctype} header.")
            else:
                raise

    def validate_csv_header(self, header: OrderedDict):
        # The expected header values
        doc_reference = self.tinyrna_sheet_fields[self.doctype]
        expected = {key.lower() for key in doc_reference.keys()}

        # The header values that were read
        read_vals = {val.lower() for key, val in header.items() if None not in (key, val)}
        read_vals.update(val.lower() for val in header.get(None, ()))  # Extra headers
        self.check_backward_compatibility(read_vals)

        # Find differences between actual and expected headers
        unknown = {col_name for col_name in read_vals if col_name not in expected}
        missing = expected - read_vals

        if len(missing):
            raise ValueError('\n\t'.join([f"The following columns are missing from your {self.doctype}:", *missing]))
        if len(unknown):
            raise ValueError('\n\t'.join([f"The following columns in your {self.doctype} are unrecognized:", *unknown]))

        doc_ref_lowercase = {key.lower(): value for key, value in doc_reference.items()}
        header_lowercase = {key: value.lower() for key, value in header.items()}

        if tuple(header_lowercase.values()) != tuple(doc_ref_lowercase.keys()):
            # Remap column order to match client's
            self.fieldnames = tuple(doc_ref_lowercase[key] for key in header_lowercase.values())

    def check_backward_compatibility(self, header_vals):
        """Catch column differences from old versions so a helpful error can be provided"""

        header_vals_lc = {val.lower() for val in header_vals}

        compat_errors = []
        if self.doctype == "Features Sheet":
            if len(header_vals_lc & {'alias by...', 'feature source'}):
                compat_errors.append('\n'.join([
                    "It looks like you're using a Features Sheet from an earlier version of",
                    "tinyRNA. Feature aliases and GFF files are now defined in the Paths File.",
                    "Please review the Paths File documentation in Configuration.md, update your",
                    'Paths File, and remove the "Alias by..." and "Feature Source" columns from',
                    "your Features Sheet to avoid this error."
                ]))

            if len(header_vals_lc & {'source filter', 'type filter'}) != 2:
                compat_errors.append('\n'.join([
                    "It looks like you're using a Features Sheet from an earlier version of",
                    "tinyRNA. Source and type filters are now defined in the Features Sheet.",
                    "They are no longer defined in the Run Config. Please review the Stage 1",
                    "section in tiny-count's documentation, then add the new columns",
                    '"Source Filter" and "Type Filter" to your Features Sheet to avoid this error.'
                ]))

            if 'tag' in header_vals_lc:
                compat_errors.append('\n'.join([
                    "It looks like you're using a Features Sheet from a version of tinyRNA",
                    'that offered "tagged counting". The "Tag" header has been repurposed as a feature',
                    "classifier and its meaning within the pipeline has changed. Additionally, feature",
                    "class is no longer determined by the Class= attribute. Please review the Stage 1",
                    'section in tiny-count\'s documentation, then rename the "Tag" column to',
                    '"Classify as..." to avoid this error.'
                ]))

            if 'mismatches' not in header_vals_lc:
                compat_errors.append('\n'.join([
                    "It looks like you're using a Features Sheet from an earlier version of",
                    'tinyRNA. An additional column, "Mismatches", is now expected. Please review',
                    "the Stage 2 section in tiny-count's documentation for more info, then add",
                    "the new column to your Features Sheet to avoid this error."
                ]))

        if self.doctype == "Samples Sheet":
            if 'fastq/sam files' in header_vals_lc:
                compat_errors.append('\n'.join([
                    "It looks like you're using a Samples Sheet from an earlier version of",
                    'tinyRNA. The "FASTQ/SAM files" column has been renamed to "Input Files"',
                    'due to the addition of BAM file support. Please rename the column in',
                    "your Samples Sheet to avoid this error."
                ]))

        if compat_errors: raise ValueError('\n\n'.join(compat_errors))


if __name__ == '__main__':
    Configuration.main()