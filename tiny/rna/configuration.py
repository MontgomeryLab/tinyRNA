import ruamel.yaml
import argparse
import shutil
import errno
import sys
import csv
import os
import re

from pkg_resources import resource_filename
from collections import Counter, OrderedDict, defaultdict
from typing import Union, Any, Optional, List
from glob import glob

from tiny.rna.counter.validation import GFFValidator
from tiny.rna.util import get_timestamp, get_r_safename


class ConfigBase:
    """Base class for basic tinyrna configuration operations

    Attributes:
        yaml: the YAML interface for reading config and writing processed config
        config: the configuration object produced by loading the config file
        inf: the filename of the configuration .yml file to process
        dir: parent directory of the input file. Used for calculating paths relative to config file.
        templates: path to the package templates directory
        dt: a date-time string for default output naming
    """

    def __init__(self, config_file: str):
        if '~' in config_file: config_file = os.path.expanduser(config_file)
        self.dir = os.path.dirname(os.path.abspath(config_file))
        self.inf = config_file
        self.basename = os.path.basename(config_file)
        self.templates = ''
        self.dt = ''

        self.yaml = ruamel.yaml.YAML()
        with open(config_file, 'r') as f:
            self.config = self.yaml.load(f)

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
        """Apply the setting if it has not been previously set"""
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

        if (
                isinstance(destination, dict) and
                isinstance(destination.get("path"), (str, bytes)) and
                len(destination['path'])
        ):
            return dict(destination, path=self.joinpath(origin, destination["path"]))
        elif isinstance(destination, (str, bytes)) and bool(destination):
            return self.joinpath(origin, destination)
        else:
            return destination

    def setup_step_inputs(self):
        """For now, only tiny-plot requires additional setup for step inputs
        This function is called at both startup and resume"""

        def setup_tiny_plot_inputs():
            cs_filter = 'plot_class_scatter_filter'
            style_req = ['include', 'exclude']
            classes = self.get(cs_filter, {}).get('classes')  # backward compatibility
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
            os.mkdir(run_dir)

        return run_dir

    def get_outfile_path(self, infile: str = None) -> str:
        """Returns the path and file for the processed run config"""
        if infile is None: infile = self.inf
        return self.joinpath(self['run_directory'], os.path.basename(infile))

    def write_processed_config(self, filename: str = None) -> str:
        """Writes the current configuration to disk"""
        if filename is None: filename = self.get_outfile_path(self.inf)

        with open(filename, 'w') as outconf:
            if 'paths_config' in self and not os.path.isabs(self['paths_config']):
                # Processed config will be written to the Run Directory
                # Ensure paths_config is an absolute path so it remains valid
                self['paths_config'] = self.from_here(self['paths_config'])

            self.yaml.dump(self.config, outconf)

        return filename


class Configuration(ConfigBase):
    """A class for processing and updating a YAML config file for CWL

    Ultimately, this class populates workflow settings and per-library settings. This
    is a convenience to the user as it is tedious to define inputs and outputs pertaining
    to each workflow step. Settings are determined by the Paths, Samples, and Features Sheets.
    Users may provide both relative and absolute paths

    IMPORTANT: Paths provided in any config file are evaluated relative to the containing config file.

    Attributes:
        paths: the configuration object from processing the paths_config file.
            This holds path info for other config files and prefixes, and is updated
            appropriately if 'run_bowtie_index' is set to 'true'
    """

    def __init__(self, config_file: str, validate_inputs=False):
        # Parse YAML configuration file
        super().__init__(config_file)

        self.paths = self.load_paths_config()
        self.assimilate_paths_file()

        self.setup_pipeline()
        self.setup_file_groups()
        self.setup_ebwt_idx()
        self.process_samples_sheet()
        self.process_features_sheet()
        self.setup_step_inputs()
        if validate_inputs: self.validate_inputs()

    def load_paths_config(self):
        """Constructs a sub-configuration object containing workflow file preferences
        self['paths_config'] is the user-facing file path (just the path string)
        self['paths_file'] is a CWL file object used as a workflow input."""
        paths_file_path = self.from_here(self['paths_config'])
        return PathsFile(paths_file_path)

    def assimilate_paths_file(self):
        for key in [*PathsFile.single, *PathsFile.groups]:
            self[key] = self.paths.as_cwl_file_obj(key)
        for key in PathsFile.prefix:
            self[key] = self.paths[key]
        self['paths_file'] = self.cwl_file(self.paths.inf)

    def process_samples_sheet(self):
        samples_sheet_path = self.paths['samples_csv']
        samples_sheet = SamplesSheet(samples_sheet_path)

        self['sample_basenames'] = samples_sheet.sample_basenames
        self['control_condition'] = samples_sheet.control_condition
        self['run_deseq'] = samples_sheet.is_compatible_df

        self['in_fq'] = [self.cwl_file(fq, verify=False) for fq in samples_sheet.fastq_files]
        self['fastp_report_titles'] = [f"{g}_rep_{r}" for g, r in samples_sheet.groups_reps]

    def process_features_sheet(self) -> List[dict]:
        """Retrieves GFF Source and Type Filter definitions for use in GFFValidator"""
        features_sheet_path = self.paths['features_csv']
        reader = CSVReader(features_sheet_path, "Features Sheet").rows()

        interests = ("Filter_s", "Filter_t")
        return [{selector: rule[selector] for selector in interests}
                for rule in reader]

    def setup_file_groups(self):
        """Configuration keys that represent lists of files"""

        self.set_default_dict({per_file_setting_key: [] for per_file_setting_key in
            ['in_fq', 'sample_basenames', 'gff_files', 'fastp_report_titles']
        })
            
    def setup_pipeline(self):
        """Overall settings for the whole pipeline"""

        self.dt = get_timestamp()
        self['run_date'], self['run_time'] = self.dt.split('_')

        default_run_name = '_'.join(x for x in [self['user'], "tinyrna"] if x)
        self['run_name'] = self.get('run_name', default=default_run_name) + "_" + self.dt

        # Create prefixed Run Directory name
        run_dir_parent, run_dir = os.path.split(self['run_directory'].rstrip(os.sep))
        self['run_directory'] = self.joinpath(run_dir_parent, self['run_name'] + "_" + run_dir)

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

        selection_rules = self.process_features_sheet()
        gff_files = {gff['path']: [] for gff in self['gff_files']}
        ebwt = self.paths['ebwt'] if not self['run_bowtie_build'] else None

        GFFValidator(
            gff_files,
            selection_rules,
            ebwt=ebwt,
            genomes=self.paths['reference_genome_files'],
            alignments=None  # Used in tiny-count standalone runs
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

        from importlib.metadata import version
        self['version'] = version('tinyrna')

        samples_sheet_name = os.path.basename(self['samples_csv']['path'])
        shutil.copyfile(self['samples_csv']['path'], f"{self['run_directory']}/{samples_sheet_name}")
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
    }.get(context, None)

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
          append to them. Instead, use the append_to() helper function.
        - Chained assignments can produce unexpected results.
    """

    # Parameter types
    required = ('samples_csv', 'features_csv', 'gff_files')
    single =   ('samples_csv', 'features_csv', 'plot_style_sheet', 'adapter_fasta')
    groups =   ('reference_genome_files', 'gff_files')
    prefix =   ('ebwt', 'run_directory', 'tmp_directory')

    # Parameters that need to be held constant between resume runs for analysis integrity
    resume_forbidden = ('samples_csv', 'run_directory', 'ebwt', 'reference_genome_files')

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

        if isinstance(value, dict) and value.get("path") is not None:
            return dict(value, path=os.path.basename(value['path']))
        elif isinstance(value, (str, bytes)):
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
            return [self.cwl_file(sub) for sub in val if sub]
        elif key in self.prefix:
            raise ValueError(f"The parameter {key} isn't meant to be a CWL file object.")
        else:
            raise ValueError(f'Unrecognized parameter: "{key}"')

    def validate_paths(self):
        assert all(self[req] for req in self.required), \
            "The following parameters are required in {selfname}: {params}" \
            .format(selfname=self.basename, params=', '.join(self.required))

        assert any(gff.get('path') for gff in self['gff_files']), \
            "At least one GFF file path must be specified under gff_files in {selfname}" \
            .format(selfname=self.basename)

        # Some entries in Paths File are omitted from tiny-count's working directory during
        #  pipeline runs. There is no advantage to checking file existence here vs. in load_*
        if self.in_pipeline: return

        for key in self.single:
            resolved_path = self[key]
            if not resolved_path: continue
            assert os.path.isfile(resolved_path), \
                "The file provided for {key} in {selfname} could not be found:\n\t{file}" \
                .format(key=key, selfname=self.basename, file=resolved_path)

        for key in self.groups:
            for entry in self[key]:
                if isinstance(entry, dict): entry = entry['path']
                assert os.path.isfile(entry), \
                    "The following file provided under {key} in {selfname} could not be found:\n\t{file}" \
                    .format(key=key, selfname=self.basename, file=entry)

    def check_backward_compatibility(self):
        assert 'gff_files' in self, \
            "The gff_files parameter was not found in your Paths File. This likely means " \
            "that you are using a Paths File from an earlier version of tinyRNA. Please " \
            "check the release notes and update your configuration files."

    # Override
    def append_to(self, key: str, val: Any):
        """Overrides method in the base class. This is necessary due to automatic
        path resolution in __getitem__ which returns a *temporary* list value, and
        items appended to the temporary list would otherwise be lost."""

        assert key in self.groups, "Tried appending to a non-list type parameter"
        target = self.config.get(key, [])
        target.append(val)
        return target


class SamplesSheet:
    def __init__(self, file):
        self.csv = CSVReader(file, "Samples Sheet")
        self.basename = os.path.basename(file)
        self.dir = os.path.dirname(file)
        self.file = file

        self.fastq_files = []
        self.groups_reps = []
        self.sample_basenames = []
        self.control_condition = None
        self.is_compatible_df = False

        self.read_csv()

    def read_csv(self):
        reps_per_group = Counter()
        for row in self.csv.rows():
            fastq_file = Configuration.joinpath(self.dir, row['File'])
            group_name = row['Group']
            rep_number = row['Replicate']
            is_control = row['Control'].lower() == 'true'
            basename = self.get_sample_basename(fastq_file)

            self.validate_fastq_filepath(fastq_file)
            self.validate_group_rep(group_name, rep_number)
            self.validate_control_group(is_control, group_name)

            self.fastq_files.append(fastq_file)
            self.sample_basenames.append(basename)
            self.groups_reps.append((group_name, rep_number))
            reps_per_group[group_name] += 1

            if is_control: self.control_condition = group_name

        self.is_compatible_df = self.validate_deseq_compatibility(reps_per_group)

    def validate_fastq_filepath(self, file: str):
        """Checks file existence, extension, and duplicate entries.
        Args:
            file: fastq file path. For  which has already been resolved relative to self.dir
        """

        root, ext = os.path.splitext(file)

        assert os.path.isfile(file), \
            "The fastq file on row {row_num} of {selfname} was not found:\n\t{file}" \
            .format(row_num=self.csv.row_num, selfname=self.basename, file=file)

        assert ext in (".fastq", ".gz"), \
            "Files in {selfname} must have a .fastq(.gz) extension (row {row_num})"\
            .format(selfname=self.basename, row_num=self.csv.row_num)

        assert file not in self.fastq_files, \
            "Fastq files cannot be listed more than once in {selfname} (row {row_num})"\
            .format(selfname=self.basename, row_num=self.csv.row_num)

    def validate_group_rep(self, group:str, rep:str):
        assert (group, rep) not in self.groups_reps, \
            "The same group and replicate number cannot appear on " \
            "more than one row in {selfname} (row {row_num})"\
            .format(selfname=self.basename, row_num=self.csv.row_num)

    def validate_control_group(self, is_control: bool, group: str):
        if not is_control: return
        assert self.control_condition in (group, None), \
            "tinyRNA does not support multiple control conditions " \
            "(row {row_num} in {selfname}).\nHowever, if the control condition " \
            "is unspecified, all possible comparisons will be made and this " \
            "should accomplish your goal."\
            .format(row_num=self.csv.row_num, selfname=self.basename)

    @staticmethod
    def validate_deseq_compatibility(sample_groups: Counter) -> bool:
        SamplesSheet.validate_r_safe_sample_groups(sample_groups)

        total_samples = sum(sample_groups.values())
        total_coefficients = len(sample_groups)
        degrees_of_freedom = total_samples - total_coefficients

        if degrees_of_freedom < 1:
            print("Your experiment design has less than one degree of freedom, which is incompatible "
                  "with DESeq2. The DGE step will be skipped and most plots will not be produced.",
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

    @staticmethod
    def get_sample_basename(filename):
        root, _ = os.path.splitext(filename)
        return os.path.basename(root)


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
        }),
        "Samples Sheet": OrderedDict({
            "FASTQ/SAM Files":   "File",
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
        with open(os.path.expanduser(self.tinyrna_file), 'r', encoding='utf-8-sig') as f:
            super().__init__(f, fieldnames=self.tinyrna_fields, delimiter=',')
            header = next(self)

            # Compatibility check. Column headers are still often changed at this stage
            # and it doesn't make sense to infer column identity
            self.validate_csv_header(header)

            for row in self:
                self.row_num += 1
                yield row

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
        compat_errors = []
        if self.doctype == "Features Sheet":
            if len(header_vals & {'alias by...', 'feature source'}):
                compat_errors.append('\n'.join([
                    "It looks like you're using a Features Sheet from an earlier version of",
                    "tinyRNA. Feature aliases and GFF files are now defined in the Paths File.",
                    "Please review the Paths File documentation in Configuration.md, update your",
                    'Paths File, and remove the "Alias by..." and "Feature Source" columns from',
                    "your Features Sheet to avoid this error."
                ]))

            if len(header_vals & {'source filter', 'type filter'}) != 2:
                compat_errors.append('\n'.join([
                    "It looks like you're using a Features Sheet from an earlier version of",
                    "tinyRNA. Source and type filters are now defined in the Features Sheet.",
                    "They are no longer defined in the Run Config. Please review the Stage 1",
                    "section in tiny-count's documentation, then add the new columns",
                    '"Source Filter" and "Type Filter" to your Features Sheet to avoid this error.'
                ]))

            if 'tag' in header_vals:
                compat_errors.append('\n'.join([
                    "It looks like you're using a Features Sheet from a version of tinyRNA",
                    'that offered "tagged counting". The "Tag" header has been repurposed as a feature',
                    "classifier and its meaning within the pipeline has changed. Additionally, feature",
                    "class is no longer determined by the Class= attribute. Please review the Stage 1",
                    'section in tiny-count\'s documentation, then rename the "Tag" column to',
                    '"Classify as..." to avoid this error.'
                ]))

        if compat_errors: raise ValueError('\n\n'.join(compat_errors))


if __name__ == '__main__':
    Configuration.main()