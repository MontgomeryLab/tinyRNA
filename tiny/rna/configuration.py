import ruamel.yaml
import argparse
import shutil
import errno
import sys
import csv
import os
import re

from pkg_resources import resource_filename
from collections import Counter, OrderedDict
from datetime import datetime
from typing import Union, Any

from tiny.rna.counter.validation import GFFValidator

timestamp_format = re.compile(r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}")


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

    def __setitem__(self, key: str, val: Union[str, list, dict, bool]) -> Union[str, list, dict, bool]:
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
    def joinpath(path1: str, path2: str) -> str:
        """Combines two relative paths intelligently"""
        path1, path2 = (os.path.expanduser(p) for p in [path1, path2])
        if os.path.isabs(path2): return path2
        return os.path.normpath(os.path.join(path1, path2))

    @staticmethod
    def cwl_file(file: str, verify=True) -> dict:
        """Returns a minimal File object as defined by CWL"""
        if '~' in file: file = os.path.expanduser(file)
        if verify and not os.path.exists(file):
            raise(FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file))
        return {'class': 'File', 'path': file}

    @staticmethod
    def cwl_dir(dir: str) -> dict:
        """Returns a minimal Directory object as defined by CWL"""
        return {'class': 'Directory', 'path': dir}

    def from_here(self, destination: str, origin: str = None):
        """Calculates paths relative to the input config file"""
        origin = self.dir if origin is None else origin
        return self.joinpath(origin, destination)

    def create_run_directory(self) -> str:
        """Create the destination directory for pipeline outputs"""
        run_dir = self["run_directory"]
        if not os.path.isdir(run_dir):
            os.mkdir(run_dir)

        return run_dir

    def get_outfile_path(self, infile: str = None) -> str:
        """Returns the path and filename for the processed run config"""
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
        self.process_paths_sheet()
        
        self.setup_pipeline()
        self.setup_file_groups()
        self.setup_ebwt_idx()
        self.process_sample_sheet()
        self.process_feature_sheet()
        if validate_inputs: self.validate_inputs()
        
    def load_paths_config(self):
        """Constructs a sub-configuration object containing workflow file preferences"""
        path_sheet = self.from_here(self['paths_config'])
        return ConfigBase(path_sheet)

    def process_paths_sheet(self):
        """Loads file paths and groups from the Paths File and resolves relative paths"""

        def to_cwl_file_class(input_file_path):
            path_to_input = self.paths.from_here(input_file_path)
            return self.cwl_file(path_to_input)

        required = ('samples_csv', 'features_csv', 'gff_files')
        single =   ('samples_csv', 'features_csv', 'plot_style_sheet', 'adapter_fasta')
        groups =   ('reference_genome_files', 'gff_files')

        # Convert file path strings to CWL File objects where necessary
        try:
            # Individual path strings
            for key in single:
                if not self.paths[key] and key not in required: continue
                self[key] = to_cwl_file_class(self.paths[key])
    
            # Lists of path strings
            for key in groups:
                self[key] = [
                    to_cwl_file_class(file) for file in self.paths[key]
                    if key in required or bool(file)
                ]
        except FileNotFoundError as e:
            msg = f"The file provided for {key} in your Paths File could not be found"
            msg += "\n\t" + e.filename
            sys.exit(msg)

        # Path strings that do not need to be converted to CWL File objects
        for file_string in ['run_directory', 'tmp_directory', 'ebwt']:
            path = self.paths[file_string]
            self[file_string] = self.paths.from_here(path) if path else ''

    def process_sample_sheet(self):
        sample_sheet = self.paths.from_here(self['samples_csv']['path'])
        sample_sheet_dir = os.path.dirname(sample_sheet)
        groups_reps = Counter()

        csv_reader = CSVReader(sample_sheet, "Samples Sheet")
        for row in csv_reader.rows():
            if not os.path.splitext(row['File'])[1] in [".fastq", ".gz"]:
                raise ValueError("Files in samples.csv must have a .fastq(.gz) extension:\n%s" % (row['File'],))

            fastq_file = self.from_here(row['File'], origin=sample_sheet_dir)
            sample_basename = self.prefix(os.path.basename(fastq_file))

            group_name = row['Group']
            rep_number = row['Replicate']
            groups_reps[group_name] += 1

            self.append_to('sample_basenames', sample_basename)
            self.append_to('fastp_report_titles', f"{group_name}_rep_{rep_number}")
            if row['Control'].lower() == 'true':
                self['control_condition'] = group_name

            try:
                self.append_to('in_fq', self.cwl_file(fastq_file))
            except FileNotFoundError:
                line = csv_reader.line_num
                sys.exit("The fastq file on line %d of your Samples Sheet was not found:\n%s" % (line, fastq_file))

        self.check_deseq_compatibility(groups_reps)

    def check_deseq_compatibility(self, sample_groups):
        total_samples = sum(sample_groups.values())
        total_coefficients = len(sample_groups)
        degrees_of_freedom = total_samples - total_coefficients

        if degrees_of_freedom < 1:
            self['run_deseq'] = False
            print("Your experiment design has less than one degree of freedom, which is incompatible "
                  "with DESeq2. The DGE step will be skipped and most plots will not be produced.",
                  file=sys.stderr)

    def process_feature_sheet(self):
        # Configuration doesn't currently do anything with Feature Sheet contents
        return

    def setup_file_groups(self):
        """Configuration keys that represent lists of files"""

        self.set_default_dict({per_file_setting_key: [] for per_file_setting_key in
            ['in_fq', 'sample_basenames', 'gff_files', 'fastp_report_titles']
        })
            
    def setup_pipeline(self):
        """Overall settings for the whole pipeline"""

        self.dt = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self['run_date'], self['run_time'] = self.dt.split('_')

        default_run_name = '_'.join(x for x in [self['user'], "tinyrna"] if x)
        self['run_name'] = self.get('run_name', default=default_run_name) + "_" + self.dt

        # Create prefixed Run Directory name
        run_dir_parent, run_dir = os.path.split(self['run_directory'].rstrip(os.sep))
        self['run_directory'] = self.joinpath(run_dir_parent, self['run_name'] + run_dir)

        self.templates = resource_filename('tiny', 'templates/')

    def setup_ebwt_idx(self):
        """Bowtie index files and prefix"""

        # Determine if bowtie-build should run, and set Bowtie index prefix accordingly
        bt_index_prefix = self.paths['ebwt']
        if self['run_bowtie_build'] and not bt_index_prefix:
            if not self['reference_genome_files']:
                raise ValueError(f"If {self.basename} contains 'run_bowtie_build: True', you "
                                 f"need to provide your reference genome files in {self.paths.basename}")

            # Outputs are saved in {run_directory}/bowtie-build, within which prefix is first genome file's basename
            first_genome_file = self.paths.from_here(self['reference_genome_files'][0]['path'])
            bt_index_prefix = self.prefix(os.path.join(
                self['run_directory'], "bowtie-build", os.path.basename(first_genome_file))
            )

            self['ebwt'] = self.paths['ebwt'] = bt_index_prefix
        else:
            # bowtie-build should only run if 'run_bowtie_build' is True AND ebwt (index prefix) is undefined
            self['run_bowtie_build'] = False
            bt_index_prefix = self.paths.from_here(bt_index_prefix)

        # Bowtie index files
        try:
            self['bt_index_files'] = [self.cwl_file(bt_index_prefix + postfix, verify=(not self['run_bowtie_build']))
                            for postfix in ['.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt']]
        except FileNotFoundError as e:
            sys.exit("The following file could not be found from the Bowtie index prefix defined in your Paths File:\n"
                     "%s" % (e.filename,))

        # When CWL copies bt_index_filex for the bowtie.cwl InitialWorkDirRequirement, it does not
        # preserve the prefix path. What the workflow "sees" is the ebwt files at working dir root
        self["ebwt"] = os.path.basename(self["ebwt"])

    def validate_inputs(self):
        """For now, only GFF files are validated here"""

        gff_files = {gff['path']: [] for gff in self['gff_files']}
        prefs = {x: self[f'{x}_filter'] for x in ['source', 'type']}
        ebwt = self.paths['ebwt'] if not self['run_bowtie_build'] else None

        GFFValidator(
            gff_files,
            prefs=prefs,
            ebwt=ebwt,
            genomes=self.paths['reference_genome_files'],
            alignments=None  # Used in tiny-count standalone runs
        ).validate()

    def save_run_profile(self, config_file_name=None) -> str:
        """Saves Samples Sheet and processed run config to the Run Directory for record keeping"""

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


class CSVReader(csv.DictReader):
    """A simple wrapper class for csv.DictReader

    This makes field labels consistent across the project and simplifies the code
    """

    # user-facing name -> internal short name
    tinyrna_sheet_fields = {
        "Features Sheet": OrderedDict({
           "Select for...":     "Key",
           "with value...":     "Value",
           "Tag":               "Tag",
           "Hierarchy":         "Hierarchy",
           "Strand":            "Strand",
           "5' End Nucleotide": "nt5end",
           "Length":            "Length",
           "Overlap":           "Overlap",
        }),
        "Samples Sheet": OrderedDict({
            "Input FASTQ Files": "File",
            "Sample/Group Name": "Group",
            "Replicate Number":  "Replicate",
            "Control":           "Control",
            "Normalization":     "Normalization"
        })
    }

    def __init__(self, filename: str, doctype: str = None):
        self.doctype = doctype
        self.tinyrna_file = filename
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
                yield row

    def validate_csv_header(self, header: OrderedDict):
        doc_reference = self.tinyrna_sheet_fields[self.doctype]
        expected = {key.lower() for key in doc_reference.keys()}
        read_vals = {val.lower() for val in header.values() if val is not None}

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


if __name__ == '__main__':
    Configuration.main()