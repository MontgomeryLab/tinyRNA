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
        self.setup_per_file()
        self.setup_ebwt_idx()
        self.process_sample_sheet()
        self.process_feature_sheet()
        if validate_inputs: self.validate_inputs()
        
    def load_paths_config(self):
        """Constructs a sub-configuration object containing workflow file preferences"""
        path_sheet = self.from_here(self['paths_config'])
        return ConfigBase(path_sheet)

    def process_paths_sheet(self):
        """Loads the paths of all related config files and workflow inputs"""

        def to_cwl_file_class(input_file_path):
            path_to_input = self.paths.from_here(input_file_path)
            return self.cwl_file(path_to_input)

        for absorb_key in ['ebwt', 'plot_style_sheet', 'adapter_fasta', 'tmp_directory']:
            self[absorb_key] = self.paths[absorb_key]
        self['run_directory'] = self.paths.from_here(self.paths['run_directory'])

        # Configurations that need to be converted from string to a CWL File object
        self['samples_csv'] = to_cwl_file_class(self.paths.from_here(self.paths['samples_csv']))
        self['features_csv'] = to_cwl_file_class(self.paths.from_here(self.paths['features_csv']))
        self['reference_genome_files'] = [
            to_cwl_file_class(self.paths.from_here(genome))
            for genome in self.paths['reference_genome_files']
        ]

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
        feature_sheet = self.paths.from_here(self['features_csv']['path'])
        feature_sheet_dir = os.path.dirname(feature_sheet)

        csv_reader = CSVReader(feature_sheet, "Features Sheet")
        for row in csv_reader.rows():
            gff_file = self.from_here(row['Source'], origin=feature_sheet_dir)
            try:
                self.append_if_absent('gff_files', self.cwl_file(gff_file))
            except FileNotFoundError:
                line = csv_reader.line_num
                sys.exit("The GFF file on line %d of your Features Sheet was not found:\n%s" % (line, gff_file))

    def setup_per_file(self):
        """Per-library settings lists to be populated by entries from samples_csv"""

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
        run_dir_resolved = self.paths.from_here(self.get('run_directory', default='run_directory'))
        run_dir_parent = os.path.dirname(run_dir_resolved)
        run_dir_withdt = self['run_name'] + '_' + os.path.basename(run_dir_resolved)
        self['run_directory'] = self.joinpath(run_dir_parent, run_dir_withdt)

        self.templates = resource_filename('tiny', 'templates/')

    def setup_ebwt_idx(self):
        """Determines Bowtie index prefix and whether bowtie-build should run"""

        # Empty values for ebwt (''/~/None) trigger bowtie-build
        self['run_bowtie_build'] = not bool(self.paths['ebwt'])

        if self['run_bowtie_build']:
            # Set the prefix to the run directory outputs. This is necessary
            # because workflow requires bt_index_files to be a populated list.
            prefix = self.get_ebwt_prefix()
        else:
            prefix = self.paths.from_here(self.paths['ebwt'])

        # verify_bowtie_build_outputs() will check if these end up being long indexes
        self['bt_index_files'] = self.get_bt_index_files(prefix)

        # When CWL copies bt_index_filex for the bowtie.cwl InitialWorkDirRequirement, it does not
        # preserve the prefix path. What the workflow "sees" is the ebwt files at working dir root
        self["ebwt"] = os.path.basename(prefix)

    def get_ebwt_prefix(self):
        """Determines the output prefix path for bowtie indexes that haven't been built yet. The basename
        of the prefix path is simply the basename of the reference genome sans file extension"""

        genome_files = self['reference_genome_files']
        if not genome_files:
            raise ValueError("If your Paths Sheet doesn't have a value for \"ebtw:\", then bowtie indexes "
                             "will be built, but you'll need to provide your reference genome files under "
                             '"reference_genome_files:" (also in your Paths Sheet)')

        genome_basename = os.path.basename(genome_files[0]['path'])
        return self.prefix(os.path.join( # prefix path:
            self['run_directory'], self['dir_name_bt_build'], genome_basename
        ))

    def get_bt_index_files(self, prefix):
        """Builds the list of expected bowtie index files from the ebwt prefix. If an index file
        doesn't exist then they will be automatically rebuilt from the user's reference genomes.
        File existence isn't checked if bowtie-build is already scheduled for this run."""

        try:
            verify_file_paths = not bool(self['run_bowtie_build'])
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
                new_prefix = self.get_ebwt_prefix()
                self['run_bowtie_build'] = True
                return self.get_bt_index_files(new_prefix)
            else:
                sys.exit(problem + userfix)

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
           "Alias by...":       "Name",
           "Tag":               "Tag",
           "Hierarchy":         "Hierarchy",
           "Strand":            "Strand",
           "5' End Nucleotide": "nt5end",
           "Length":            "Length",
           "Overlap":           "Overlap",
           "Feature Source":    "Source"
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