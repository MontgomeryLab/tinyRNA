import ruamel.yaml
import argparse
import csv
import os
import re

from ruamel.yaml.comments import CommentedOrderedMap
from pkg_resources import resource_filename
from datetime import datetime
from typing import Union, Any


class ConfigBase:
    """Base class for basic aquatx configuration operations

    Attributes:
        yaml: the YAML interface for reading config and writing processed config
        config: the configuration object produced by loading the config file
        inf: the filename of the configuration .yml file to process
        dir: parent directory of the input file. Used for calculating paths relative to config file.
        extras: path to the package extras directory
        dt: a date-time string for default output naming
    """

    def __init__(self, config_file: str):
        self.dir = os.path.dirname(os.path.abspath(config_file))
        self.inf = config_file
        self.basename = os.path.basename(config_file)
        self.extras = ''
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
        if os.path.isabs(path2): return path2
        return os.path.normpath(os.path.join(path1, path2))

    @staticmethod
    def cwl_file(file: str) -> dict:
        """Returns a minimal File object as defined by CWL"""
        # Todo: validate that file exists at this step
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

    def get_outfile_path(self, infile: str) -> str:
        """Prepend date+time to the Config File name. Change its path to reside in the Run Directory."""
        return self.joinpath(self['run_directory'], os.path.basename(infile))

    def write_processed_config(self, filename: str = None) -> str:
        """Writes the current configuration """
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

    def __init__(self, config_file: str):
        # Parse YAML configuration file
        super().__init__(config_file)

        self.paths = self.load_paths_config()
        self.process_paths_sheet()
        
        self.setup_pipeline()
        self.setup_per_file()
        self.setup_ebwt_idx()
        self.process_sample_sheet()
        self.process_feature_sheet()
        
    def load_paths_config(self):
        """Constructs a sub-configuration object containing workflow file preferences"""
        path_sheet = self.from_here(self['paths_config'])
        return ConfigBase(path_sheet)

    def process_paths_sheet(self):
        """Loads the paths of all related config files and workflow inputs"""

        def to_cwl_file_class(input_file_path):
            path_to_input = self.paths.from_here(input_file_path)
            return self.cwl_file(path_to_input)

        self['ebwt'] = self.paths['ebwt']
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

        with open(sample_sheet, 'r', encoding='utf-8-sig') as sf:
            fieldnames = ("File", "Group", "Replicate")
            csv_reader = csv.DictReader(sf, fieldnames=fieldnames, delimiter=',')

            next(csv_reader)  # Skip header line
            for row in csv_reader:
                if not os.path.splitext(row['File'])[1] in [".fastq", ".gz"]:
                    raise ValueError("Files in samples.csv must have a .fastq(.gz) extension:\n%s" % (row['File'],))
                fastq_file = self.from_here(row['File'], origin=sample_sheet_dir)
                sample_basename = self.prefix(os.path.basename(fastq_file))
                group_name = row['Group']
                rep_number = row['Replicate']

                self.append_to('report_title', f"{group_name}_rep_{rep_number}")
                self.append_to('in_fq', self.cwl_file(fastq_file))

                self.append_to('out_fq', sample_basename + '_cleaned.fastq')
                self.append_to('outfile', sample_basename + '_aligned_seqs.sam')
                self.append_to('logfile', sample_basename + '_alignment_log.log')
                self.append_to('un', sample_basename + '_unaligned_seqs.fa')
                self.append_to('json', sample_basename + '_qc.json')
                self.append_to('html', sample_basename + '_qc.html')
                self.append_to('uniq_seq_prefix', sample_basename)

    def process_feature_sheet(self):
        feature_sheet = self.paths.from_here(self['features_csv']['path'])
        feature_sheet_dir = os.path.dirname(feature_sheet)

        with open(feature_sheet, 'r', encoding='utf-8-sig') as ff:
            fieldnames = ("ID", "Key", "Value", "Hierarchy", "Strand", "nt5", "Length", "Strict", "Source")
            csv_reader = csv.DictReader(ff, fieldnames=fieldnames, delimiter=',')

            next(csv_reader) # Skip header line
            for row in csv_reader:
                gff_file = self.from_here(row['Source'], origin=feature_sheet_dir)
                self.append_if_absent('gff_files', self.cwl_file(gff_file))
            
    def setup_per_file(self):
        """Per-library settings lists to be populated by entries from samples_csv"""

        self.set_default_dict({per_file_setting_key: [] for per_file_setting_key in
            ['un', 'in_fq', 'out_fq', 'uniq_seq_prefix', 'gff_files', 'outfile', 'logfile', 'report_title', 'json', 'html']
        })
            
    def setup_pipeline(self):
        """Overall settings for the whole pipeline"""

        self.dt = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.set_default_dict({
            'run_date': self.dt.split('_')[0],
            'run_time': self.dt.split('_')[1]
        })

        default_run_name = '_'.join(x for x in [self['user'], "aquatx"] if x)
        self['run_name'] = self.get('run_name', default=default_run_name) + "_" + self.dt

        # Create prefixed Run Directory name
        run_dir_resolved = self.paths.from_here(self.get('run_directory', default='run_directory'))
        run_dir_parent = os.path.dirname(run_dir_resolved)
        run_dir_withdt = self['run_name'] + '_' + os.path.basename(run_dir_resolved)
        self['run_directory'] = self.joinpath(run_dir_parent, run_dir_withdt)

        self.extras = resource_filename('aquatx', 'extras/')

    # Todo: better heuristics for determining if the prefix outputs actually exist
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

            # Finally, update user's paths file with the new prefix
            self.paths.write_processed_config(self.paths.inf)
        else:
            # bowtie-build should only run if 'run_bowtie_build' is True AND ebwt (index prefix) is undefined
            self['run_bowtie_build'] = False
            bt_index_prefix = self.paths.from_here(bt_index_prefix)

        # Bowtie index files
        self['bt_index_files'] = [self.cwl_file(bt_index_prefix + postfix)
                        for postfix in ['.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt']]

        # When CWL copies bt_index_filex for the bowtie.cwl InitialWorkDirRequirement, it does not
        # preserve the prefix path. What the workflow "sees" is the ebwt files at working dir root
        self["ebwt"] = os.path.basename(self["ebwt"])

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
        Configuration(args.input_file).write_processed_config()

    if __name__ == '__main__':
        main()


class ResumeConfig(ConfigBase):
    """A class for modifying the workflow and config to resume a run at Counter

    The modified workflow document is written to the package resource directory.
    Directory output names have a timestamp appended to them to keep outputs
    separate between runs. Timestamped run prefixes are also updated to the
    date and time of the resume run. A copy of the Features Sheet is included
    in the counts output directory to maintain a record of the selection rules
    that were used for the resume run.
    """

    def __init__(self, processed_config, workflow):
        # Parse the pre-processed YAML configuration file
        super().__init__(processed_config)

        # Load the CWL workflow YAML for modification
        with open(workflow, 'r') as f:
            self.workflow: CommentedOrderedMap
            self.workflow = self.yaml.load(f)

        self.dt = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.steps = ["counter", "dge"]
        self.new_counter_inputs = {
            'aligned_seqs': "resume_sams",
            'fastp_logs': "resume_fastp_logs",
            'collapsed_fa': "resume_collapsed_fas"
        }

        self.create_truncated_workflow()
        self.rebuild_counter_inputs()
        self.add_timestamps()

    def rebuild_counter_inputs(self):
        """Set the new path inputs for the counts step

        Normally, the counts step receives inputs from the current run's previous step outputs.
        When resuming, these outputs should already exist on disk. We need to populate the new
        File[] arrays with their corresponding pipeline outputs on disk.
        """

        def cwl_file_resume(subdir, file):
            return self.cwl_file('/'.join([subdir, file]))

        self['resume_sams'] = [cwl_file_resume(self['dir_name_bowtie'], sam) for sam in self['outfile']]
        self['resume_fastp_logs'] = [cwl_file_resume(self['dir_name_fastp'], log) for log in self['json']]
        self['resume_collapsed_fas'] = [cwl_file_resume(self['dir_name_collapser'], prefix + "_collapsed.fa")
                                        for prefix in self['uniq_seq_prefix']]

    def create_truncated_workflow(self):
        """Modifies the workflow to begin running at the counts step

        New input variables must be created, which point to the files produced by the
        prior pipeline run. The counts step is modified to receive these inputs directly
        rather than expecting CWL to run the previous workflow steps. Parameters for Counter outputs and
        DEG inputs do not require modification. The subdirs step is removed and replaced with modified
        copies of its own calls to make-subdir.cwl, for the creation of Counter and DEG outputs.
        The Features Sheet is also included in the new Counter folder.
        """

        # Remove incompatible steps
        for step in ['bt_build_optional', 'counts-prep', 'subdirs']:
            del self.workflow['steps'][step]

        # Setup new inputs
        for inp, new_input in self.new_counter_inputs.items():
            # Create WorkflowInputParameters for new "resume_" File arrays
            self.workflow['inputs'][new_input] = "File[]"
            # Update WorkflowStepInputs for these new variables in counts
            self.workflow['steps']['counts']['in'][inp] = new_input

        # Remove subdir WorkflowOutputParameters
        for outdir in ['bt_build', 'fastp', 'collapser', 'bowtie', 'counter', 'dge']:
            del self.workflow['outputs'][f'{outdir}_out_dir']

        # Replace the subdirs subworkflow step with calls to its underlying CommandLineTool.
        with open(resource_filename('aquatx', 'cwl/workflows/organize-outputs.cwl')) as f:
            organizer_sub_wf = self.yaml.load(f)

            sources = {"counter":
                           ["counts/feature_counts", "counts/other_counts", "counts/alignment_stats",
                            "counts/summary_stats", "counts/intermed_out_files", "features_csv"],
                       "dge":
                            ['dge/norm_counts', 'dge/comparisons']}

            for step in self.steps:
                step_name = f'organize_{step}'
                # Copy relevant steps from organize-outputs.cwl
                context = self.workflow['steps'][step_name] = organizer_sub_wf['steps'][step_name]
                # Update WorkflowStepInputs
                context['in']['dir_name'] = f'dir_name_{step}'
                context['in']['dir_files'] = {'source': sources[step], 'pickValue': 'all_non_null',
                                              'linkMerge': 'merge_flattened'}
                # Update WorkflowOutputParameter
                self.workflow['outputs'][f'{step}_out_dir'] = {
                    'type': "Directory",
                    'outputSource': f'{step_name}/subdir'
                }

    def add_timestamps(self):
        """Adds timestamps to subdirectory outputs, as well as a copy of Features Sheet

        This is done for record keeping and to keep outputs organized between runs.
        """

        # Rename output directories with timestamp
        for subdir in self.steps:
            step_dir = "dir_name_" + subdir
            self[step_dir] = self[step_dir] + "_" + self.dt

        # Update run_name output prefix variable for the current date and time
        self['run_name'] = re.sub(r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}", self.dt, self['run_name'])


    def write_workflow(self, workflow_outfile: str) -> None:
        with open(workflow_outfile, "w") as wf:
            self.yaml.dump(self.workflow, wf)
