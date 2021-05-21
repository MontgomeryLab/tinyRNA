import ruamel.yaml
import argparse
import csv
import os

from pkg_resources import resource_filename
from datetime import datetime
from shutil import copyfile
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
        self.dir = os.path.dirname(config_file)
        self.inf = config_file
        self.extras = ''
        self.dt = ''

        self.yaml = ruamel.yaml.YAML()
        with open(config_file, 'r') as f:
            self.config = self.yaml.load(f)

    def __getitem__(self, key: str) -> Any:
        return self.get(key)

    def get(self, key: str) -> Any:
        return self.config.get(key, None)

    def set(self, key: str, val: Union[str, list, dict, bool]) -> Union[str, list, dict, bool]:
        self.config[key] = val
        return val

    def set_if_not(self, key: str, val: Union[str, list, dict, bool]) -> Any:
        """Apply the setting if it has not been previously set"""
        if not self.get(key):
            return self.set(key, val)
        else:
            return self.get(key)

    def set_default_dict(self, setting_dict: dict) -> None:
        """Apply each setting in setting_dict if it has not been previously set"""
        for key, val in setting_dict.items():
            # Can't use config.setdefault(), it considers None and [] "already set"
            self.set_if_not(key, val)

    def append_to(self, key: str, val: Any) -> list:
        """Append a list-type setting (per-file settings)"""
        target = self.get(key)
        if type(target) is list:
            target.append(val)
            return target
        else:
            raise ValueError(f"Tried appending to a non-existent key: {key}")

    def append_if_absent(self, key: str, val: Any) -> list:
        """Append to list-type setting if the value is not already present"""
        target = self.get(key)
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
        """Returns a file input/output specification for the CWL config"""
        # Todo: validate that file exists at this step
        return {'class': 'File', 'path': file}

    @staticmethod
    def cwl_dir(dir: str) -> dict:
        """Returns a directory specificaion for the CWL config"""
        return {'class': 'Directory', 'path': dir}

    def create_run_directory(self) -> str:
        """Create the destination directory for pipeline outputs"""
        run_dir = self.get("run_directory")
        if not os.path.isdir(run_dir):
            os.mkdir(run_dir)

        return run_dir

    def get_outfile_name(self, infile: str) -> str:
        """If the user's config file was named run_config_template.yml, copy and rename.
        This will likely be changed in the near future"""
        if os.path.basename(infile) == 'run_config_template.yml':
            output_name = self.dt + '_run_config.yml'
            copyfile(infile, output_name)
            return output_name
        else:
            return infile

    def write_processed_config(self, filename: str = None) -> str:
        """Writes the current configuration """
        if filename is None: filename = self.get_outfile_name(self.inf)

        with open(filename, 'w') as outconf:
            self.yaml.dump(self.config, outconf)

        return filename


class Configuration(ConfigBase):
    """A class for processing and updating a YAML config file for CWL

    Ultimately, this class populates pipeline settings and per-file settings for pipeline steps.
    Per-file settings are determined by samples.csv, features.csv, and the input config file.
    Paths provided in these three files are evaluated relative to the given file.
    Absolute paths may also be supplied.

    Attributes:
        paths: the configuration object from processing the paths_config file.
            This holds path info for other config files and prefixes, and is updated
            appropriately if 'run_bowtie_index' is set to 'true'
    """

    def __init__(self, config_file: str):
        # Parse YAML run configuration file
        super().__init__(config_file)

        self.paths = self.load_paths_config()
        self.process_paths_sheet()
        
        self.setup_pipeline()
        self.setup_per_file()
        self.setup_ebwt_idx()
        self.process_sample_sheet()
        self.process_feature_sheet()
        
    def load_paths_config(self):
        path_sheet = self.joinpath(self.dir, self.get('paths_config'))
        return ConfigBase(path_sheet)

    def process_paths_sheet(self):
        from_here = self.paths.dir

        def to_cwl_file_class(input_file_path):
            return self.cwl_file(
                self.joinpath(from_here, input_file_path))

        # No conversion required, simply copy
        self.set('ebwt', self.paths['ebwt'])
        self.set('run_directory', self.paths['run_directory'])

        # Configurations that need to be converted from string to a CWL File class object
        self.set('samples_csv', to_cwl_file_class(self.paths['samples_csv']))
        self.set('features_csv', to_cwl_file_class(self.paths['features_csv']))
        self.set('reference_genome_files', [to_cwl_file_class(genome) for genome in self.paths['reference_genome_files']])

    def process_sample_sheet(self):
        sample_sheet = self.joinpath(self.dir, self.get('samples_csv')['path'])
        from_here = os.path.dirname(sample_sheet)

        with open(sample_sheet, 'r', encoding='utf-8-sig') as sf:
            fieldnames = ("File", "Group", "Replicate")
            csv_reader = csv.DictReader(sf, fieldnames=fieldnames, delimiter=',')

            next(csv_reader)  # Skip header line
            for row in csv_reader:
                if not os.path.splitext(row['File'])[1] in [".fastq", ".gz"]:
                    raise ValueError("Files in samples.csv must have a .fastq(.gz) extension:\n%s" % (row['File'],))
                sample_basename = self.prefix(os.path.basename(row['File']))
                fastq_file = self.joinpath(from_here, row['File'])
                group_name = row['Group']
                rep_number = row['Replicate']

                self.append_to('report_title', f"{group_name}_rep_{rep_number}")
                self.append_to('in_fq', self.cwl_file(fastq_file))

                self.append_to('out_fq', sample_basename + '_cleaned.fastq')
                self.append_to('outfile', sample_basename + '_aligned_seqs.sam')
                self.append_to('un', sample_basename + '_unaligned_seqs.fa')
                self.append_to('json', sample_basename + '_qc.json')
                self.append_to('html', sample_basename + '_qc.html')
                self.append_to('uniq_seq_prefix', sample_basename)

    def process_feature_sheet(self):
        feature_sheet = self.joinpath(self.dir, self.get('features_csv')['path'])
        from_here = os.path.dirname(feature_sheet)

        with open(feature_sheet, 'r', encoding='utf-8-sig') as ff:
            fieldnames = ("ID", "Key", "Value", "Hierarchy", "Strand", "nt5", "Length", "Strict", "Source")
            csv_reader = csv.DictReader(ff, fieldnames=fieldnames, delimiter=',')

            next(csv_reader) # Skip header line
            for row in csv_reader:
                gff_file = self.joinpath(from_here, row['Source'])
                self.append_if_absent('gff_files', self.cwl_file(gff_file))
            
    def setup_per_file(self):
        """Per-file settings lists to be populated by entries from samples_csv"""

        self.set_default_dict({per_file_setting_key: [] for per_file_setting_key in
            ['un', 'in_fq', 'out_fq', 'uniq_seq_prefix', 'gff_files', 'outfile', 'report_title', 'json', 'html']
        })
            
    def setup_pipeline(self):
        """Overall settings for the whole pipeline"""

        self.dt = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        default_prefix = '_'.join(x for x in [self.dt, self.get('user'), "aquatx"] if x)
        self.set_default_dict({
            'run_directory': default_prefix,
            'run_prefix': default_prefix,
            'output_prefix': default_prefix,
            'run_date': self.dt.split('_')[0],
            'run_time': self.dt.split('_')[1]
        })

        self.extras = resource_filename('aquatx', 'extras/')

    def setup_ebwt_idx(self):
        """Bowtie index files and prefix"""

        # Determine if bowtie-build should run, and set Bowtie index prefix accordingly
        bt_index_prefix = self.paths.get('ebwt')
        if self.get('run_bowtie_build') and bt_index_prefix == '':
            if not self.get('reference_genome_files'):
                config_file = os.path.basename(self.inf)
                paths_file = os.path.basename(self.paths.inf)
                raise ValueError(f"If {config_file} contains 'run_bowtie_build: true', you "
                                 f"need to provide your reference genome files in {paths_file}")

            # Outputs are saved in run_directory, prefix is simply the first genome file's basename
            ref_genome_base = os.path.basename(self.get('reference_genome_files')[0]['path'])
            bt_index_prefix = self.joinpath(self.get('run_directory'), ref_genome_base)
            self.set('ebwt', bt_index_prefix)

            # Finally, update user's paths file with the new prefix
            self.paths.set('ebwt', bt_index_prefix)
            self.paths.write_processed_config()
        else:
            # bowtie-build should only run if 'run_bowtie_build' is True AND ebwt is ''
            self.set('run_bowtie_build', False)

        # Bowtie index files
        self.set('bt_index_files', [self.cwl_file(bt_index_prefix + postfix)
                        for postfix in ['.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt']])

        # When CWL copies bt_index_filex for the bowtie.cwl InitialWorkDirRequirement, it does not
        # preserve the prefix path. What the workflow "sees" is the ebwt files at working dir root
        self.set("ebwt", os.path.basename(self.get("ebwt")))

    """========== COMMAND LINE =========="""

    @staticmethod
    def main():
        """Main routine to process the run information."""

        # Get input config file
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input-file', metavar='CONFIG', required=True,
                            help="Input file")

        args = parser.parse_args()
        Configuration(args.input_file).write_processed_config()

        # TODO: need to specify the non-model organism run when no reference genome is given

    if __name__ == '__main__':
        main()
