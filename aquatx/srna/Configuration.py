import ruamel.yaml
import argparse
import csv
import os
import sys

from pkg_resources import resource_filename
from datetime import datetime
from shutil import copyfile
from typing import Union

from ruamel.yaml.comments import CommentedMap


class Configuration:
    def __init__(self, input_file: str):
        self.dir = os.path.dirname(input_file) + os.sep
        self.inf = input_file

        # Parse YAML run configuration file
        self.yaml = ruamel.yaml.YAML()
        with open(input_file, 'r') as conf:
            self.config: CommentedMap = self.yaml.load(conf)

        self.setup()
        self.process_sample_sheet()
        self.process_reference_sheet()

    def process_sample_sheet(self):
        sample_sheet = self.joinpath(self.dir, self.get('samples_csv'))
        from_here = os.path.dirname(sample_sheet)

        with open(sample_sheet, 'r', encoding='utf-8-sig') as sf:
            csv_reader = csv.DictReader(sf, delimiter=',')
            ext = 'gz' if self.get('compress') else ''
            for row in csv_reader:
                sample_basename = self.prefix(os.path.basename(row['Input FastQ/A Files']))
                fastq_file = self.joinpath(from_here, row['Input FastQ/A Files'])
                group_name = row['Sample/Group Name']
                rep_number = row['Replicate number']

                self.append_to('report_title', f"{group_name}_replicate_{rep_number}_fastp_report")
                self.append_to('out_prefix', f"{group_name}_replicate_{rep_number}")
                self.append_to('in_fq', self.cwl_file(fastq_file))

                self.append_to('out_fq', sample_basename + '_cleaned.fastq')
                self.append_to('outfile', sample_basename + '_aligned_seqs.sam')
                self.append_to('un', sample_basename + '_unaligned_seqs.fa')
                self.append_to('json', sample_basename + '_qc.json')
                self.append_to('html', sample_basename + '_qc.html')
                self.append_to('uniq_seq_prefix', sample_basename)

    def process_reference_sheet(self):
        reference_sheet = self.joinpath(self.dir, self.get('features_csv'))
        from_here = os.path.dirname(reference_sheet)

        with open(reference_sheet, 'r', encoding='utf-8-sig') as rf:
            csv_reader = csv.DictReader(rf, delimiter=',')
            for row in csv_reader:
                gff_file = self.joinpath(from_here, row['Feature Source'])
                self.append_to('identifier', row['Identifier'])
                self.append_to('srna_class', row['Class'])
                self.append_to('strand', row['Strand (sense/antisense/both)'])
                self.append_to('ref_annotations', self.cwl_file(gff_file))
                self.append_to('hierarchy', row['Hierarchy'])
                self.append_to('5end_nt', row["5' End Nucleotide"])
                self.append_to('length', row['Length'])

    def setup(self):
        """Populates default values and prepares per-file configuration lists"""

        self.dt = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        default_prefix = '_'.join(x for x in [self.dt, self.get('user'), "aquatx"] if x)
        self.set_default_dict({
            'run_directory': default_prefix,
            'run_prefix': default_prefix,
            'output_prefix': default_prefix,
            'run_date': self.dt.split('_')[0],
            'run_time': self.dt.split('_')[1]
        })

        # Per-file settings lists to be populated by process_sample_sheet() and process_reference_sheet()
        self.set_default_dict({per_file_setting_key: [] for per_file_setting_key in
            ['identifier', 'srna_class', 'strand', 'hierarchy', '5end_nt', 'length', 'ref_annotations', 'un',
             'in_fq', 'out_fq', 'uniq_seq_prefix', 'out_prefix', 'outfile', 'report_title', 'json', 'html']
        })

        self.extras = resource_filename('aquatx', 'extras/')
        self.set('output_file_stats', self.get('output_prefix') + '_run_stats.csv')
        self.set('output_file_counts', self.get('output_prefix') + '_raw_counts.csv')

        # Bowtie index file prefix
        bt_idx = (self.set('ebwt', self.prefix(self.get('ref_genome')))
                  if self.get('run_idx') and not self.get('ebwt')
                  else self.get('ebwt'))

        # Bowtie index files
        bt_idx_files = [self.cwl_file(bt_idx + postfix)
                        for postfix in ['.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt']]

        self.set('bt_index_files', bt_idx_files)

        # Discard these configurations if they are set as such
        if self.get('adapter_sequence') == 'auto_detect':
            self.config.pop('adapter_sequence')

    """========== GETTERS AND SETTERS =========="""

    def get(self, key: str) -> Union[str, list, dict, None]:
        return self.config.get(key, None)

    def set(self, key: str, val: Union[str, list, dict]) -> Union[str, list, dict]:
        self.config[key] = val
        return val

    def set_if_not(self, key: str, val: Union[str, list, dict]) -> Union[str, list, dict]:
        """Apply the setting if it has not been previously set"""
        if not self.get(key):
            return self.set(key, val)
        else:
            return self.get(key)

    def set_default_dict(self, setting_dict: dict) -> None:
        """Apply each setting in the input dictionary if it has not been previously set"""
        for key, val in setting_dict.items():
            self.set_if_not(key, val)

    def append_to(self, key: str, val: Union[str, list, dict]) -> list:
        """Append a file setting to a per-file settings list"""
        target = self.get(key)
        if type(target) is list:
            target.append(val)
            return target
        else:
            print(f"Tried appending to a non-existent key: {key}", file=sys.stderr)

    """========== HELPERS =========="""

    def cwl_file(self, file: str) -> dict:
        """Returns a file input/output specification for the CWL config"""
        return {'class': 'File', 'path': file}

    def prefix(self, path: str) -> str:
        """Returns everything from path except the file extension"""
        return os.path.splitext(path)[0]

    def joinpath(self, path1: str, path2: str) -> str:
        """Combines two relative paths intelligently."""
        if os.path.isabs(path2): return path2
        return os.path.normpath(os.path.join(path1, path2))

    def get_outfile_name(self, infile: str) -> str:
        """This will likely be changed in the near future"""
        if os.path.basename(infile) == 'run_config_template.yml':
            output_name = self.dt + '_run_config.yml'
            copyfile(infile, output_name)
            return output_name
        else:
            return infile

    def write_processed_config(self, filename: str = None) -> str:
        if filename is None: filename = self.get_outfile_name(self.inf)

        with open(filename, 'w') as outconf:
            self.yaml.dump(self.config, outconf)

        return filename

    def create_run_directory(self) -> str:
        run_dir = self.get("run_directory")
        if not os.path.isdir(run_dir):
            os.mkdir(run_dir)

        return run_dir

    """========== COMMAND LINE =========="""

    def get_args(self):
        """Get the input arguments"""

        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input-file', metavar='CONFIG', required=True,
                            help="Input file")

        args = parser.parse_args()
        return args

    def main(self):
        """
        Main routine to process the run information.
        """

        # Get input config file
        args = self.get_args()
        Configuration(args.input_file)

        # TODO: need to specify the non-model organism run when no reference genome is given

    if __name__ == '__main__':
        main()