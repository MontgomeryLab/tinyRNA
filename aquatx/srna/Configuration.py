import ruamel.yaml
import argparse
import sys
import csv
import os

from typing import Union
from shutil import copyfile
from datetime import datetime
from pkg_resources import resource_filename

class Configuration:
    def __init__(self, input_file: str):
        self.dir = os.path.dirname(input_file) + os.sep

        # Parse YAML run configuration file
        with open(input_file) as conf:
            self.config: dict = ruamel.yaml.YAML().load(conf)

        self.setup()
        self.process_sample_sheet()
        self.process_reference_sheet()



    def process_sample_sheet(self):
        sample_sheet = self.dir + self.get('sample_sheet_file')
        from_here = os.path.dirname(sample_sheet) + os.sep

        with open(sample_sheet, 'r') as sf:
            csv_reader = csv.DictReader(sf, delimiter=',')
            for row in csv_reader:
                sample_basename = self.prefix(os.path.basename(row['Input FastQ/A Files']))
                group_name = row['Sample/Group Name']
                rep_number = row['Replicate number']

                self.append_to('report_title', f"{group_name}_replicate_{rep_number}_fastp_report'")
                self.append_to('out_prefix', f"{group_name}_replicate_{rep_number}")
                self.append_to('in_fq', self.cwl_file_def(f"{from_here}{row['Input FastQ/A Files']}"))

                self.append_to('out_fq', sample_basename + '_cleaned.fastq')
                self.append_to('uniq_seq_file', sample_basename + '_unique_seqs_collapsed.fa')
                self.append_to('keep_low_counts', sample_basename + '_low_count_uniq_seqs.fa')
                self.append_to('outfile', sample_basename + '_aligned_seqs.sam')
                self.append_to('un', sample_basename + '_unaligned_seqs.fa')
                self.append_to('json', sample_basename + '_qc.json')
                self.append_to('html', sample_basename + '_qc.html')

        if not self.get('keep_low_counts'):
            self.config.pop('keep_low_counts')

    def process_reference_sheet(self):
        reference_sheet = self.dir + self.get('reference_sheet_file')
        from_here = os.path.dirname(reference_sheet) + os.sep

        with open(reference_sheet) as rf:
            csv_reader = csv.DictReader(rf, delimiter=',')
            for row in csv_reader:
                anno = row['Reference Annotation Files']
                mask = row['Reference Mask Annotation Files']
                anti = row['Also count antisense?'].lower()

                self.append_to('ref_annotations', self.cwl_file_def(from_here + anno))

                if mask.lower() in ('none', ''):
                    self.append_to('mask_annotations', self.cwl_file_def(self.extras + '_empty_maskfile_aquatx.gff'))
                else:
                    self.append_to('mask_annotations', self.cwl_file_def(from_here + mask))

                if anti in ('true', 'false'):
                    self.append_to('antisense', anti)
                else:
                    raise ValueError('The value associated with reference file %s for '
                                     'antisense counting is not true/false.'
                                     % row['Reference Annotation Files'])

    def setup(self):
        """Populates default values and prepares per-file configuration lists"""

        dt = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.set_default_dict({
            'run_directory': dt + self.get('user') + "_aquatx",
            'run_prefix':    dt + self.get('user') + "_aquatx",
            'output_prefix': self.get('run_prefix'),
            'run_date':      dt.split('_')[0],
            'run_time':      dt.split('_')[1]
        })

        self.extras = resource_filename('aquatx', 'extras/')
        self.set('output_file_stats', self.get('output_prefix')+ '_run_stats.csv')
        self.set('output_file_counts', self.get('output_prefix') + '_raw_counts.csv')

        # Per-file settings lists to be populated by process_sample_sheet() and process_reference_sheet()
        for settings_list in ['ref_annotations', 'mask_annotations', 'antisense', 'json', 'html', 'un',
                              'outfile', 'uniq_seq_file', 'report_title', 'out_fq', 'in_fq', 'out_prefix']:
            self.set(settings_list, [])

        # Bowtie index files
        bt_idx = (self.set('ebwt', self.prefix(self.get('ref_genome')))
                  if self.get('run_idx') and not self.get('ebwt')
                  else self.get('ebwt'))
        self.set('bt_index_files',
                 [self.cwl_file_def(bt_idx + fpath) for fpath in
                  ['.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt']])

        if self.get('adapter_sequence') == 'auto_detect':
            self.config.pop('adapter_sequence')

    """========== GETTERS AND SETTERS =========="""

    def get(self, key: str) -> Union[str,list,dict]:
        return self.config.get(key, None)

    def set(self, key: str, val: Union[str,list,dict]) -> Union[str,list,dict]:
        self.config[key] = val
        return val

    def set_default(self, key: str, val: str) -> str:
        """Apply the setting if it has not been previously set"""
        if key not in self.config:
            return self.set(key, val)
        else:
            return self.get(key)

    def set_default_dict(self, setting_dict: dict) -> None:
        """Apply all settings in the input dictionary if they have not been previously set"""
        for key,val in setting_dict.items():
            self.set_default(key,val)

    def append_to(self, key: str, val: Union[str,list,dict]) -> list:
        """Append a file setting to a per-file settings list"""
        target = self.get(key)
        if key:
            target.append(val)
            return target
        else: print("Tried appending to a non-existent key.", file=sys.stderr)

    """========== HELPERS =========="""

    def cwl_file_def(self, file: str) -> dict:
        """Returns a file input/output specification for the CWL config"""
        return {'class': 'File', 'path': file}

    def prefix(self, path: str) -> str:
        """Returns everything from path except the file extension"""
        return os.path.splitext(path)[0]