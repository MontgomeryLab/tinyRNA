"""
This script uses the simple configuration file and spreadsheets to create the
configuration YAML that can be used as input in the CWL workflow.
"""
import os
import csv
import argparse
import ruamel.yaml
from shutil import copyfile
from datetime import datetime
from pkg_resources import resource_filename

"""
Remaining issues observed in the processed (output) config file:
    - ebwt relative path is dropped from processed config file BUT
    it is preserved under the bt_index_files option
"""

def get_args():
    """
    Get the input arguments
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input-file', metavar='CONFIG', required=True,
                        help="Input file")

    args = parser.parse_args()

    return args

def process_sample_sheet(sample_file, config_settings):
    """
    Function to process the input sample sheet.
    """
    sample_files = list()
    sample_identifiers = list()
    sample_out_fq = list()
    sample_report = list()
    sample_uniq_seq = list()
    sample_sams = list()
    sample_unaligned = list()
    sample_json = list()
    sample_html = list()

    sample_sheet_dir = os.path.dirname(sample_file) + os.sep
    collapser_extension = '.fa.gz' if config_settings['compress'] else '.fa'

    with open(sample_file) as sf:
        csv_reader = csv.DictReader(sf, delimiter=',')
        for row in csv_reader:
            sample_basename = os.path.splitext(os.path.basename(row['\ufeffInput FastQ/A Files']))[0]
            sample_files.append({'class': 'File',
                                 'path': sample_sheet_dir + row['\ufeffInput FastQ/A Files']})
            sample_identifiers.append(row['Sample/Group Name'] + '_replicate_'
                                      + row['Replicate number'])
            sample_out_fq.append(sample_basename + '_cleaned.fastq')
            sample_report.append(row['Sample/Group Name'] + '_replicate_'
                                 + row['Replicate number'] + '_fastp_report')
            sample_uniq_seq.append(sample_basename + '_collapsed' + collapser_extension)
            sample_sams.append(sample_basename + '_aligned_seqs.sam')
            sample_unaligned.append(sample_basename + '_unaligned_seqs.fa')
            sample_json.append(sample_basename + '_qc.json')
            sample_html.append(sample_basename + '_qc.html')

            # Todo: find out if this was necessary for plotter... ultimately it is never added to config...
            # sample_compare.add(row['Sample/Group Name'])

    config_settings['json'] = sample_json
    config_settings['html'] = sample_html
    config_settings['un'] = sample_unaligned
    config_settings['outfile'] = sample_sams
    config_settings['uniq_seq_file'] = sample_uniq_seq
    config_settings['report_title'] = sample_report
    config_settings['out_fq'] = sample_out_fq
    config_settings['in_fq'] = sample_files
    config_settings['out_prefix'] = sample_identifiers

    return config_settings

def process_reference_sheet(ref_file, config_settings):
    """
    Function to process the input reference sheet.
    """
    ref_files = list()
    mask_files = list()
    antisense = list()

    aquatx_extras_path = resource_filename('aquatx', 'extras/')
    ref_sheet_dir = os.path.dirname(ref_file) + os.sep

    with open(ref_file) as rf:
        csv_reader = csv.DictReader(rf, delimiter=',')
        for row in csv_reader:
            if config_settings['run_idx'] and config_settings['ebwt'] == '':
                bt_idx = os.path.splitext(config_settings['ref_genome'])[0]
                config_settings['ebwt'] = bt_idx
            else:
                bt_idx = config_settings['ebwt']
            config_settings['bt_index_files'] = list(({'class': 'File',
                                                        'path': bt_idx + '.1.ebwt'},
                                                       {'class': 'File',
                                                        'path': bt_idx + '.2.ebwt'},
                                                       {'class': 'File',
                                                        'path': bt_idx + '.3.ebwt'},
                                                       {'class': 'File',
                                                        'path': bt_idx + '.4.ebwt'},
                                                       {'class': 'File',
                                                        'path': bt_idx + '.rev.1.ebwt'},
                                                       {'class': 'File',
                                                        'path': bt_idx + '.rev.2.ebwt'}))

            ref_files.append({'class': 'File', 'path': ref_sheet_dir + row['Reference Annotation Files']})
            if row['Reference Mask Annotation Files'].lower() in ('none', ''):
                mask_files.append({'class': 'File', 'path': aquatx_extras_path + '_empty_maskfile_aquatx.gff'})
            else:
                mask_files.append({'class': 'File', 'path': ref_sheet_dir + row['Reference Mask Annotation Files']})

            if row['Also count antisense?'].lower() in ('true', 'false'):
                antisense.append(row['Also count antisense?'].lower())
            else:
                raise ValueError('The value associated with reference file %s for '
                                 'antisense counting is not true/false.'
                                 % row['Reference Annotation Files'])

    config_settings['ebwt'] = os.path.splitext(os.path.basename(bt_idx))[0]
    config_settings['ref_annotations'] = ref_files
    config_settings['mask_annotations'] = mask_files
    config_settings['antisense'] = antisense

    return config_settings

def setup_config(input_file):
    """
    Function to set up the configuration file from template
    """

    # Parse YAML run configuration file
    with open(input_file) as conf:
        config_settings = ruamel.yaml.YAML().load(conf)

    # Run time information
    time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    config_settings['run_date'] = time.split('_')[0]
    config_settings['run_time'] = time.split('_')[1]

    if config_settings['run_directory'] == '':
        config_settings['run_directory'] = time + config_settings['user'] + '_aquatx'
    if config_settings['run_prefix'] == '':
        config_settings['run_prefix'] = time + config_settings['user'] + '_aquatx'
    if config_settings['output_prefix'] == '':
        config_settings['output_prefix'] = config_settings['run_prefix']

    # Path specification should be relative to the config file, not relative to the script's CWD
    sample_sheet_path = os.path.dirname(input_file) + os.sep + config_settings['sample_sheet_file']
    config_settings = process_sample_sheet(sample_sheet_path, config_settings)

    config_settings['output_file_stats'] = config_settings['output_prefix'] + '_run_stats.csv'
    config_settings['output_file_counts'] = config_settings['output_prefix'] + '_raw_counts.csv'

    # Path specification should be relative to the config file, not relative to the script's CWD
    reference_sheet_path = os.path.dirname(input_file) + os.sep + config_settings['reference_sheet_file']
    config_settings = process_reference_sheet(reference_sheet_path, config_settings)

    if config_settings.get('adapter_sequence', None) == 'auto_detect':
        config_settings.pop('adapter_sequence')

    # TODO
    if os.path.basename(input_file) == 'run_config_template.yml':
        input_name = time + '_run_config.yml'
        copyfile(input_file, input_name)
    else:
        input_name = input_file

    with open(input_name, 'w') as outconf:
        ruamel.yaml.round_trip_dump(config_settings, outconf,
                                    default_flow_style=False, indent=4, block_seq_indent=2)
    
    # output filenames to stdout for running cwltool
    return input_name



def main():
    """
    Main routine to process the run information.
    """

    # Get input config file
    args = get_args()
    setup_config(args.input_file)

    #TODO: need to specify the non-model organism run when no reference genome is given



if __name__ == '__main__':
    main()
