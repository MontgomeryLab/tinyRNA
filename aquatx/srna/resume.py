import os
import re
import sys

from ruamel.yaml.comments import CommentedOrderedMap
from pkg_resources import resource_filename
from abc import ABC, abstractmethod
from datetime import datetime
from glob import glob

from aquatx.srna.Configuration import ConfigBase, timestamp_format


class ResumeConfig(ConfigBase, ABC):
    """A class for modifying the workflow and config to resume at a given step

    The modified workflow document must be written to the package resource directory
    for CWL workflows in order to preserve the relative paths copied from the original
    workflows. Directory output names have a timestamp appended to them to keep outputs
    separate between runs. Timestamped run prefixes are also updated to reflect the
    date and time of the resume run.
    """

    def __init__(self, processed_config, workflow, steps, entry_inputs):
        # Parse the pre-processed YAML configuration file
        super().__init__(processed_config)
        self.check_dir_requirements(processed_config)

        # Load the CWL workflow YAML for modification
        with open(workflow, 'r') as f:
            self.workflow: CommentedOrderedMap
            self.workflow = self.yaml.load(f)

        self.dt = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.entry_inputs = entry_inputs
        self.steps = steps

        self._create_truncated_workflow()
        self._rebuild_entry_inputs()
        self._add_timestamps()

    def check_dir_requirements(self, config_file):
        """Ensure that user has invoked this command from within the target run directory"""

        target_run_dir = os.path.basename(self['run_directory'])
        config_file_dir = os.path.dirname(config_file)
        invocation_dir = os.path.basename(os.getcwd())

        if config_file_dir != '' or invocation_dir != target_run_dir:
            sys.exit("Resume commands must be executed from within the target Run Directory.")

    @abstractmethod
    def _rebuild_entry_inputs(self):
        """This method must be implemented by the inheriting class"""
        pass

    def _create_truncated_workflow(self):
        """Modifies the workflow to run only the specified steps

        New input variables must be created which point to the files produced by the
        prior pipeline run. The entry/ step is modified to receive these inputs directly
        rather than expecting CWL to run the previous workflow steps. Parameters between
        resumed steps don't require modification. The subdirs step is removed and replaced
        with modified copies of its own calls to make-subdir.cwl.
        """

        wf_steps = self.workflow["steps"]
        wf_inputs = self.workflow["inputs"]
        wf_outputs = self.workflow["outputs"]

        # Remove all steps except those defined for the resume run
        for step in list(wf_steps.keys()):
            if step not in self.steps:
                del wf_steps[step]

        # Remove all WorkflowOutputParameters in preparation for the new
        wf_outputs.clear()

        # Setup new inputs at the workflow level and entry step
        for param, new_input in self.entry_inputs.items():
            # Create WorkflowInputParameters (workflow level)
            wf_inputs[new_input['var']] = new_input['type']
            # Update WorkflowStepInputs (entry step)
            wf_steps[self.steps[0]]['in'][param] = new_input['var']

        # Load the organize-outputs subworkflow so that we may copy relevant steps
        with open(resource_filename('aquatx', 'cwl/workflows/organize-outputs.cwl')) as f:
            organizer_sub_wf = self.yaml.load(f)

        for step in self.steps:
            # Copy relevant steps from organize-outputs.cwl
            step_name = f'organize_{step}'
            context = wf_steps[step_name] = organizer_sub_wf['steps'][step_name]

            # Update WorkflowStepInputs
            context['in']['dir_name'] = f'dir_name_{step}'
            context['in']['dir_files'] = {
                'source': [f"{step}/{output}" for output in wf_steps[step]['out']],
                'pickValue': 'all_non_null',
                'linkMerge': 'merge_flattened'
            }

            # Update WorkflowOutputParameter
            wf_outputs[f'{step}_out_dir'] = {
                'type': "Directory",
                'outputSource': f'{step_name}/subdir'
            }

    def _add_timestamps(self):
        """Differentiates resume-run output subdirs by adding a timestamp to them"""

        # Rename output directories with timestamp
        for subdir in self.steps:
            step_dir = "dir_name_" + subdir
            self[step_dir] = self[step_dir] + "_" + self.dt

        # Update run_name output prefix variable for the current date and time
        self['run_name'] = re.sub(timestamp_format, self.dt, self['run_name'])

    def write_workflow(self, workflow_outfile: str) -> None:
        with open(workflow_outfile, "w") as wf:
            self.yaml.dump(self.workflow, wf)


class ResumeCounterConfig(ResumeConfig):
    """A class for modifying the workflow and config to resume a run at Counter"""

    def __init__(self, processed_config, workflow):
        steps = ["counter", "dge", "plotter"]

        inputs = {
            'aligned_seqs': {'var': "resume_sams", 'type': "File[]"},
            'fastp_logs': {'var': "resume_fastp_logs", 'type': "File[]"},
            'collapsed_fa': {'var': "resume_collapsed_fas", 'type': "File[]"},
        }

        # Parse the pre-processed YAML configuration file
        super().__init__(processed_config, workflow, steps, inputs)

    def _rebuild_entry_inputs(self):
        """Set the new path inputs for the Counter step

        Normally, the Counter step receives inputs from the current run's previous step outputs.
        When resuming, these outputs should already exist on disk. We need to populate the new
        File[] arrays with their corresponding pipeline outputs on disk.
        """

        def cwl_file_resume(subdir, file):
            try:
                return self.cwl_file('/'.join([subdir, file]))
            except FileNotFoundError as e:
                sys.exit("The following pipeline output could not be found:\n%s" % (e.filename,))

        self['resume_sams'] = [cwl_file_resume(self['dir_name_bowtie'], sam) for sam in self['outfile']]
        self['resume_fastp_logs'] = [cwl_file_resume(self['dir_name_fastp'], log) for log in self['json']]
        self['resume_collapsed_fas'] = [cwl_file_resume(self['dir_name_collapser'], prefix + "_collapsed.fa")
                                        for prefix in self['uniq_seq_prefix']]


class ResumePlotterConfig(ResumeConfig):
    """A class for modifying the workflow and config to resume a run at Plotter"""

    def __init__(self, processed_config, workflow):
        steps = ["plotter"]

        inputs = {
            'raw_counts': {'var': "resume_raw", 'type': "File"},
            'norm_counts': {'var': "resume_norm", 'type': "File"},
            'deg_tables': {'var': "resume_deg", 'type': "File[]"},
            'len_dist': {'var': "resume_len_dist", 'type': "File[]"}
        }

        # Parse the pre-processed YAML configuration file
        super().__init__(processed_config, workflow, steps, inputs)

    def _rebuild_entry_inputs(self):
        """Set the new path inputs for the Plotter step

        Normally, the Plotter step receives inputs from the current run's previous step outputs.
        When resuming, these outputs should already exist on disk. We need to populate the new
        File[] arrays with their corresponding pipeline outputs on disk.
        """

        counter = self['dir_name_counter']
        dge = self['dir_name_dge']

        try:
            self['resume_raw'] = self.cwl_file(glob(counter + "/*_feature_counts.csv")[0])
            self['resume_norm'] = self.cwl_file(glob(dge + "/*_norm_counts.csv")[0])
            self['resume_len_dist'] = list(map(self.cwl_file, glob(counter + "/*_nt_len_dist.csv")))
            self['resume_deg'] = list(map(self.cwl_file, glob(dge + "/*_deseq_table.csv")))
        except FileNotFoundError as e:
            sys.exit("The following pipeline output could not be found:\n%s" % (e.filename,))
