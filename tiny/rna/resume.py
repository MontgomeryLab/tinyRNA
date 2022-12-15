import os
import re
import sys

import ruamel.yaml
from ruamel.yaml.comments import CommentedOrderedMap
from pkg_resources import resource_filename
from abc import ABC, abstractmethod
from glob import glob

from tiny.rna.configuration import ConfigBase, PathsFile
from tiny.rna.util import timestamp_format, get_timestamp


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
        self.check_dir_requirements()

        # Load the workflow YAML for modification
        if '~' in workflow: os.path.expanduser(workflow)
        with open(workflow, 'r') as f:
            self.workflow: CommentedOrderedMap
            self.workflow = self.yaml.load(f)

        self.dt = get_timestamp()
        self.entry_inputs = entry_inputs
        self.steps = steps + [f"organize_{s}" for s in steps]

        self.paths = self.load_paths_config()
        self.assimilate_paths_file()
        self.setup_step_inputs()

        self._create_truncated_workflow()
        self._rebuild_entry_inputs()
        self._add_timestamps(steps)

    def check_dir_requirements(self):
        """Ensure that user has invoked this command from within the target run directory"""

        target_run_dir = os.path.basename(self['run_directory'])
        invocation_dir = os.path.basename(os.getcwd())

        if invocation_dir != target_run_dir:
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

        # Remove all unassociated WorkflowOutputParameters
        for step in list(wf_outputs.keys()):
            if step[:step.rindex("_out_dir")] not in self.steps:
                del wf_outputs[step]

        # Setup new inputs at the workflow level and entry step
        for param, new_input in self.entry_inputs.items():
            # Create WorkflowInputParameters (workflow level)
            wf_inputs[new_input['var']] = new_input['type']
            # Update WorkflowStepInputs (entry step)
            wf_steps[self.steps[0]]['in'][param] = new_input['var']

    def load_paths_config(self):
        """Constructs a sub-configuration object containing workflow file preferences"""
        paths_file_path = self.from_here(self['paths_config'])
        return PathsFile(paths_file_path)

    def assimilate_paths_file(self):
        """Updates the processed workflow with resume-safe Paths File parameters"""
        for key in [*PathsFile.single, *PathsFile.groups]:
            if key in PathsFile.resume_forbidden: continue
            self[key] = self.paths.as_cwl_file_obj(key)
        for key in PathsFile.prefix:
            if key in PathsFile.resume_forbidden: continue
            self[key] = self.paths[key]

    def _add_timestamps(self, steps):
        """Differentiates resume-run output subdirs by adding a timestamp to them"""

        # Rename output directories with timestamp
        for subdir in steps:
            step_dir = "dir_name_" + subdir
            self[step_dir] = self[step_dir] + "_" + self.dt

        # Update run_name output prefix variable for the current date and time
        self['run_name'] = re.sub(timestamp_format, self.dt, self['run_name'])

    # Override
    def get_outfile_path(self, infile: str = None) -> str:
        if infile is None: infile = self.inf
        root, ext = os.path.splitext(os.path.basename(infile))
        return '_'.join(["resume", root, self.dt]) + ext

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

        # Build resume config from the previously-processed Run Config
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

        resume_file_lists = ['resume_sams', 'resume_fastp_logs', 'resume_collapsed_fas']
        self.set_default_dict({key: [] for key in resume_file_lists})

        for sample in self['sample_basenames']:
            self['resume_sams'].append(cwl_file_resume(self['dir_name_bowtie'], sample + '_aligned_seqs.sam'))
            self['resume_fastp_logs'].append(cwl_file_resume(self['dir_name_fastp'], sample + '_qc.json'))
            self['resume_collapsed_fas'].append(cwl_file_resume(self['dir_name_collapser'], sample + '_collapsed.fa'))


class ResumePlotterConfig(ResumeConfig):
    """A class for modifying the workflow and config to resume a run at Plotter"""

    def __init__(self, processed_config, workflow):
        steps = ["plotter"]

        inputs = {
            'raw_counts': {'var': "resume_raw", 'type': "File"},
            'rule_counts': {'var': "resume_rule", 'type': "File"},
            'summ_stats': {'var': "resume_stat", 'type': "File"},
            'len_dist_tables': {'var': "resume_len_dist", 'type': "File[]"}
        }

        dge_outputs = {
            'norm_counts': {'var': "resume_norm", 'type': "File"},
            'dge_tables': {'var': "resume_dge", 'type': "File[]"}
        }

        if self._dge_subdir_exists(processed_config):
            inputs.update(dge_outputs)

        # Build resume config from the previously-processed Run Config
        super().__init__(processed_config, workflow, steps, inputs)

    def _dge_subdir_exists(self, processed_config):
        with open(processed_config, 'r') as f:
            config = ruamel.yaml.YAML().load(f)
            exists = os.path.isdir(config['dir_name_dge'])

        self.dge_ran = exists
        return exists

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
            self['resume_stat'] = self.cwl_file(glob(counter + "/*_summary_stats.csv")[0])
            self['resume_rule'] = self.cwl_file(glob(counter + "/*_counts_by_rule.csv")[0])
            self['resume_len_dist'] = list(map(self.cwl_file, glob(counter + "/*_nt_len_dist.csv")))

            if self.dge_ran:
                self['resume_dge'] = list(map(self.cwl_file, glob(dge + "/*_deseq_table.csv")))
                self['resume_norm'] = self.cwl_file(glob(dge + "/*_norm_counts.csv")[0])

        except FileNotFoundError as e:
            sys.exit("The following pipeline output could not be found:\n%s" % (e.filename,))
        except IndexError:
            msg = "Expected outputs for feature counts or norm counts could not be found. "
            if counter is None or not os.path.isdir(counter): msg += f"The directory {counter} was not found. "
            if dge is None or not os.path.isfile(dge): msg += f"The directory {dge} was not found. "
            sys.exit(msg)

        steps = self.workflow['steps']
        plotter_inputs = steps['plotter']['in']

        # Remove the PCA plot input since dge is not a step in this resume workflow
        steps['organize_plotter']['in']['dir_files']['source'].remove('dge/pca_plot')
        if not self.dge_ran:
            prune = set()
            for input_var, value in plotter_inputs.items():
                if type(value) is str and value.startswith('dge/'):
                    prune.add(input_var)
            for input_var in prune:
                del plotter_inputs[input_var]
