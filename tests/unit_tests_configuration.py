import contextlib
import csv
import io
import os
import random
import unittest
from unittest.mock import patch, mock_open, call

from tiny.rna.configuration import Configuration, SamplesSheet, PathsFile, get_templates
from unit_test_helpers import csv_factory, paths_template_file, make_paths_file
from tiny.rna.util import CSVReader, r_reserved_keywords


def patch_isfile():
    return patch('tiny.rna.configuration.os.path.isfile', return_value=True)


def patch_open(read_data):
    """Constructs a mock_open that supports rewinding read_data via seek(0)."""

    def rewind(offset):
        if offset == 0:
            # Same as calling mock_open's function-scoped reset_data()
            mopen.side_effect()
        else:
            # In order to support this we need access to _state in mock_open's function scope.
            # There's a way to do that, but it's ugly and the patch isn't persistent.
            raise ValueError("Patched mock_open.seek() only supports seek(0).")

    mopen = mock_open(read_data=read_data)
    mopen.return_value.seek.side_effect = rewind
    return patch('builtins.open', mopen)


class BowtieIndexesTests(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.root_cfg_dir = os.path.abspath("./testdata/config_files")
        self.run_config = self.root_cfg_dir + "/run_config_template.yml"
        self.paths = self.root_cfg_dir + "/paths.yml"

        self.maxDiff = 1522

    """============ Helper functions ============"""

    def config_with(self, prefs):
        config = Configuration(self.run_config)
        for key, val in prefs.items():
            config[key] = val
        return config

    def bt_idx_files_from_prefix(self, prefix):
        if not os.path.abspath(prefix):
            prefix = os.path.normpath(os.path.join(self.root_cfg_dir, prefix))

        return [
            {'path': f"{prefix}.{subext}.ebwt", 'class': 'File'}
            for subext in ['1', '2', '3', '4', 'rev.1', 'rev.2']
        ]

    def get_default_prefix(self, config):
        return os.path.join(
            self.root_cfg_dir,
            config['run_directory'],
            "bowtie-build/ram1"
        )

    """================ Tests =================="""

    """Does get_ebwt_prefix() produce the expected prefix path?"""

    def test_get_ebwt_prefix(self):
        config = Configuration(self.run_config)
        actual_prefix = config.get_ebwt_prefix()
        expected_prefix = self.get_default_prefix(config)

        self.assertEqual(actual_prefix, expected_prefix)

    """Does get_ebwt_prefix() throw an error if reference genome files aren't provided?"""

    def test_get_ebwt_prefix_no_genome(self):
        config = Configuration(self.run_config)
        config['reference_genome_files'] = None

        with self.assertRaises(ValueError):
            config.get_ebwt_prefix()

    """Does get_bt_index_files() output the paths of indexes that have already been built?"""

    def test_get_bt_index_files_prebuilt_indexes(self):
        config = self.config_with({'run_bowtie_build': False})
        prefix = os.path.abspath("./testdata/counter/validation/ebwt/ram1")
        config.paths['ebwt'] = prefix
        expected = self.bt_idx_files_from_prefix(prefix)
        self.assertListEqual(config.get_bt_index_files(), expected)

    """Does get_bt_index_files() output the paths of the index files that are expected
    to be built from the reference genome?"""

    def test_get_bt_index_files_unbuilt_indexes_with_genome(self):
        config = self.config_with({'run_bowtie_build': True})
        config.paths['ebwt'] = "mock_prefix"
        prefix = config.paths['ebwt']
        expected = self.bt_idx_files_from_prefix(prefix)
        self.assertListEqual(config.get_bt_index_files(), expected)

    """Does get_bt_index_files() produce an error and quit when index files are
    missing and a reference genome has not been provided?"""

    def test_get_bt_index_files_missing_indexes_without_genome(self):
        config = self.config_with({'run_bowtie_build': False, 'reference_genome_files': None})
        prefix = "missing"
        config.paths['ebwt'] = prefix
        errmsg = '\n'.join([
            "The following Bowtie index file couldn't be found:",
            "\t" + f"{self.root_cfg_dir}/{prefix}.1.ebwt",
            "\nPlease either correct your ebwt prefix or add reference genomes in the Paths File."
        ])

        with self.assertRaisesRegex(SystemExit, errmsg):
            config.get_bt_index_files()

    """Does get_bt_index_files() produce an error without quitting when index files
    are missing but a reference genome was provided, and does it return the list of
    index files that will be built from the genome?"""

    def test_get_bt_index_files_missing_indexes_with_genome(self):
        config = self.config_with({'run_bowtie_build': False})
        bad_prefix = "missing"
        config.paths['ebwt'] = bad_prefix
        genome_prefix = self.get_default_prefix(config)

        expected_files = self.bt_idx_files_from_prefix(genome_prefix)
        expected_error = '\n'.join([
            "The following Bowtie index file couldn't be found:",
            "\t" + f"{self.root_cfg_dir}/{bad_prefix}.1.ebwt",
            "\nIndexes will be built from your reference genome files during this run.",
            ""
        ])

        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            actual = config.get_bt_index_files()

        self.assertEqual(stderr.getvalue(), expected_error)
        self.assertListEqual(actual, expected_files)

    """Does verify_bowtie_build_outputs() update the paths in ["bt_index_files"] and rewrite 
    these changes to the processed Run Config if long indexes were produced? Does it also
    write to the Paths File to update the new ebwt prefix?"""

    def test_verify_bowtie_build_outputs(self):
        ebwt_short = ["1.ebwt", "2.ebwt", "3.ebwt"]
        ebwt_long = ["1.ebwtl", "2.ebwtl", "3.ebwtl"]
        run_conf_ebwt = [Configuration.cwl_file(f, verify=False) for f in ebwt_short]
        expected_ebwt = [Configuration.cwl_file(f, verify=False) for f in ebwt_long]

        config = self.config_with({'bt_index_files': run_conf_ebwt})

        with patch('tiny.rna.configuration.open', mock_open()) as mo, \
                patch('tiny.rna.configuration.glob', return_value=ebwt_long) as g:
            config.verify_bowtie_build_outputs()

        expected_writes = [
            call(self.paths, 'w'),
            call(os.path.join(self.root_cfg_dir, config['run_directory'], os.path.basename(self.run_config)), 'w')
        ]

        self.assertListEqual(config['bt_index_files'], expected_ebwt)
        self.assertListEqual(mo.call_args_list, expected_writes)


class SamplesSheetTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.contexts = ["Pipeline Start", "Pipeline Step", "Standalone Run"]

    """Does SamplesSheet catch multi-assignment of control condition?"""

    def test_validate_control_group(self):
        rows = [
            {'File': '1', 'Group': 'G1', 'Replicate': '1', 'Control': True, 'Normalization': ''},  # Good
            {'File': '2', 'Group': 'G1', 'Replicate': '2', 'Control': True, 'Normalization': ''},  # Good
            {'File': '3', 'Group': 'G2', 'Replicate': '1', 'Control': True, 'Normalization': ''}   # Bad
        ]                        # ^^^^                               ^^^^

        start_sheet = csv_factory("samples.csv", [dict(r, File=f"{i}.fastq") for i, r in enumerate(rows)])
        step_sheet =  csv_factory("samples.csv", [dict(r, File=f"{i}.sam")   for i, r in enumerate(rows)])
        run_sheet = step_sheet

        with patch_isfile():
            # Control condition should not be evaluated in Pipeline Step or Standalone Run context
            with patch_open(step_sheet):
                SamplesSheet('mock_filename', context="Pipeline Step")
            with patch_open(run_sheet):
                SamplesSheet('mock_filename', context="Standalone Run")

            # Control condition should be evaluated in Pipeline Start context
            exp_contains = r".*(multiple control conditions).*"
            with self.assertRaisesRegex(AssertionError, exp_contains), patch_open(start_sheet):
                SamplesSheet('mock_filename', context="Pipeline Start")

    """Does SamplesSheet catch duplicate entries for the same group and rep?"""

    def test_validate_group_rep(self):
        rows = [
            {'File': '1', 'Group': 'G1', 'Replicate': '1', 'Control': True, 'Normalization': ''},  # Good
            {'File': '2', 'Group': 'G1', 'Replicate': '2', 'Control': True, 'Normalization': ''},  # Good
            {'File': '3', 'Group': 'G1', 'Replicate': '2', 'Control': True, 'Normalization': ''}   # Bad
        ]                        # ^^^^               ^^^

        start_sheet = csv_factory("samples.csv", [dict(r, File=f"{i}.fastq") for i, r in enumerate(rows)])
        step_sheet =  csv_factory("samples.csv", [dict(r, File=f"{i}.sam")   for i, r in enumerate(rows)])
        run_sheet = step_sheet
        contexts = list(zip(self.contexts, [start_sheet, step_sheet, run_sheet]))

        exp_contains = r".*(same group and replicate).*"

        # Validation should take place in all contexts
        for context, sheet in contexts:
            with self.assertRaisesRegex(AssertionError, exp_contains), patch_open(sheet), patch_isfile():
                SamplesSheet('mock_filename', context=context)

    """Does SamplesSheet catch bad fastq entries in Pipeline Start context?"""

    def test_validate_fastq_filepath(self):
        context = "Pipeline Start"
        csv_rows = [
            {'File': '1.fastq', 'Group': 'G1', 'Replicate': '1', 'Control': True, 'Normalization': ''},  # Good
            {'File': '1.fastq', 'Group': 'G1', 'Replicate': '2', 'Control': True, 'Normalization': ''},  # Bad
            {'File': '2.fasta', 'Group': 'G2', 'Replicate': '1', 'Control': True, 'Normalization': ''}   # Bad
        ]           # ^^^^^^^
        sheet = csv_factory("samples.csv", csv_rows)

        # File doesn't exist
        exp_contains = r".*(fastq file on row 1 of mock_filename was not found).*"
        with self.assertRaisesRegex(AssertionError, exp_contains), patch_open(sheet):
            SamplesSheet('mock_filename', context=context)

        # Duplicate filename
        exp_contains = r".*(listed more than once).*\(row 2\)"
        with self.assertRaisesRegex(AssertionError, exp_contains), patch_open(sheet), patch_isfile():
            SamplesSheet('mock_filename', context=context)

        # Bad file extension
        sheet = csv_factory("samples.csv", csv_rows[1:])  # remove duplicate entry
        exp_contains = r".*(\.fastq\(\.gz\) extension).*\(row 2\)"
        with self.assertRaisesRegex(AssertionError, exp_contains), patch_open(sheet), patch_isfile():
            SamplesSheet('mock_filename', context=context)

    """Does SamplesSheet catch bad alignment file entries in Standalone Run context?"""

    def test_validate_alignments_filepath(self):
        context = "Standalone Run"
        csv_rows = [
            {'File': '1.sam', 'Group': 'G1', 'Replicate': '1', 'Control': True, 'Normalization': ''},  # Good
            {'File': '1.sam', 'Group': 'G1', 'Replicate': '2', 'Control': True, 'Normalization': ''},  # Good
            {'File': '1.bam', 'Group': 'G1', 'Replicate': '3', 'Control': True, 'Normalization': ''},  # Bad
            {'File': '2.xyz', 'Group': 'G2', 'Replicate': '1', 'Control': True, 'Normalization': ''}   # Bad
        ]          # ^^^^^^^
        sheet = csv_factory("samples.csv", csv_rows)

        # File doesn't exist
        exp_contains = r".*(file on row 1 of mock_filename was not found).*"
        with self.assertRaisesRegex(AssertionError, exp_contains), patch_open(sheet):
            SamplesSheet('mock_filename', context=context)

        # Duplicate filename
        exp_contains = r".*(listed more than once).*\(row 2\)"
        with self.assertRaisesRegex(AssertionError, exp_contains), patch_open(sheet), patch_isfile():
            SamplesSheet('mock_filename', context=context)

        # Bad file extension
        sheet = csv_factory("samples.csv", csv_rows[1:])  # remove duplicate entry
        exp_contains = r".*(\.sam or \.bam extension).*\(row 3\)"
        with self.assertRaisesRegex(AssertionError, exp_contains), patch_open(sheet), patch_isfile():
            SamplesSheet('mock_filename', context=context)

    """Does validate_r_safe_sample_groups() detect group names that will cause namespace collisions in R?"""

    def test_validate_r_safe_sample_groups(self):
        non_alphanum_chars = [bad.join(('a', 'b')) for bad in "~!@#$%^&*()+-=`<>?/,:;\"'[]{}\| \t\n\r\f\v"]
        leading_dot_number = [".0", "X.0"]

        for bad in [non_alphanum_chars, leading_dot_number]:
            msg = " ≈ ".join(bad)
            with self.assertRaisesRegex(AssertionError, msg):
                SamplesSheet.validate_r_safe_sample_groups(dict.fromkeys(bad))

        for kwd in r_reserved_keywords:
            bad = (kwd, kwd + '.')
            msg = " ≈ ".join(bad)
            with self.assertRaisesRegex(AssertionError, msg):
                SamplesSheet.validate_r_safe_sample_groups(dict.fromkeys(bad))

    """Does validate_normalization() do what it should?"""

    def test_validate_normalization(self):
        good = [
            {'File': '1.fastq', 'Group': 'G1', 'Replicate': '1', 'Control': False, 'Normalization': 'rpm'},      # Good
            {'File': '2.fastq', 'Group': 'G1', 'Replicate': '2', 'Control': False, 'Normalization': 'RPM'},      # Good
            {'File': '3.fastq', 'Group': 'G1', 'Replicate': '3', 'Control': False, 'Normalization': ' RPM '},    # Good
            {'File': '4.fasta', 'Group': 'G2', 'Replicate': '1', 'Control': False, 'Normalization': '  1'},      # Good
            {'File': '5.fasta', 'Group': 'G2', 'Replicate': '2', 'Control': False, 'Normalization': '1.1'},      # Good
            {'File': '6.fasta', 'Group': 'G2', 'Replicate': '3', 'Control': False, 'Normalization': ''},         # Good
        ]                                                                                          # ^^^^

        start_sheet = csv_factory("samples.csv", [dict(r, File=f"{i}.fastq") for i, r in enumerate(good)])
        step_sheet = csv_factory("samples.csv", [dict(r, File=f"{i}.sam") for i, r in enumerate(good)])
        run_sheet = step_sheet

        # These SHOULD NOT throw an error
        for context, sheet in zip(self.contexts, [start_sheet, step_sheet, run_sheet]):
            with patch_open(sheet), patch_isfile():
                SamplesSheet("mock_filename", context=context)

        bad = [
            {'File': '1.fastq', 'Group': 'G1', 'Replicate': '1', 'Control': False, 'Normalization': 'abc'},      # Bad
            {'File': '2.fastq', 'Group': 'G1', 'Replicate': '2', 'Control': False, 'Normalization': '123.rpm'},  # Bad
            {'File': '3.fastq', 'Group': 'G1', 'Replicate': '3', 'Control': False, 'Normalization': '.'},        # Bad
            {'File': '1.fastq', 'Group': 'G2', 'Replicate': '1', 'Control': False, 'Normalization': '_'},        # Bad
        ]

        start_sheet = csv_factory("samples.csv", [dict(r, File=f"{i}.fastq") for i, r in enumerate(bad)])
        step_sheet = csv_factory("samples.csv", [dict(r, File=f"{i}.sam") for i, r in enumerate(bad)])
        run_sheet = step_sheet

        exp_contains = r".*(Invalid normalization value).*"

        # These SHOULD throw an error
        for context, sheet in zip(self.contexts, [start_sheet, step_sheet, run_sheet]):
            with self.assertRaisesRegex(AssertionError, exp_contains), patch_open(sheet), patch_isfile():
                SamplesSheet("mock_filename", context=context)


class PathsFileTests(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.template_dir = os.path.dirname(paths_template_file)

    """Does PathsFile automatically resolve paths when queried?"""

    def test_getitem_single(self):
        config = make_paths_file()
        config['mock_parameter'] = "./some/../file"
        self.assertEqual(config['mock_parameter'], os.path.join(self.template_dir, "file"))

    """Does PathsFile automatically resolve lists of paths when queried? What if some list items are
    strings, others are mappings with a "path" key, and others are None/empty?"""

    def test_getitem_group(self):
        config = make_paths_file()
        config.groups = ('mock_parameter',)

        mapping_1 = {'path': "./some/../file", "other_key": "irrelevant"}
        mapping_2 = {'path': "../config_files/another_file"}
        path_string = "../../../START_HERE/reference_data/ram1.gff3"

        config['mock_parameter'] = [mapping_1, mapping_2, path_string, None, '']

        # Notice: only 'path' is modified in mappings and empty entries are returned unmodified
        expected = [
            dict(mapping_1, path=os.path.join(self.template_dir, "file")),
            dict(mapping_2, path=os.path.join(self.template_dir, "another_file")),
            os.path.normpath(os.path.join(self.template_dir, path_string)),
            None,
            ''
        ]

        self.assertListEqual(config['mock_parameter'], expected)

    """Does PathsFile check for required parameters?"""

    def test_validate_required_parameters(self):
        config = make_paths_file()
        for key in PathsFile.required:
            oldval = config[key]

            with self.assertRaisesRegex(AssertionError, r".*(parameters are required).*"):
                config[key] = ''
                config.validate_paths()

            with self.assertRaisesRegex(AssertionError, r".*(parameters are required).*"):
                config[key] = None
                config.validate_paths()

            with self.assertRaisesRegex(AssertionError, r".*(parameters are required).*"):
                del config.config[key]
                config.validate_paths()

            config[key] = oldval

    """Does PathsFile notify the user of the change in counting style 
    when no GFF files have been provided?"""

    def test_no_gff_files(self):
        config = make_paths_file()
        stdout = io.StringIO()

        gff_configs = [
            [{'path': "", 'alias': []}],
            [{'irrelevant': "value"}],
            None,
            []
        ]

        with contextlib.redirect_stdout(stdout):
            # Skip countdown timer for each config
            with patch('tiny.rna.configuration.time.sleep') as countdown:
                for gff_config in gff_configs:
                    config['gff_files'] = gff_config
                    config.validate_paths()
                    config_return = config.get_gff_config()

                    self.assertRegex(stdout.getvalue(), r".*(No GFF files).*")
                    self.assertEqual(config_return, {})
                    stdout.truncate(0)

        countdown.assert_has_calls([call(1)] * 6 * len(gff_configs))

    """Does PathsFile check for missing files for single entry parameters?"""

    def test_validate_missing_file_single(self):
        config = make_paths_file()
        key = random.choice(PathsFile.single)
        config[key] = '/dev/null/file_dne'

        with self.assertRaisesRegex(AssertionError, f".*(file provided for {key}).*"):
            config.validate_paths()

    """Does PathsFile check for missing files under list-type parameters?"""

    def test_validate_missing_file_group(self):
        config = make_paths_file()
        bad_path = "/dev/null/file_dne"
        config.append_to('gff_files', {'path': bad_path, 'alias': []})

        # Config now has one good entry and one bad entry; make sure the bad entry was specifically reported
        with self.assertRaisesRegex(AssertionError, f".*(file provided under gff_files).*\n\t{bad_path}"):
            config.validate_paths()

    """Does PathsFile detect a backward compatibility issue?"""

    def test_backward_incompatibility(self):
        config = make_paths_file()
        del config.config['gff_files']

        with self.assertRaisesRegex(AssertionError, r".*(check the release notes).*"):
            config.check_backward_compatibility()

    """Does PathsFile map all paths to the working directory when context=True?"""

    def test_pipeline_mapping(self):
        # There are entries in the template Paths File that will cause FileNotFound
        #  when context=True (they're not in the template directory). Since
        #  file existence isn't relevant for this test, we patch that out.
        with patch('tiny.rna.configuration.os.path.isfile', return_value=True):
            config = make_paths_file(in_pipeline=True)

        config['mock_mapping'] = {'path': "/dev/null/file_dne", 'other': "irrelevant"}
        config['mapping_no_path'] = {'nopath': True}
        config['mock_path'] = "/a/very/long/path.gz"
        config['empty_path'] = ""
        config['none_path'] = None

        self.assertDictEqual(config['mock_mapping'], {'path': "file_dne", 'other': "irrelevant"})
        self.assertDictEqual(config['mapping_no_path'], {'nopath': True})
        self.assertEqual(config['mock_path'], "path.gz")
        self.assertEqual(config['empty_path'], "")
        self.assertEqual(config['none_path'], None)

    """Does get_gff_config properly screen for "ID" alias attributes?"""

    def test_get_gff_config_id_alias_attr(self):
        config = make_paths_file()
        config['gff_files'] = [{'path': "/irrelevant", 'alias': ['ID']}]

        gff_files = config.get_gff_config()
        expected = {"/irrelevant": []}

        self.assertDictEqual(gff_files, expected)

    """Does get_gff_config merge aliases if the same GFF file is listed more than once? 
    Are duplicates also removed and original order preserved?"""

    def test_get_gff_config_merge_alias_attr(self):
        config = make_paths_file()
        config['gff_files'] = [
            {'path': "/same/path", 'alias': ['attr1', 'attr2']},
            {'path': "/same/path", 'alias': ['attr1', 'attrN']},
            {'path': "/different/path", 'alias': ['attr3']}
        ]

        gff_files = config.get_gff_config()

        expected = {
            "/same/path": ['attr1', 'attr2', 'attrN'],
            "/different/path": ['attr3']
        }

        self.assertDictEqual(gff_files, expected)


class ConfigurationTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.file = './testdata/config_files/run_config_template.yml'

    """Does get_templates copy the expected number of files for each context?"""

    def test_get_templates_contexts(self):
        context_file_count = {
            'tiny': 5,
            'tiny-count': 3,
            'tiny-plot': 1
        }

        with patch('tiny.rna.configuration.shutil.copyfile') as cf:
            for context, count in context_file_count.items():
                cf.reset_mock()
                get_templates(context)
                # Check against unique calls to copyfile incase there are duplicate calls for some reason
                self.assertEqual(len(set(call.args[0] for call in cf.call_args_list)), count)

    """Does get_templates properly handle cases where template files already exist in the CWD?"""

    def test_get_templates_conflicts(self):
        tiny_count = ['paths.yml', 'samples.csv', 'features.csv']
        tiny_plot = ['tinyrna-light.mplstyle']
        tiny = ['run_config_template.yml', *tiny_count, *tiny_plot]

        contexts = {
            'tiny': tiny,
            'tiny-count': tiny_count,
            'tiny-plot': tiny_plot
        }

        for context, files in contexts.items():
            with patch('tiny.rna.configuration.os.path.exists') as pe:
                pe.return_value = True

                try:
                    get_templates(context)
                except SystemExit as e:
                    err_msg = e.args[0]
                    err_files = err_msg.splitlines()[1:-1]
                    exp_set = set(files)
                    act_set = set(map(str.strip, err_files))

                    self.assertSetEqual(exp_set, act_set)

    """Does setup_pipeline() create the proper run_name prefix?"""

    def test_setup_pipeline_run_name(self):
        # user_vals = ('"user"', None, "''", 0)
        # run_vals =  ('"run"',  None, "''", 0)
        # for a, b in itertools.product(user_vals, run_vals):
        #   print("({'user': %6s, 'run_name': %5s}, )" % (a, b))

        cases = [
            ({'user': "user", 'run_name': "run"}, "run_{ts}"),
            ({'user': "user", 'run_name':  None}, "user_tinyrna_{ts}"),
            ({'user': "user", 'run_name':    ''}, "user_tinyrna_{ts}"),
            ({'user': "user", 'run_name':     0}, "0_{ts}"),
            ({'user':   None, 'run_name': "run"}, "run_{ts}"),
            ({'user':   None, 'run_name':  None}, "tinyrna_{ts}"),
            ({'user':   None, 'run_name':    ''}, "tinyrna_{ts}"),
            ({'user':   None, 'run_name':     0}, "0_{ts}"),
            ({'user':     '', 'run_name': "run"}, "run_{ts}"),
            ({'user':     '', 'run_name':  None}, "tinyrna_{ts}"),
            ({'user':     '', 'run_name':    ''}, "tinyrna_{ts}"),
            ({'user':     '', 'run_name':     0}, "0_{ts}"),
            ({'user':      0, 'run_name': "run"}, "run_{ts}"),
            ({'user':      0, 'run_name':  None}, "0_tinyrna_{ts}"),
            ({'user':      0, 'run_name':    ''}, "0_tinyrna_{ts}"),
            ({'user':      0, 'run_name':     0}, "0_{ts}"),
        ]

        for inputs, output in cases:
            config = Configuration(self.file, skip_setup=True)
            config.config.update(inputs)
            config.setup_pipeline()

            actual = config['run_name']
            expected = output.format(ts=config.dt)
            self.assertEqual(actual, expected)

    """Does setup_pipeline() create the proper run_directory suffix?"""

    def test_setup_pipeline_run_directory(self):
        run_name = "run"
        cases = [
            ({'run_directory':       "dir"}, "run_{ts}_dir"),
            ({'run_directory':  "path/dir"}, "path/run_{ts}_dir"),
            ({'run_directory': "path/dir/"}, "path/run_{ts}_dir"),
            ({'run_directory':         "/"}, "run_{ts}"),
            ({'run_directory':        None}, "run_{ts}"),
            ({'run_directory':          ''}, "run_{ts}"),
            ({'run_directory':           0}, "run_{ts}_0"),
        ]

        for inputs, output in cases:
            config = Configuration(self.file, skip_setup=True)
            config.config.update(inputs)
            config['run_name'] = run_name
            config.setup_pipeline()

            actual = config['run_directory']
            expected = output.format(ts=config.dt)
            self.assertEqual(actual, expected)

    """Does CSVReader handle CSVs with other delimiters?
    Adherence to RFC 4180 isn't strict in the wild. Delimiter can also be locale-dependent, 
    e.g. if a locale's decimal separator is a comma, then the delimiter is likely to be 
    converted to a semicolon when the user saves the file."""

    def test_csv_reader_delimiter(self):
        for doctype in CSVReader.tinyrna_sheet_fields.keys():
            fieldnames = CSVReader.tinyrna_sheet_fields[doctype].keys()

            for delimiter in csv.Sniffer().preferred:
                csv_body = io.StringIO()
                csv_writer = csv.DictWriter(csv_body, fieldnames=fieldnames, delimiter=delimiter)
                test_row = {col: '.' for col in fieldnames}

                csv_writer.writeheader()
                csv_writer.writerow(test_row)

                with patch_open(read_data=csv_body.getvalue()):
                    csv_reader = CSVReader("/dev/null", doctype)
                    read = list(csv_reader.rows())
                    self.assertEqual(list(read[0].values()), list(test_row.values()))

    """Does CSVReader raise an error if the CSV has a header but no body?"""

    def test_csv_no_body(self):
        for doctype in CSVReader.tinyrna_sheet_fields.keys():
            fieldnames = CSVReader.tinyrna_sheet_fields[doctype].keys()

            csv_body = io.StringIO()
            csv_writer = csv.DictWriter(csv_body, fieldnames=fieldnames, delimiter=';')
            csv_writer.writeheader()

            with self.assertRaisesRegex(csv.Error, rf"{doctype} is empty"):
                with patch_open(read_data=csv_body.getvalue()):
                    csv_reader = CSVReader("/dev/null", doctype)
                    next(csv_reader.rows())

    """Does CSVReader raise an error if the CSV has a row with an unexpected column count?"""

    def test_csv_column_count_consistency(self):
        for doctype in CSVReader.tinyrna_sheet_fields.keys():
            fieldnames = list(CSVReader.tinyrna_sheet_fields[doctype].keys())
            test_rows = [
                ['.' for _ in fieldnames[:-1]],         # Missing last column
                ['.' for _ in [*fieldnames, "Extra"]]   # Extra column
            ]

            for test_row in test_rows:
                csv_body = io.StringIO()
                csv_writer = csv.writer(csv_body)

                csv_writer.writerow(fieldnames)  # Header
                csv_writer.writerow(test_row)

                with self.assertRaisesRegex(csv.Error, rf"Inconsistent column count"):
                    with patch_open(read_data=csv_body.getvalue()):
                        csv_reader = CSVReader("/dev/null", doctype)
                        list(csv_reader.rows())


if __name__ == '__main__':
    unittest.main()
